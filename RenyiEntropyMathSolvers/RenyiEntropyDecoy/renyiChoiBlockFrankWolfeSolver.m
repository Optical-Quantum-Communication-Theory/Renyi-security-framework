function [relEntLowerBound,modParser] = renyiChoiBlockFrankWolfeSolver(params,options,debugInfo)
% FW2StepSolver A 2 step solver to calculate the minimum relative entropy
% between the key and Eve, given a set of linear equality and inequality
% constraints on Alice and Bob's density matrix. Step 1 uses the Frank
% Wolfe algorithm to find the minimum by iteratively solving linearizations
% of the problem. Step 2 of the solver uses the linearization at solution
% of step 1 (the final iteration) and solves the dual at this point. While
% also taking into account for constraint violations, due to numerical
% imprecission, step 2 produces a valid lower bound on the relative
% entropy. See "Reliable numerical key rates for quantum key distribution",
% https://quantum-journal.org/papers/q-2018-07-26-77/, for more details.
%
% Input parameters:
% * krausOps: A cell array of matrices. The cell array contains the Kraus
%   operators that form the G map on Alice and Bob's joint system. These
%   should form a completely positive trace non-increasing linear map. Each
%   Kraus operator must be the same size.
% * keyProj: A cell array of projection operators that extract the key from
%   G(\rho). These projection operators should sum to <= identity.
% * equalityConstraints (EqualityConstraint.empty(0,1)): Array of class 
%   EqualityConstraint. Represents constraints of the form Tr[operator*rho]
%   = scalar. Up to the linearConstraintTolerance.
% * inequalityConstraints (InequalityConstraint.empty(0,1)): Same idea as
%   the equalityConstraints but uses the class InequalityConstraint.
%   Represents constraints of the form lowerBound <= Tr[operator*rho] <=
%   upperBound.
% * vectorOneNormConstraints (VectorOneNormConstraint.empty(0,1)): Array of
%   class VectorOneNormConstraint. Represents constraints of the form
%   ||sum_i Tr[operators_i rho] e_i - vector||_1 <= scalar which are often
%   used in the finite-size analysis.
% * matrixOneNormConstraints (MatrixOneNormConstraint.empty(0,1)): Array of
%   class MatrixOneNormConstraints. Represents constraints of the form
%   ||Phi(rho) - operator||_1 <= scalar, where Phi is a hermitian
%   preserving superoperator and operator is a hermitian operator and the
%   norm is the trace norm.
% * blockDimsA: If the blockDiagonal option is true, this is a list that
%   holds the numbers of dimensions of each block of Alice's system.
% * blockDimsB: If the blockDiagonal option is true, this a list that holds
%   the numbers of dimensions of each block of Bob's system. For example,
%   if Bob is a qubit with an additional loss dimension, blockDimsB = [2,1]
% * rhoA (nan): The density matrix of Alice, if known to be unchanged (for
%   example, in the source replacement scheme of prepare-and-measure
%   schemes). Optional, as entanglement based protocols do not use rhoA,
%   but including rhoA when possible can significantly improve key rate.
%   This produces a set of equality constraints based on the information
%   from rhoA. If your proof technique has a more nuanced relation, then
%   you may have to remove this and add it by hand to your other
%   constraints.
% Output:
% * relEntLowerBound: The lower bound on the relative entropy between the
%   key and Eve, given the constraints.
% Options:
% * cvxSolver (global option): See makeGlobalOptionsParser for details.
% * cvxPrecision (global option): See makeGlobalOptionsParser for details.
% * verboseLevel (global option): See makeGlobalOptionsParser for details.
% * maxIter (20): maximum number of Frank Wolfe iteration steps taken to
%   minimize the relative entropy.
% * maxGap (1e-6): Exit condition for the Frank Wolfe algithm. When the
%   relative gap between the current and previous iteration is small
%   enough, the Frank wolfe algorithm exits and returns the current point.
%   The gap must be a postive scalar.
% * linearSearchPrecision (1e-20): Precision the fminbnd tries to achieve
%   when searching along the line between the current point and the points
%   along the gradient line. See fminbnd and optimset for more details.
% * linearSearchMinStep (1e-3): Minimum step size fminbnd must take
%   during the Frank Wolf algorithm. Initially, this can help with faster
%   convergence, but can also just prevent convergence. See fminbnd for
%   more details (the second argument, x1, in the function). This is only
%   used for the "vanilla" Frank Wolfe method.
% * tolVert (1e-12): Tolerance for determining if two steps in the pairwise
%   Frank Wolfe agorithm should be concidered equal. Must be a positive
%   scalar. This is only used for the "pairwise" Frank Wolfe method.
% * linearConstraintTolerance (1e-10): constraint tolerance on the
%   equalityConstraints, inequalityConstraints, vectorOneNormConstraints
%   and matrixOneNormConstraints. A bit broader than the name suggests.
% * frankWolfeMethod ("vanilla"): Choice between using the vanilla or
%   pairwise Frank Wolfe Algorithms. Select "vanilla" or "pairwise"
%   respectively.
% * initMethod (1): Integer selected from {1,2,3}. For the Frank Wolfe
%   algorithm, the initial point must satisfy the constraints. This selects
%   which technique is used, from 3 choices:
%    #1 Closest point in the set to the maximally mixed state.
%    #2 Point in the set with the largest minimum eigenvalue.
%    #3 No objective function, just let CVX find a feasible point.
% * blockDiagonal (false): Tells the solver whether to set up rhoAB to be
%   block diagonal. If true, the solver also requires two parameters
%   blockDimsA and blockDimsB, which tell the solver what the block
%   dimensions of A and B are respectively.

% DebugInfo:
% * relEntStep1: Value of the relative entropy achieved during the
%   Frank-Wolfe aproximate minimization. It can help with determining if
%   Frank-Wolfe could not converge, or just barely converges.
% * relEntStep2Linearization: The relative entropy at the start of the
%   linearization routine used for step 2. The initial point may be
%   perturbed slightly from step 1.
% * relEntLowerBound: Lower bound on the relative entropy from step 2 of
%   the solver. The lower bound returned by the solver takes the maximum
%   between this and 0.
% * closestDensityMatrixStatus: String containing the status of the CVX
%   solver after attempting to solve for the closest density matrix in step
%   1 of the solver.
% * subproblemStatus: Array of strings containing the status of the CVX
%   solver as it find the direction to move along for iterations of step 1.
% * submaxproblemStatus: CVX status for the solver when calculating the
%   dual solution (a maximization) used in step 2's linearization.
% * dualSolution: Optimal value achieved from the subMaxProblem mentioned
%   above.
%
% See also: QKDMathSolverModule, makeGlobalOptionsParser, fminbnd, optimset
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options

%start with the global parser and add on the extra options
optionsParser = makeGlobalOptionsParser(mfilename);

optionsParser.addOptionalParam("maxIter",20,@mustBeInteger);
optionsParser.addAdditionalConstraint(@(x) x>0, "maxIter");
optionsParser.addOptionalParam("maxGap",1e-6,@(x) x>0);
optionsParser.addOptionalParam("linearSearchPrecision",1e-20,@(x) x>0);
optionsParser.addOptionalParam("linearSearchMinStep",0,@(x) 0<=x && x<=1);
optionsParser.addOptionalParam("tolVert",1e-12,@(x)x>0);
optionsParser.addOptionalParam("linearConstraintTolerance",1e-10,@(x) x>0);
optionsParser.addOptionalParam("frankWolfeMethod",  "vanilla",@(x)mustBeMember(x,["vanilla","pairwise"]));
optionsParser.addOptionalParam("initMethod",1,@(x) mustBeMember(x,[1,2,3]));
optionsParser.addOptionalParam("blockDiagonal", false, @(x) mustBeMember(x, [true, false]));
optionsParser.addOptionalParam("numericPerturbation",1e-8,@(x) mustBeInRange(x,0,1))
optionsParser.addOptionalParam("numericPerturbationAffine",0,@(x) mustBeInRange(x,0,1))
optionsParser.addOptionalParam("QESMode", false, @(x) mustBeMember(x, [true, false]));
optionsParser.addOptionalParam("stopNegGap",true,@(x) mustBeMember(x, [true, false]));
optionsParser.parse(options);

options = optionsParser.Results;

%% module parser

modParser = moduleParser(mfilename);

% kraus operators for G map and key projection operators for Z map
modParser.addRequiredParam("krausOps",@mustBeNonempty);
modParser.addAdditionalConstraint(@(x) allCells(x, @isCPTNIKrausOps),"krausOps");
modParser.addRequiredParam("keyProj",@mustBeNonempty);
modParser.addAdditionalConstraint(@(x) allCellsMust(x, @mustBeAKeyProj),"keyProj");
modParser.addRequiredParam("probTest",@(x) x>=0);
modParser.addRequiredParam("renyiAlpha",@(x)mustBeInRange(x,1,2,"exclude-lower"));
modParser.addRequiredParam("dimAPrime",@mustBeInteger);
modParser.addAdditionalConstraint(@mustBePositive,"dimAPrime");

modParser.addRequiredParam("stateTest",@(x) allCells(x,@isDensityOperator));
modParser.addRequiredParam("stateGen",@(x) allCells(x,@isDensityOperator));

modParser.addRequiredParam("probSignalConTestAndBlock",@(x) eachRowMustBeAProbDist(x.'));

modParser.addRequiredParam("probBlockCondTest",@(x) sum(x,2)<=1); %p(n|test) vector $n \in \{0,1\}$.
modParser.addRequiredParam("probBlockCondGen",@(x) sum(x,2)<=1);  %p(n|gen) vector $n \in \{0,1\}$.

modParser.addRequiredParam("blockPhotonNum", @(x) mustBeInteger(x));
modParser.addAdditionalConstraint(@(x) mustBeNonnegative(x),"blockPhotonNum");

modParser.addRequiredParam("probSignalConTest",@(x)mustBeProbDist(x)); % <--- check against statics size

%corrections due to phase imperfections
modParser.addRequiredParam("continuityBnd",@mustBeNonnegative);
modParser.addRequiredParam("epsilonBlockCondTest",@mustBeNonnegative);
modParser.addRequiredParam("epsilonProbTest",@mustBeNonnegative);


%same dimension for matrix multiplication
% modParser.addAdditionalConstraint(@(x,y) size(x{1},1) == size(y{1},2),["krausOps","keyProj"])

% Add constraints from test rounds
modParser.addRequiredParam("testConstraints",@(x)mustBeA(x,"EqualityConstraintChoiBlock")); % must sum to identity and 1. a|test x b|test x mu|test.

% Add prob of not being in blocks for each observation
modParser.addRequiredParam("probRemaining",@mustBeNonnegative);
modParser.addAdditionalConstraint(@(x) mustBeInRange(x,0,1),"probRemaining");


% Add constraints from squashing
modParser.addOptionalParam("squashingConstraintsGen",...
    InequalityConstraint.empty(0,1),@(x)mustBeA(x,"InequalityConstraint"));
modParser.addOptionalParam("squashingConstraintsTest",...
    InequalityConstraint.empty(0,1),@(x)mustBeA(x,"InequalityConstraint"));

% make sure the constraints all act on a rhoAB of the same size.
% modParser.addAdditionalConstraint(@(testConstraints,krausOps)...
%     all([testConstraints.rhoDim] == size(krausOps{1},2),"all"),...
%     ["testConstraints","krausOps"]);

% block diagonal constraints (if enabled)
if options.blockDiagonal
    modParser.addRequiredParam("blockDimsAPrime", @(x) allCells(x, @isBlockDimsWellFormated));
    modParser.addOptionalParam("blockDimsB", nan, @isBlockDimsWellFormated);

    % make sure they sum to the total number of dimensions % fix this for
    % choi matrix
    % modParser.addAdditionalConstraint(@(blockDimsAPrime,blockDimsB,krausOps)...
    %     sum(blockDimsAPrime,"all")*sum(blockDimsB,"all") == size(krausOps{1},2),...
    %     ["blockDimsAPrime","blockDimsB","krausOps"]);
end

%Apply modParser to input parameters
modParser.parse(params,"warnUnusedParams",true);
params = modParser.Results;

% get all the dimensions
dimAPrime = params.dimAPrime;
dimAB = cellfun(@(x) size(x{1},2), params.krausOps);
dimAAPrime = cellfun(@(x) size(x,1), params.stateTest);
dimA = dimAAPrime./dimAPrime;
dimB = dimAB./dimA;
dimAPrimeB = dimAPrime.*dimB;

%Test and squashing constraints (split by test and generation rounds)
testCons = params.testConstraints;
squashingConsTest = params.squashingConstraintsTest;
squashingConsGen = params.squashingConstraintsGen;

%Testing probability
probTest = params.probTest;

%version with alphaHat
% alphaHat = 1/(2-params.renyiAlpha);
% entropyAlpha = alphaHat;

%version with alpha
entropyAlpha = params.renyiAlpha;


%% start the solver

%generate maximally mixed state, later used as initial guess for FW iteration
% choi0 = eye(dimAPrime*dimB)/dimB;

% generate block diagonal transformation matrices
if options.blockDiagonal
    numBlocks = numel(params.blockDimsAPrime);
    options.blockP = cell(1,numBlocks);
    options.newDims = cell(1,numBlocks);
    for index = 1:numel(params.blockDimsAPrime)
        [options.blockP{index}, options.newDims{index}] = blockRearrange(params.blockDimsAPrime{index}, params.blockDimsB);
    end
end

% curry the FW functions
perturbation = options.numericPerturbation;
perturbationAff = options.numericPerturbationAffine;

func =     @(vecFW) RenyiTildeDownChoiBlockPA.funcFW(vecFW,probTest,entropyAlpha,perturbation,perturbationAff, ...
    dimAPrime,params.stateGen,params.keyProj,params.krausOps,params.probBlockCondGen);
gradFunc = @(vecFW) RenyiTildeDownChoiBlockPA.gradFuncFW(vecFW,probTest,entropyAlpha,perturbation,perturbationAff, ...
    dimAPrime,params.stateGen,params.keyProj,params.krausOps,params.probBlockCondGen);

subProblemFunc = @(vecFW,vecFWGrad) RenyiTildeDownChoiBlockPA.subproblemUniqueDecoy(vecFW,vecFWGrad,...
    params.stateTest,params.stateGen,dimAPrime,params.probTest,testCons,squashingConsTest,squashingConsGen,...
    params.probBlockCondTest,params.probSignalConTestAndBlock,params.probRemaining,...
    params.epsilonProbTest,params.epsilonBlockCondTest,options);

% run step 1 routines
numInit = dimAPrimeB*dimAPrimeB.' + 1;
vec0 = zeros(numInit,1);
grad0 = zeros(numInit,1);

[vecFWInit,exitFlag] = RenyiTildeDownChoiBlockPA.subproblemUniqueDecoy(vec0,grad0,...
    params.stateTest,params.stateGen,dimAPrime,params.probTest,testCons,squashingConsTest,squashingConsGen,...
    params.probBlockCondTest,params.probSignalConTestAndBlock,params.probRemaining,...
    params.epsilonProbTest,params.epsilonBlockCondTest,options,true);
if exitFlag == SubProblemExitFlag.failed
    error("Could not find an initial Point for Frank Wolfe solver.")
end

switch options.frankWolfeMethod
    case "vanilla"
        [vecFW,relEntStep1,exitFlag,output] = FrankWolfe.vanilla(vecFWInit, func,...
            gradFunc, subProblemFunc, debugInfo, options.verboseLevel>=1,...
            "linearSearchMinStep",options.linearSearchMinStep, ...
            "linearSearchPrecision",options.linearSearchPrecision,...
            "maxGap",options.maxGap,...
            "maxIter",options.maxIter,...
            "stopNegGap",options.stopNegGap,...
            "storePoints",true);
    case "pairwise"
        [vecFW,relEntStep1,exitFlag,output] = FrankWolfe.pairwise(vecFWInit, func,...
            gradFunc, subProblemFunc, debugInfo, options.verboseLevel>=1,...
            "tolVert",options.tolVert,...
            "linearSearchPrecision",options.linearSearchPrecision,...
            "maxGap",options.maxGap,...
            "maxIter",options.maxIter,...
            "stopNegGap",options.stopNegGap,...
            "storePoints",true);
end

% Warn the user if FW had any problems.
if options.verboseLevel >= 1
    switch exitFlag
        case FrankWolfeExitFlag.exceededMaxIter
            warning("FW2StepSolver:ExceededMaxIter", ...
                "Exceeded maximum number of Frank-Wolfe iterations.")
        case FrankWolfeExitFlag.subproblemFailed
            warning("FW2StepSolver:SubproblemFailed",...
                "Subproblem function could not determine step direction. " + ...
                "Using point from last iteration.")
    end
end

debugInfo.storeInfo("relEntStep1",relEntStep1);

%% step 2

% using the final point from FW
[deltaVec,exitFlagStep2FixedLength,extraOut] = subProblemFunc(vecFW,gradFunc(vecFW));
if exitFlagStep2FixedLength == SubProblemExitFlag.failed
    error("Could not find an initial Point for Frank Wolfe solver.")
end
gap = -FrankWolfe.inProdR(deltaVec,gradFunc(vecFW));
relEntLowerBound = func(vecFW) - abs(gap);

% using the best lowerbound point from FW (currently we trust our func and
% gradfunc enough to use the value output by FW directly)
if output.lowerBoundFWVal > relEntLowerBound
    fprintf("!!A better point was found!!\nold %e, new %e\n",relEntLowerBound,output.lowerBoundFWVal)
end
relEntLowerBound = max(relEntLowerBound,output.lowerBoundFWVal);

debugInfo.storeInfo("relEntLowerBound",relEntLowerBound);

if options.QESMode == true

    %% Step 1 for QES

    %optimal linearization point from fixed length calculation
    optimalLinPoint = output.lowerBoundFWPoint;
    optimalLinPoint = QESFinderChoiBlock.perturbEigvalsChoi(optimalLinPoint,dimAPrimeB,false);
    
    % find gradient of QES at linearization point
    [qesGrad,tolGrad,exitFlag] = QESFinderChoiBlock.getGradQES(optimalLinPoint,gradFunc(optimalLinPoint), ...
        params.stateTest,params.stateGen, dimAPrime, probTest, testCons,...
        squashingConsTest, squashingConsGen, params.probBlockCondTest, ...
        params.probSignalConTestAndBlock, params.probRemaining, ...
        params.epsilonProbTest,params.epsilonBlockCondTest, options);
    if exitFlag == SubProblemExitFlag.failed
        error("Could not get gradient for QES.")
    end

    %for phase imperfections continuity bnd needs to be added to the gen component (last entry) of qes gradient
    qesGrad(end) = qesGrad(end) + params.continuityBnd;

    %Define functions for finding constant term in Qes (required for FW
    %algorithm)
    funcQes = @(vecQesFW) QESFinderChoiBlock.funcFW(vecQesFW,params.stateGen,params.keyProj,params.krausOps,testCons, ...
        probTest,qesGrad,entropyAlpha,dimAPrime,params.probBlockCondGen,perturbation,perturbationAff);

    gradFuncQes = @(vecQesFW) QESFinderChoiBlock.gradFuncFW(vecQesFW,params.stateGen,params.keyProj,params.krausOps,testCons, ...
        probTest,qesGrad,entropyAlpha,dimAPrime,params.probBlockCondGen,perturbation,perturbationAff);

    subProblemFuncQes = @(vecQesFW,vecQesGrad) QESFinderChoiBlock.subproblemQes(vecQesFW,vecQesGrad,...
        params.stateTest,params.stateGen,dimAPrime,testCons,squashingConsTest,squashingConsGen, ...
        params.probBlockCondTest,params.probSignalConTestAndBlock,params.probRemaining,...
        params.epsilonProbTest,params.epsilonBlockCondTest,options);
    
    %Initial point for FW iteration
    [vecQesInit,exitFlag] = QESFinderChoiBlock.initPoint(optimalLinPoint(1:end-1),params.stateTest,params.stateGen,dimAPrime,...
        testCons,squashingConsTest,squashingConsGen,params.probBlockCondTest,params.probSignalConTestAndBlock,...
        params.probRemaining,params.epsilonProbTest,params.epsilonBlockCondTest,options);
    if exitFlag == SubProblemExitFlag.failed
        error("Could not find an initial Point for QES Frank Wolfe solver.")
    end
            
    debugQES = debugInfo.addLeaves("QES");
    %Run FW solver
    switch options.frankWolfeMethod
        case "vanilla"
            [vecQes,qesConst,exitFlagQes,outputQes] = FrankWolfe.vanilla(vecQesInit, funcQes,...
                gradFuncQes, subProblemFuncQes, debugQES, options.verboseLevel>=1,...
                "linearSearchMinStep",options.linearSearchMinStep, ...
                "linearSearchPrecision",options.linearSearchPrecision,...
                "maxGap",options.maxGap,...
                "maxIter",options.maxIter,...
                "stopNegGap",options.stopNegGap,...
                "storePoints",false);
        case "pairwise"
            [vecQes,qesConst,exitFlagQes,outputQes] = FrankWolfe.pairwise(vecQesInit, funcQes,...
                gradFuncQes, subProblemFuncQes, debugQES, options.verboseLevel>=1,...
                "tolVert",options.tolVert,...
                "linearSearchPrecision",options.linearSearchPrecision,...
                "maxGap",options.maxGap,...
                "maxIter",options.maxIter,...
                "stopNegGap",options.stopNegGap,...
                "storePoints",false);
    end
    debugInfo.storeInfo("qesConstStep1",qesConst);

    % Warn the user if FW had any problems.
    if options.verboseLevel >= 1
        switch exitFlagQes
            case FrankWolfeExitFlag.exceededMaxIter
                warning("FW2StepSolver:ExceededMaxIter", ...
                    "Exceeded maximum number of Frank-Wolfe iterations in QES.")
            case FrankWolfeExitFlag.subproblemFailed
                warning("FW2StepSolver:SubproblemFailed",...
                    "Subproblem function could not determine step direction in QES. " + ...
                    "Using point from last iteration.")
        end
    end

    %% step 2 for QES    
    %perturb result from step 1 if necessary
    vecQes = QESFinderChoiBlock.perturbEigvalsChoi(vecQes,dimAPrimeB,true);
    
    % Run one more FW to find reliable lower bound
    [deltaVecQes,exitFlagStep2Qes,extraOut] = subProblemFuncQes(vecQes,gradFuncQes(vecQes));
    if exitFlagStep2Qes == SubProblemExitFlag.failed
        error("Could not compute gap for QES step 2.")
    end
    gapQes = -FrankWolfe.inProdR(deltaVecQes,gradFuncQes(vecQes));
    qesConstLowerBound = funcQes(vecQes) - abs(gapQes);
    % qesConstLowerBound = -1/(alphaHat-1)*log2(-(funcQes(vecQes) - abs(gapQes)));
    
    %Pick best lowerbound
    if outputQes.lowerBoundFWVal > qesConstLowerBound
        fprintf("!!A better QES constant was found!!\nold %e, new %e\n",qesConstLowerBound,outputQes.lowerBoundFWVal);
    end
    qesConstLowerBound = max(qesConstLowerBound,outputQes.lowerBoundFWVal);

    %Store QES as a cell array
    qesCell = {qesGrad,qesConstLowerBound};
    debugInfo.storeInfo("qesCell",qesCell);
end
end

%% validation functions
function rowsMustSumTo1(conExp)
if~all(ismembertol(full(sum(conExp,2)),1))
    throwAsCaller(MException("deocyAnalysisIndependentLP:RowDontSumTo1",...
        "All rows must sum to 1."))
end
end

function eachRowMustBeAProbDist(expectationsConditional)

    % get the dimensions of the conditional expectations. Then based on that
    % pick a strategy to handle it
    dimExpCon = size(expectationsConditional);
    
    errorID ="BasicBB84WCPDecoyKeyRateFunc:InvalidRowsAreNotProbDists";
    errorTXT = "A row in the conditional distribution is not a valid probability distribution.";
    
    if numel(dimExpCon) == 2 % Matlab's minimum number of dimensions is 2.
        % The array is 2d and the slicing is easy
        for index = 1:dimExpCon(1)
            if~isProbDist(expectationsConditional(index,:))
               throwAsCaller(MException(errorID,errorTXT));
            end
        end
    else
        % We have some tricky slicing to do for 3 plus dimensions.
        % We need to index the first dimension and the combination of
        % dimensions 3 and up. The second dimension will just use :.
        maskedDims = [dimExpCon(1),prod(dimExpCon(3:end))];
    
        for index = 1:prod(maskedDims)
            vecIndex = ind2subPlus(maskedDims,index);
            if ~isProbDist(expectationsConditional(vecIndex(1),:,vecIndex(2)))
                throwAsCaller(MException(errorID,errorTXT));
            end
        end
    end
end