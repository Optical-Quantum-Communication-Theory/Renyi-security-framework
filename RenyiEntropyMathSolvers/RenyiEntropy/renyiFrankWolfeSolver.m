function [relEntLowerBound,modParser] = renyiFrankWolfeSolver(params,options,debugInfo)
% A 2 step Frank-Wolfe solver for the Renyi sandwich conditional entropy
% restricted to QKD (Using the Renyi Petz relative entropy).
%
% * Finite size
% * Variable length keys via QES
% * Formulated as an optimization over Eve's attack channel.
% 
% 
% This is the simplest of the Renyi entropy solvers and a good place to
% start learning their numerical representation. The conditional entropy of
% the key register conditioned on public announcements and Eve is computed
% by taking the dual to the Renyi Petz Conditional entropy of the key
% register conditioned on public announcements and Bob's measurement
% results. Simplifications are then performed to get the final expression
% as seen in https://arxiv.org/abs/2504.12248.
%
% Input parameters:
% * krausOps: A cell array of matrices. The cell array contains the Kraus
%   operators that form the G map on Alice and Bob's joint system. These
%   should form a completely positive trace non-increasing linear map. Each
%   Kraus operator must be the same size.
% * keyProj: A cell array of projection operators that extract the key from
%   G(\rho). These projection operators should sum to <= identity.
% * probTest: Scalar double representing the probability of sending a test
%   state.
% * renyiAlpha: The Renyi entropy parameter. A scalar double in the range
%   (1,2].
% * dimAPrime: Positive integer representing the dimension of the signal
%   state Alice sends to Bob.
% * stateTest: State from the source replacement scheme reduced to the
%   systems A, Alice's and A' forwarded to Bob, and conditioned on a test
%   round. The systems are ordered AA'.
% * stateGen: State from the source replacement scheme reduced to the
%   systems A, Alice's and A' forwarded to Bob, and conditioned on a
%   generation round.The systems are ordered AA'.
% * testConstraints: An array of class EqualityConstraint. Represents the
%   constraints from test rounds. For each con in testConstraints,
%   con.scalar represents the accepted frequency conditioned on test
%   rounds, and con.operator is the corresponding operator of Alice and
%   Bob's joint POVM element conditioned on test rounds. The sum of all
%   operators should be identity and the sum of all scalars should be 1.
% * blockDimsAPrime: If the blockDiagonal option is true, this is a list
%   that holds the numbers of dimensions for each block of the signal
%   system Alice sends.
% * blockDimsB: If the blockDiagonal option is true, this a list that holds
%   the numbers of dimensions of each block of Bob's system. For example,
%   if Bob is a qubit with an additional loss dimension, blockDimsB = [2,1]
%
% Output:
% * relEntLowerBound: The lower bound on the Renyi sandwich down arrow
%   relative entropy between the secret key register and Eve, given the
%   constraints.
%
% Options:
% * cvxSolver (global option): See makeGlobalOptionsParser for details.
% * cvxPrecision (global option): See makeGlobalOptionsParser for details.
% * verboseLevel (global option): See makeGlobalOptionsParser for details.
% * maxIter (20): maximum number of Frank Wolfe iteration steps taken to
%   minimize the relative entropy.
% * maxGap (1e-6): Exit condition for the Frank Wolfe algorithm. When the
%   relative gap between the current and previous iteration is small
%   enough, the Frank Wolfe algorithm exits and returns the current point.
%   The gap must be a positive scalar.
% * linearSearchPrecision (1e-20): Precision the fminbnd tries to achieve
%   when searching along the line between the current point and the points
%   along the gradient line. See fminbnd and optimset for more details.
% * linearSearchMinStep (0): Minimum step size fminbnd must take during the
%   Frank Wolf algorithm. Initially, this can help with faster convergence,
%   but can also just prevent convergence. See fminbnd for more details
%   (the second argument, x1, in the function). This is only used for the
%   "vanilla" Frank Wolfe method.
% * tolVert (1e-12): Tolerance for determining if two steps in the pairwise
%   Frank Wolfe algorithm should be considered equal. Must be a positive
%   scalar. This is only used for the "pairwise" Frank Wolfe method.
% * linearConstraintTolerance (1e-10): constraint tolerance on the
%   equalityConstraints, inequalityConstraints, vectorOneNormConstraints
%   and matrixOneNormConstraints. A bit broader than the name suggests.
% * frankWolfeMethod ("vanilla"): Choice between using the vanilla or
%   pairwise Frank Wolfe Algorithms. Select "vanilla" or "pairwise"
%   respectively.
% * blockDiagonal (false): Tells the solver whether to set up rhoAB to be
%   block diagonal. If true, the solver also requires two parameters
%   blockDimsA and blockDimsB, which tell the solver what the block
%   dimensions of A and B are respectively.
% * stopNegGap (false): During the direction finding step, all gaps should
%   be positive, but if the gradient is close to zero then numerical
%   inaccuracy may cause the subproblem routine to give an impossible
%   negative gap. If stopNegGap is true, the Frank Wolfe routine will exit
%   early. Otherwise, it will try and continue.
% * numericPerturbation (1e-8): Linear perturbation of the subnormalized
%   state after the G-map. Removes eigenvalues of zero and improves the
%   stability of the algorithm at the cost of some key rate. Will not work
%   for the single case, G(rho) = 0.
% * numericPerturbationAffine (0): Affine perturbation of the subnormalized
%   state after the G-map and linear perturbation numericPertrubation is
%   applied. Unlike numericPerturbation, this affine perturbation works for
%   G(rho) = 0, but the relative perturbation and key rate cost increases
%   for smaller subnormalized states. We recommend setting this to a very
%   low value and tweaking numericPertubation first.
% * QESMode (false): Toggles the Quantum Entropic Score Frank Wolfe routine
%   on and off. When False, the solver stops after a lowerbound for the
%   Renyi entropy is found. When True, The solver continues to determine
%   the kappa QES constant used for variable length security proofs.
%
% DebugInfo:
% * relEntStep1: Value of the relative entropy achieved during the
%   Frank-Wolfe approximate minimization. It can help with determining if
%   Frank-Wolfe could not converge, or just barely converged.
% * relEntLowerBound: Lower bound on the relative entropy from step 2 of
%   the solver. The lower bound returned by the solver takes the maximum
%   between this and 0.
% * pointsQueue: Stores the sequence of points from each Frank Wolfe step
%   in a SimpleQueue.
% * QES: Debug information returned from running the QES Frank Wolfe
%   routine.
% * qesConstStep1: Value of the QES constant kappa returned from step 1 of
%   Frank Wolfe routine. The Value is an upper bound on the true minimum.
% * qesCell: A 2 element cell array containing the gradient of the QES
%   function at the optimal point, and the lower bound on the QES constant
%   kappa for variable length security proofs.
% * subproblemStatus: Array of strings containing the status of the CVX
%   solver as it find the direction to move along for iterations of step 1.
%
% See also: QKDMathSolverModule, makeGlobalOptionsParser, FrankWolfe,
% optimset, SimpleQueue
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
optionsParser.addOptionalParam("blockDiagonal", false, @(x) mustBeMember(x, [true, false]));
optionsParser.addOptionalParam("numericPerturbation",1e-8,@(x) mustBeInRange(x,0,1))
optionsParser.addOptionalParam("numericPerturbationAffine",0,@(x) mustBeInRange(x,0,1))
optionsParser.addOptionalParam("QESMode", false, @(x) mustBeMember(x, [true, false]));
optionsParser.addOptionalParam("stopNegGap",false,@(x) mustBeMember(x, [true, false]));
optionsParser.parse(options);

options = optionsParser.Results;

%% module parser

modParser = moduleParser(mfilename);

% Kraus operators for G map and key projection operators for Z map
modParser.addRequiredParam("krausOps",@mustBeNonempty);
modParser.addAdditionalConstraint(@isCPTNIKrausOps,"krausOps");
modParser.addRequiredParam("keyProj",@mustBeNonempty);
modParser.addAdditionalConstraint(@mustBeAKeyProj,"keyProj");
modParser.addRequiredParam("probTest",@(x) x>=0);
modParser.addRequiredParam("renyiAlpha",@(x)mustBeInRange(x,1,2,"exclude-lower"));
modParser.addRequiredParam("dimAPrime",@mustBeInteger);
modParser.addAdditionalConstraint(@mustBePositive,"dimAPrime");

modParser.addRequiredParam("stateTest",@isDensityOperator);
modParser.addRequiredParam("stateGen",@isDensityOperator);


%same dimension for matrix multiplication
modParser.addAdditionalConstraint(@(x,y)size(x{1},1) == size(y{1},2),["krausOps","keyProj"])

% Add constraints
modParser.addRequiredParam("testConstraints",@(x)mustBeA(x,"EqualityConstraint")); % must sum to identity and 1.

% make sure the constraints all act on a rhoAB of the same size.
modParser.addAdditionalConstraint(@(testConstraints,krausOps)...
    all([testConstraints.rhoDim] == size(krausOps{1},2),"all"),...
    ["testConstraints","krausOps"]);

% block diagonal constraints (if enabled)
if options.blockDiagonal
    modParser.addRequiredParam("blockDimsAPrime", @isBlockDimsWellFormated);
    modParser.addRequiredParam("blockDimsB", @isBlockDimsWellFormated);

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
dimA = size(params.stateGen,1)/dimAPrime;
dimB = size(params.krausOps{1},2)/dimA;

testCons = params.testConstraints;
probTest = params.probTest;

%version with alphaHat
% alphaHat = 1/(2-params.renyiAlpha);
% entropyAlpha = alphaHat;

%version with alpha
entropyAlpha = params.renyiAlpha;

%% start the solver

%generate maximally mixed state, later used as initial guess for FW iteration
choi0 = eye(dimAPrime*dimB)/dimB;

% generate block diagonal transformation matrices
if options.blockDiagonal
    [options.blockP, options.newDims] = blockRearrange(params.blockDimsAPrime, params.blockDimsB);
end

% curry the FW functions
perturbation = options.numericPerturbation;
perturbationAff = options.numericPerturbationAffine;

func     = @(vecFW)     RenyiTildeDownPA.funcFW(vecFW,params.stateGen,...
    params.keyProj,params.krausOps,probTest,entropyAlpha,dimAPrime,...
    perturbation,perturbationAff);
gradFunc = @(vecFW) RenyiTildeDownPA.gradFuncFW(vecFW,params.stateGen,...
    params.keyProj,params.krausOps,probTest,entropyAlpha,dimAPrime,...
    perturbation,perturbationAff);

subProblemFunc = @(vecFW,vecFWGrad) RenyiTildeDownPA.subproblemUnique(vecFW,vecFWGrad,...
    params.stateTest,dimAPrime,probTest,testCons,options);

% run step 1 routines
[vecFWInit,exitFlag] = RenyiTildeDownPA.closestVecFW(choi0,params.stateTest,...
    probTest,testCons,dimA,dimAPrime,dimB,options);
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
[deltaVec,exitFlag] = subProblemFunc(vecFW,gradFunc(vecFW));
if exitFlag == SubProblemExitFlag.failed
    error("Could not compute gap for step 2.")
end
gap = -FrankWolfe.inProdR(deltaVec,gradFunc(vecFW));
relEntLowerBound = func(vecFW) - abs(gap);

% using the best lowerbound point from FW (currently we trust our func and
% gradfunc enough to use the value output by FW directly)
if options.verboseLevel >= 1 && output.lowerBoundFWVal > relEntLowerBound
    fprintf("!!A better point was found!!\nold %e, new %e\n",relEntLowerBound,output.lowerBoundFWVal)
end
relEntLowerBound = max(relEntLowerBound,output.lowerBoundFWVal);

% %% testing the new step 2
% % An experimental step 2 solver that uses the intersect of multiple
% % linearizations for a better bound. We can revisit this idea at a later
% % date.
%
% vecFWs = debugInfo.info.pointsQueue.toCell;%{vecFW;output.lowerBoundFWPoint};
% vecGrads = cellfun(@(x) gradFunc(x), vecFWs,"UniformOutput",false);
% constOffsets = cellfun(@(x,y) func(x)-FrankWolfe.inProdR(x,y),vecFWs,vecGrads);
% [testLowerBound,exitFlag] = RenyiTildeDownPA.simpleStep2(constOffsets,vecGrads,params.stateTest,dimAPrime,probTest,testCons,options);
% %disp(exitFlag)
% 
% if testLowerBound > relEntLowerBound
%     fprintf("!!A better point was found using new step 2!!\nold %e, new %e\n",relEntLowerBound,testLowerBound)
% end
% relEntLowerBound = max(relEntLowerBound,testLowerBound);

%Store best lower bound
debugInfo.storeInfo("relEntLowerBound",relEntLowerBound);

%% Find QES for variable length
if options.QESMode
    
    %Optimal linearization point
    optimalLinPoint = output.lowerBoundFWPoint;

    [qesGrad,tolGrad,exitFlag] = QESFinder.getGradQES(optimalLinPoint, ...
        gradFunc(optimalLinPoint),params.stateTest,dimAPrime,probTest, ...
        testCons,options);
    if exitFlag == SubProblemExitFlag.failed
        error("Could not get gradient for QES.")
    end

    
    funcQes     = @(vecQesFW)     QESFinder.funcFW(vecQesFW,params.stateGen,...
        params.keyProj,params.krausOps,testCons,probTest,qesGrad,entropyAlpha,...
        dimAPrime,perturbation,perturbationAff);
    gradFuncQes = @(vecQesFW) QESFinder.gradFuncFW(vecQesFW,params.stateGen,...
        params.keyProj,params.krausOps,testCons,probTest,qesGrad,entropyAlpha,...
        dimAPrime,perturbation,perturbationAff);
    subProblemFuncQes = @(vecQesFW,vecQesGrad) QESFinder.subproblemQes(vecQesFW,...
        vecQesGrad,params.stateTest,dimAPrime,testCons,options);
    
    %Initial point for FW iteration    
    vecQesInit = QESFinder.initPoint(optimalLinPoint(1:end-1),params.stateTest,...
        testCons,dimAPrime);
    
    debugQES = debugInfo.addLeaves("QES");
    switch options.frankWolfeMethod
        case "vanilla"
            [vecQes,qesConst,exitFlagQes,outputQes] = FrankWolfe.vanilla(vecQesInit, funcQes,...
                gradFuncQes, subProblemFuncQes, debugQES, options.verboseLevel>=1,...
                "linearSearchMinStep",options.linearSearchMinStep, ...
                "linearSearchPrecision",options.linearSearchPrecision,...
                "maxGap",options.maxGap,...
                "maxIter",options.maxIter,...
                "storePoints",false);
        case "pairwise"
            [vecQes,qesConst,exitFlagQes,outputQes] = FrankWolfe.pairwise(vecQesInit, funcQes,...
                gradFuncQes, subProblemFuncQes, debugQES, options.verboseLevel>=1,...
                "tolVert",options.tolVert,...
                "linearSearchPrecision",options.linearSearchPrecision,...
                "maxGap",options.maxGap,...
                "maxIter",options.maxIter,...
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
    
    %step 2 for QES
    [deltaVecQes,exitFlagStep2Qes,extraOut,tolKappa] = subProblemFuncQes(vecQes,gradFuncQes(vecQes));
    if exitFlagStep2Qes == SubProblemExitFlag.failed
        error("Could not compute gap for QES step 2.")
    end
    gapQes = -FrankWolfe.inProdR(deltaVecQes,gradFuncQes(vecQes));
    qesConstLowerBound = funcQes(vecQes) - abs(gapQes);
    % qesConstLowerBound = -1/(alphaHat-1)*log2(-(funcQes(vecQes) - abs(gapQes)));
    
    %Pick best lowerbound
    if options.verboseLevel >= 1 && outputQes.lowerBoundFWVal > qesConstLowerBound
        fprintf("!!A better QES constant was found!!\nold %e, new %e\n",qesConstLowerBound,outputQes.lowerBoundFWVal);
    end
    qesConstLowerBound = max(qesConstLowerBound,outputQes.lowerBoundFWVal);

    %Store QES as a cell array
    qesCell = {qesGrad,qesConstLowerBound};
    debugInfo.storeInfo("qesCell",qesCell);
end

end

