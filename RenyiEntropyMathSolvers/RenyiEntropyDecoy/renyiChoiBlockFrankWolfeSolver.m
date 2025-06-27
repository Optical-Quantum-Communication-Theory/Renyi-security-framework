function [relEntLowerBound,modParser] = renyiChoiBlockFrankWolfeSolver(params,options,debugInfo)

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
if options.verboseLevel >= 1 && output.lowerBoundFWVal > relEntLowerBound
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
    if options.verboseLevel >= 1 && outputQes.lowerBoundFWVal > qesConstLowerBound
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