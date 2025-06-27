function [relEntLowerBound,modParser] = renyiDecoyIntImpFrankWolfeSolver(params,options,debugInfo)

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

modParser.addRequiredParam("photonDistributionGenLower",@(x) sum(x,2)<=1); %lower bound p(n|gen) vector $n \in \{0,1\}$.
modParser.addRequiredParam("photonDistributionGenUpper",@(x) sum(x,2)<=1); %upper bound p(n|gen) vector $n \in \{0,1\}$.

modParser.addRequiredParam("probDistPhotonConMuTestLower", @(x) all(sum(x,1)<=1));
modParser.addRequiredParam("probDistPhotonConMuTestUpper", @(x) mustBeInRange(x,0,1));

modParser.addRequiredParam("probRemaining", @mustBeNonnegative);
modParser.addAdditionalConstraint(@(x) x<=1 ,"probRemaining");

modParser.addRequiredParam("blockPhotonNum", @(x) mustBeInteger(x));
modParser.addAdditionalConstraint(@(x) mustBeNonnegative(x),"blockPhotonNum");

modParser.addRequiredParam("decoyTest",@(x)allCells(x,@(y) y>=0)); %cell array vector
modParser.addRequiredParam("probDecoyConTest",@(x)mustBeProbDist(x)); % vector
modParser.addAdditionalConstraint(@(x,y) mustBeEqualSize(x,y), ["decoyTest","probDecoyConTest"]);
modParser.addRequiredParam("probSignalConTest",@(x)mustBeProbDist(x)); % <--- check against statics size


%same dimension for matrix multiplication
% modParser.addAdditionalConstraint(@(x,y)size(x{1},1) == size(y{1},2),["krausOps","keyProj"])

% Add constraints from test rounds
modParser.addRequiredParam("testConstraints",@(x)mustBeA(x,"EqualityConstraintDecoy")); % must sum to identity and 1. a|test x b|test x mu|test.

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


%SubProblem func equal for both cases
subProblemFunc = @(vecFW,vecFWGrad) RenyiTildeDownDecoyPA.subproblemUniqueDecoyIntImperfect(vecFW,vecFWGrad,...
    params.stateTest,params.stateGen,dimAPrime,probTest,testCons,squashingConsTest,squashingConsGen,...
    params.probDistPhotonConMuTestLower, params.probDistPhotonConMuTestUpper, params.probRemaining, ...
    params.probDecoyConTest,params.probSignalConTest,params.blockPhotonNum,options);


%initial point
numInit = dimAPrimeB*dimAPrimeB.' + 1;
vec0 = zeros(numInit,1);
grad0 = zeros(numInit,1);

[vecFWInit, exitFlag] = RenyiTildeDownDecoyPA.subproblemUniqueDecoyIntImperfect(vec0,grad0,...
    params.stateTest,params.stateGen,dimAPrime,probTest,testCons,squashingConsTest,squashingConsGen,...
    params.probDistPhotonConMuTestLower, params.probDistPhotonConMuTestUpper, params.probRemaining,...
    params.probDecoyConTest,params.probSignalConTest,params.blockPhotonNum,options,true);
if exitFlag == SubProblemExitFlag.failed
    error("Could not find an initial Point for Frank Wolfe solver.")
end


%% Run FW with lower bound on probability in generation rounds
funcLower =     @(vecFW) RenyiTildeDownDecoyPA.funcFW(vecFW,probTest,entropyAlpha,perturbation,perturbationAff, ...
    dimAPrime,params.stateGen,params.keyProj,params.krausOps,params.photonDistributionGenLower);
gradFuncLower = @(vecFW) RenyiTildeDownDecoyPA.gradFuncFW(vecFW,probTest,entropyAlpha,perturbation,perturbationAff, ...
    dimAPrime,params.stateGen,params.keyProj,params.krausOps,params.photonDistributionGenLower);


% run step 1 routines
switch options.frankWolfeMethod
    case "vanilla"
        [vecFWLower,relEntStep1Lower,exitFlag,outputLower] = FrankWolfe.vanilla(vecFWInit, funcLower,...
            gradFuncLower, subProblemFunc, debugInfo, options.verboseLevel>=1,...
            "linearSearchMinStep",options.linearSearchMinStep, ...
            "linearSearchPrecision",options.linearSearchPrecision,...
            "maxGap",options.maxGap,...
            "maxIter",options.maxIter,...
            "stopNegGap",options.stopNegGap,...
            "storePoints",true);
    case "pairwise"
        [vecFWLower,relEntStep1Lower,exitFlag,outputLower] = FrankWolfe.pairwise(vecFWInit, funcLower,...
            gradFuncLower, subProblemFunc, debugInfo, options.verboseLevel>=1,...
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

debugInfo.storeInfo("relEntStep1Lower",relEntStep1Lower);

% step 2

% using the final point from FW
[deltaVecLower,exitFlagStep2FixedLength,extraOut] = subProblemFunc(vecFWLower,gradFuncLower(vecFWLower));
if exitFlagStep2FixedLength == SubProblemExitFlag.failed
    error("Could not compute gap for step 2.")
end
gapLower = -FrankWolfe.inProdR(deltaVecLower,gradFuncLower(vecFWLower));
relEntLowerBoundL = funcLower(vecFWLower) - abs(gapLower);

% using the best lowerbound point from FW (currently we trust our func and
% gradfunc enough to use the value output by FW directly)
if outputLower.lowerBoundFWVal > relEntLowerBoundL
    fprintf("!!A better point was found!!\nold %e, new %e\n",relEntLowerBoundL,outputLower.lowerBoundFWVal)
end

relEntLowerBoundL = max(relEntLowerBoundL,outputLower.lowerBoundFWVal);
debugInfo.storeInfo("relEntLowerBoundL",relEntLowerBoundL);

%% Run FW with upper bound on probability in generation rounds
funcUpper =     @(vecFW) RenyiTildeDownDecoyPA.funcFW(vecFW,probTest,entropyAlpha,perturbation,perturbationAff, ...
    dimAPrime,params.stateGen,params.keyProj,params.krausOps,params.photonDistributionGenUpper);
gradFuncUpper = @(vecFW) RenyiTildeDownDecoyPA.gradFuncFW(vecFW,probTest,entropyAlpha,perturbation,perturbationAff, ...
    dimAPrime,params.stateGen,params.keyProj,params.krausOps,params.photonDistributionGenUpper);


% run step 1 routines
switch options.frankWolfeMethod
    case "vanilla"
        [vecFWUpper,relEntStep1Upper,exitFlag,outputUpper] = FrankWolfe.vanilla(vecFWInit, funcUpper,...
            gradFuncUpper, subProblemFunc, debugInfo, options.verboseLevel>=1,...
            "linearSearchMinStep",options.linearSearchMinStep, ...
            "linearSearchPrecision",options.linearSearchPrecision,...
            "maxGap",options.maxGap,...
            "maxIter",options.maxIter,...
            "stopNegGap",options.stopNegGap,...
            "storePoints",true);
    case "pairwise"
        [vecFWUpper,relEntStep1Upper,exitFlag,outputUpper] = FrankWolfe.pairwise(vecFWInit, funcUpper,...
            gradFuncUpper, subProblemFunc, debugInfo, options.verboseLevel>=1,...
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

debugInfo.storeInfo("relEntStep1Upper",relEntStep1Upper);

% step 2

% using the final point from FW
[deltaVecUpper,exitFlagStep2FixedLength,extraOut] = subProblemFunc(vecFWUpper,gradFuncUpper(vecFWUpper));
if exitFlagStep2FixedLength == SubProblemExitFlag.failed
    error("Could not find an initial Point for Frank Wolfe solver.")
end
gapUpper = -FrankWolfe.inProdR(deltaVecUpper,gradFuncUpper(vecFWUpper));
relEntLowerBoundU = funcUpper(vecFWUpper) - abs(gapUpper);

% using the best lowerbound point from FW (currently we trust our func and
% gradfunc enough to use the value output by FW directly)
if options.verboseLevel >= 1 && outputUpper.lowerBoundFWVal > relEntLowerBoundU
    fprintf("!!A better point was found!!\nold %e, new %e\n",relEntLowerBoundU,outputUpper.lowerBoundFWVal)
end

relEntLowerBoundU = max(relEntLowerBoundU,outputUpper.lowerBoundFWVal);
debugInfo.storeInfo("relEntLowerBoundU",relEntLowerBoundU);

%% final lower bound on relative entropy is minimum of both branches
[relEntLowerBound,indexMin] = min([relEntLowerBoundL,relEntLowerBoundU]);
debugInfo.storeInfo("relEntLowerBound",relEntLowerBound);

if options.QESMode == true
    %add debug info
    debugQES = debugInfo.addLeaves("QES");

    %% Step 1 for QES

    %optimal linearization point from fixed length calculation
    outputList = [outputLower,outputUpper];
    optimalLinPoint = outputList(indexMin).lowerBoundFWPoint;
    optimalLinPoint = QESFinderDecoy.perturbEigvalsChoi(optimalLinPoint,dimAPrimeB,false);
    
    % find gradient of QES at linearization point
    [qesGrad,tolGrad,exitFlag] = QESFinderDecoy.getGradQESIntImperfect(optimalLinPoint,gradFuncLower(optimalLinPoint),...
        params.stateTest,params.stateGen,dimAPrime,probTest,testCons,squashingConsTest,squashingConsGen,...
        params.probDistPhotonConMuTestLower, params.probDistPhotonConMuTestUpper, params.probRemaining, ...
        params.probDecoyConTest,params.probSignalConTest,params.blockPhotonNum,options);
    if exitFlag == SubProblemExitFlag.failed
        error("Could not get gradient for QES.")
    end
    
    %Subproblem func equal for both cases
    subProblemFuncQes = @(vecQesFW,vecQesGrad) QESFinderDecoy.subproblemQesIntImperfect(vecQesFW,vecQesGrad,...
        params.stateTest,params.stateGen,dimAPrime,testCons,squashingConsTest,squashingConsGen,...
        params.probDistPhotonConMuTestLower, params.probDistPhotonConMuTestUpper, params.probRemaining,...
        params.probDecoyConTest,params.probSignalConTest,params.blockPhotonNum,options);

    %Initial point for FW iteration
    [vecQesInit,exitFlag] = QESFinderDecoy.initPointIntImperfect(optimalLinPoint(1:end-1),params.stateTest,params.stateGen,...
        dimAPrime,testCons,squashingConsTest,squashingConsGen,...
        params.probDistPhotonConMuTestLower, params.probDistPhotonConMuTestUpper, params.probRemaining, ...
        params.probDecoyConTest,params.probSignalConTest,params.blockPhotonNum,options);
    if exitFlag == SubProblemExitFlag.failed
        error("Could not find an initial Point for QES Frank Wolfe solver.")
    end

    %% Run FW for QES with lower bound on probability in generation rounds

    %Define functions for finding constant term in Qes (required for FW
    %algorithm)
    funcQesLower =   @(vecQesFW) QESFinderDecoy.funcFW(vecQesFW,params.stateGen,params.keyProj,params.krausOps,testCons, ...
        probTest,qesGrad,entropyAlpha,dimAPrime,params.probDecoyConTest,params.photonDistributionGenLower,perturbation,perturbationAff);

    gradFuncQesLower =  @(vecQesFW) QESFinderDecoy.gradFuncFW(vecQesFW,params.stateGen,params.keyProj,params.krausOps,testCons,...
        probTest,qesGrad,entropyAlpha,dimAPrime,params.probDecoyConTest,params.photonDistributionGenLower,perturbation,perturbationAff);
    
    %Run FW solver
    switch options.frankWolfeMethod
        case "vanilla"
            [vecQesLower,qesConstLower,exitFlagQes,outputQesLower] = FrankWolfe.vanilla(vecQesInit, funcQesLower,...
                gradFuncQesLower, subProblemFuncQes, debugQES, options.verboseLevel>=1,...
                "linearSearchMinStep",options.linearSearchMinStep, ...
                "linearSearchPrecision",options.linearSearchPrecision,...
                "maxGap",options.maxGap,...
                "maxIter",options.maxIter,...
                "storePoints",false);
        case "pairwise"
            [vecQesLower,qesConstLower,exitFlagQes,outputQesLower] = FrankWolfe.pairwise(vecQesInit, funcQesLower,...
                gradFuncQesLower, subProblemFuncQes, debugQES, options.verboseLevel>=1,...
                "tolVert",options.tolVert,...
                "linearSearchPrecision",options.linearSearchPrecision,...
                "maxGap",options.maxGap,...
                "maxIter",options.maxIter,...
                "storePoints",false);
    end
    debugInfo.storeInfo("qesConstStep1L",qesConstLower);

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

    % step 2 for QES    
    %perturb result from step 1 if necessary
    vecQesLower = QESFinderDecoy.perturbEigvalsChoi(vecQesLower,dimAPrimeB,true);
    
    % Run one more FW to find reliable lower bound
    [deltaVecQesLower,exitFlagStep2Qes,extraOut] = subProblemFuncQes(vecQesLower,gradFuncQesLower(vecQesLower));
    if exitFlagStep2Qes == SubProblemExitFlag.failed
        error("Could not compute gap for QES step 2.")
    end
    gapQes = -FrankWolfe.inProdR(deltaVecQesLower,gradFuncQesLower(vecQesLower));
    qesConstLowerBoundL = funcQesLower(vecQesLower) - abs(gapQes);
    % qesConstLowerBound = -1/(alphaHat-1)*log2(-(funcQes(vecQes) - abs(gapQes)));
    
    %Pick best lowerbound
    if options.verboseLevel >= 1 && outputQesLower.lowerBoundFWVal > qesConstLowerBoundL
        fprintf("!!A better QES constant was found!!\nold %e, new %e\n",qesConstLowerBoundL,outputQesLower.lowerBoundFWVal);
    end
    qesConstLowerBoundL = max(qesConstLowerBoundL,outputQesLower.lowerBoundFWVal);

    %% Run FW for QES with Upper bound on probability in generation rounds
    %Define functions for finding constant term in Qes (required for FW
    %algorithm)
    funcQesUpper =   @(vecQesFW) QESFinderDecoy.funcFW(vecQesFW,params.stateGen,params.keyProj,params.krausOps,testCons, ...
                probTest,qesGrad,entropyAlpha,dimAPrime,params.probDecoyConTest,params.photonDistributionGenUpper,perturbation,perturbationAff);

    gradFuncQesUpper =  @(vecQesFW) QESFinderDecoy.gradFuncFW(vecQesFW,params.stateGen,params.keyProj,params.krausOps,testCons,...
                probTest,qesGrad,entropyAlpha,dimAPrime,params.probDecoyConTest,params.photonDistributionGenUpper,perturbation,perturbationAff);
    
    %Run FW solver
    switch options.frankWolfeMethod
        case "vanilla"
            [vecQesUpper,qesConstUpper,exitFlagQes,outputQesUpper] = FrankWolfe.vanilla(vecQesInit, funcQesUpper,...
                gradFuncQesUpper, subProblemFuncQes, debugQES, options.verboseLevel>=1,...
                "linearSearchMinStep",options.linearSearchMinStep, ...
                "linearSearchPrecision",options.linearSearchPrecision,...
                "maxGap",options.maxGap,...
                "maxIter",options.maxIter,...
                "storePoints",false);
        case "pairwise"
            [vecQesUpper,qesConstUpper,exitFlagQes,outputQesUpper] = FrankWolfe.pairwise(vecQesInit, funcQesUpper,...
                gradFuncQesUpper, subProblemFuncQes, debugQES, options.verboseLevel>=1,...
                "tolVert",options.tolVert,...
                "linearSearchPrecision",options.linearSearchPrecision,...
                "maxGap",options.maxGap,...
                "maxIter",options.maxIter,...
                "storePoints",false);
    end
    debugInfo.storeInfo("qesConstStep1U",qesConstUpper);

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

    % step 2 for QES    
    %perturb result from step 1 if necessary
    vecQesUpper = QESFinderDecoy.perturbEigvalsChoi(vecQesUpper,dimAPrimeB,true);
    
    % Run one more FW to find reliable lower bound
    [deltaVecQesUpper,exitFlagStep2Qes,extraOut] = subProblemFuncQes(vecQesUpper,gradFuncQesUpper(vecQesUpper));
    if exitFlagStep2Qes == SubProblemExitFlag.failed
        error("Could not compute gap for QES step 2.")
    end
    gapQes = -FrankWolfe.inProdR(deltaVecQesUpper,gradFuncQesUpper(vecQesUpper));
    qesConstLowerBoundU = funcQesUpper(vecQesUpper) - abs(gapQes);
    % qesConstLowerBound = -1/(alphaHat-1)*log2(-(funcQes(vecQes) - abs(gapQes)));
    
    %Pick best lowerbound
    if options.verboseLevel >= 1 && outputQesUpper.lowerBoundFWVal > qesConstLowerBoundU
        fprintf("!!A better QES constant was found!!\nold %e, new %e\n",qesConstLowerBoundL,outputQesLower.lowerBoundFWVal);
    end
    qesConstLowerBoundU = max(qesConstLowerBoundU,outputQesUpper.lowerBoundFWVal);
    
    %% Store Final QES as a cell array
    qesConstLowerBound = min(qesConstLowerBoundL,qesConstLowerBoundU);
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
