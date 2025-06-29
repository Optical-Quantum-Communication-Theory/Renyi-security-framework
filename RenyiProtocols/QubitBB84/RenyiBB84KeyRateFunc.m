function [keyRate, modParser] = RenyiBB84KeyRateFunc(params,options,mathSolverFunc,debugInfo)
arguments
    params (1,1) struct
    options (1,1) struct
    mathSolverFunc (1,1) function_handle
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;


%% modParser
modParser = moduleParser(mfilename);

modParser.addRequiredParam("observablesJointTest",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("expectationsJointTest",@mustBeProbDist);
modParser.addAdditionalConstraint(@isEqualSize,["observablesJointTest","expectationsJointTest"]);

modParser.addRequiredParam("observablesJointGen",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("expectationsJointGen",@mustBeProbDist);
modParser.addAdditionalConstraint(@isEqualSize,["observablesJointGen","expectationsJointGen"]);

modParser.addRequiredParam("krausOps", @isCPTNIKrausOps);
modParser.addRequiredParam("keyProj", @(x) mustBeAKeyProj(x));

modParser.addRequiredParam("dimA",@mustBeInteger);
modParser.addRequiredParam("dimAPrime",@mustBeInteger);
modParser.addRequiredParam("dimB", @mustBeInteger);
modParser.addAdditionalConstraint(@mustBePositive,"dimA");
modParser.addAdditionalConstraint(@mustBePositive,"dimAPrime");
modParser.addAdditionalConstraint(@mustBePositive,"dimB");
modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJointTest","dimA","dimB"]);
modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJointGen","dimA","dimB"]);

modParser.addRequiredParam("announcementsA")
modParser.addRequiredParam("announcementsB")
modParser.addRequiredParam("keyMap",@(x)mustBeA(x,"KeyMapElement"))
modParser.addAdditionalConstraint(@mustBeSizedLikeAnnouncements,["expectationsJointTest","announcementsA","announcementsB"])
modParser.addAdditionalConstraint(@mustBeSizedLikeAnnouncements,["expectationsJointGen","announcementsA","announcementsB"])

modParser.addRequiredParam("fEC", @(x) mustBeGreaterThanOrEqual(x,1));

modParser.addRequiredParam("probTest",@(x) x>=0);
% modParser.addRequiredParam("renyiAlpha",@(x)mustBeInRange(x,1,2,"exclude-lower"));
modParser.addRequiredParam("logrenyiAlpha",@(x)mustBeInRange(x,-10,0,"exclude-lower"));
modParser.addOptionalParam("epsilon", @(x) mustBeInRange(x, 0, 1));
modParser.addRequiredParam("Ntot",@(x) x>0);

modParser.addRequiredParam("stateTest",@isDensityOperator);
modParser.addRequiredParam("stateGen",@isDensityOperator);

modParser.addOptionalParam("blockDimsAPrime", nan, @isBlockDimsWellFormated);
modParser.addOptionalParam("blockDimsB", nan, @isBlockDimsWellFormated);
modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsAPrime","dimA"]);
modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsB","dimB"]);
modParser.addAdditionalConstraint(@(blockDimsAPrime,blockDimsB) ~xor(isequaln(blockDimsAPrime,nan),isequaln(blockDimsB,nan)),["blockDimsAPrime","blockDimsB"]);
modParser.parse(params);

params = modParser.Results;

%% simple setup
debugMathSolver = debugInfo.addLeaves("mathSolver");
mathSolverInput = struct();

%% Error correction

deltaLeak = errorCorrectionCost(params.announcementsA,params.announcementsB,...
    params.expectationsJointGen,params.keyMap,params.fEC);
debugInfo.storeInfo("deltaLeak",deltaLeak);


%% translate for the math solver
%now we have all the parts we need to get a key rate from the a math
%solver, but we need to put it into a form it can understand.
%first we give it the kraus operators for the G map and the projection
%operators for the key map (Z).
mathSolverInput.krausOps = params.krausOps;
mathSolverInput.keyProj = params.keyProj;

%Recast logAlpha to alpha
params.renyiAlpha = 1 + 10^(params.logrenyiAlpha);

mathSolverInput.renyiAlpha = params.renyiAlpha;
mathSolverInput.dimAPrime = params.dimAPrime;

mathSolverInput.stateTest = params.stateTest;
mathSolverInput.stateGen = params.stateGen;


numObs = numel(params.observablesJointTest);

mathSolverInput.testConstraints = arrayfun(@(x)...
    EqualityConstraint(params.observablesJointTest{x},params.expectationsJointTest(x)),1:numObs);
mathSolverInput.probTest = params.probTest;

% if block diag information was give, then pass it to the solver.
if ~isequaln(params.blockDimsAPrime,nan)
    mathSolverInput.blockDimsAPrime = params.blockDimsAPrime;
    mathSolverInput.blockDimsB = params.blockDimsB;
end

% now we call the math solver function on the formulated inputs, with
% it's options.
[relEnt,~] = mathSolverFunc(mathSolverInput,debugMathSolver);

%store the key rate (even if negative)
keyRateFixed = finiteKeyRate(relEnt, deltaLeak, params.renyiAlpha, params.epsilon, params.Ntot, 1-params.probTest, debugInfo);
debugInfo.storeInfo("keyRateFixed",keyRateFixed);

%extract QES from debugInfo if toggled
if isfield(debugInfo.leaves.mathSolver.info,"qesCell")
    %extract qes
    qesCell = debugInfo.leaves.mathSolver.info.qesCell;
    %test constraints
    testCons = [mathSolverInput.testConstraints.scalar];
    %Lower bound on relative entropy from QES
    QESRelEnt = qesCell{1}.'*([params.probTest*testCons.';1-params.probTest]) + qesCell{2};

    %Key rate from QES
    QESRate = finiteKeyRate(QESRelEnt, deltaLeak, params.renyiAlpha, params.epsilon, params.Ntot, 1-params.probTest, debugInfo);

    %Store QES Rate in debuginfo
    debugInfo.storeInfo("QESRate",QESRate);
    debugInfo.storeInfo("qesCell",qesCell);

    %set key rate to QES rate
    keyRate = QESRate;
else
    keyRate = keyRateFixed;
end

if options.verboseLevel>=1
    %ensure that we cut off at 0 when we display this for the user.
    fprintf("Key rate: %e\n",max(keyRateFixed,0));

    if isfield(debugInfo.leaves.mathSolver.info,"qesCell")
        fprintf("QES rate = %e and difference from fixed length = %e \n", max(QESRate,0), keyRateFixed - QESRate);
    end
end

%set the linearization estimate key rate as well for debuging
if isfield(debugMathSolver.info,"relEntStep2Linearization")
    keyRateStep2Linearization = debugMathSolver.info.relEntStep2Linearization - deltaLeak; 
    debugInfo.storeInfo("keyRateRelEntStep2Linearization",keyRateStep2Linearization)

    if options.verboseLevel>=2
        fprintf("Key rate using step 2 linearization intial value: %e\n",max(keyRateStep2Linearization,0))
    end
end

%%%%%%%%%%%  FINITE Key Rate Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finite-size key rate function combining results from 2 works
% 
% Insert eq. (20) into Theorem 1 from: 
% Finite-size analysis of prepare-and-measure and decoy-state QKD via entropy accumulation 
% http://arxiv.org/abs/2406.10198
%
% and then use the lower bound on the Renyi entropy from Lemma 10 of:
% Generalized R\'enyi entropy accumulation theorem and generalized quantum
% probability estimation
% http://arxiv.org/abs/2405.05912

    function [keyRate, debugInfo] = finiteKeyRate(relEnt, deltaLeak, renyiAlpha, epsilon, Ntotal, pGen, debugInfo)
    %computes the finite size keyrate.

    %Error correction including error verification
    ECLeakage = pGen*deltaLeak + ceil(log2(1/epsilon.EC))/Ntotal;
    
    %Privacy Amplification
    privacyAmplification = renyiAlpha/(renyiAlpha-1)*log2(1/epsilon.PA)/Ntotal;

    keyRate = relEnt - ECLeakage - privacyAmplification + 2/Ntotal;
end

end