function qkdInput = RenyiBB84LossyPreset()
% BasicBB84LossyPreset a preset for a simple BB84 protocol with a lossy
% channel. Loss is included as a third dimension orthogonal to Here,
% Schmidt decomposition was used to shrink Alice from a 4d space to a 2d
% space.

qkdInput = QKDSolverInput();

%% Channel Parameters

%Add misalignment
qkdInput.addFixedParameter("misalignmentAngle",0);
%qkdInput.addScanParameter("misalignmentAngle",num2cell(linspace(0,pi/4,11)));

%Add depolarization
qkdInput.addFixedParameter("depolarization", 0.03); %depol = 0.03

% %Add loss
% lossdB = linspace(0,50,26);
% transmittance = 10.^(-lossdB/10);
% % transmittance = transmittance(1);
% qkdInput.addScanParameter("transmittance", num2cell(transmittance));
% qkdInput.addFixedParameter("transmittance", 1);

%% Define basis choice probability of Alice and Bob (Z basis)
%fixed probTest
% qkdInput.addFixedParameter("probTest",0.26);

%scan probTest
% qkdInput.addScanParameter("probTest", num2cell(linspace(0.0001,0.05,8)))

% %optimzed probTest
% probTest.lowerBound = 0.01;
% probTest.upperBound = 0.5;
% probTest.initVal = 0.1;
% qkdInput.addOptimizeParameter("probTest", probTest);

%Error correction efficiency f >= 1, f=1 is at Shanon limit
qkdInput.addFixedParameter("fEC",1.1);

%% Alpha parameter of Renyi entropy
%fixed alpha
% qkdInput.addFixedParameter("logrenyiAlpha",-1.5);

%scan alpha
% qkdInput.addScanParameter("logrenyiAlpha", num2cell(linspace(-5,-0.6,20)))

% % %optimzed alpha
% logrenyiAlpha.lowerBound = -5;
% logrenyiAlpha.upperBound = -0.5;
% logrenyiAlpha.initVal = -4;
% qkdInput.addOptimizeParameter("logrenyiAlpha", logrenyiAlpha);

%%  Finite-Size Correction parameters

% qkdInput.addFixedParameter("Ntot", 1e10); %total number of signals sent

epsilon.EC = (1/2)*1e-80; %failure probability for error-correction
epsilon.PA = (1/2)*1e-80; %failure probability for privacy amplification
qkdInput.addFixedParameter("epsilon", epsilon);


%% Modules
% description 
descriptionModule = QKDDescriptionModule(@RenyiBB84LossyDescriptionFunc);
qkdInput.setDescriptionModule(descriptionModule);

% channel model
channelModule = QKDChannelModule(@RenyiB84LossyChannelFunc);
qkdInput.setChannelModule(channelModule);

% key rate function
keyModule = QKDKeyRateModule(@RenyiBB84KeyRateFunc);
qkdInput.setKeyRateModule(keyModule);

% optimization

% coordinate descent
optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0,"maxIterations",3,"linearResolution",8),...
    struct("verboseLevel",0));

% %direct search
% optimizerMod = QKDOptimizerModule(@directSearchOptimization, ...
%     struct("verboseLevel",0,"linearResolution",50,"meshSize",0.5),struct("verboseLevel",0));

qkdInput.setOptimizerModule(optimizerMod);

% math solver options
mathSolverOptions = struct();

%Set QES mode = true for variable length
mathSolverOptions.QESMode = true;

mathSolverOptions.initMethod = 1;
mathSolverOptions.maxIter = 30;
mathSolverOptions.maxGap = 1e-10;
mathSolverOptions.blockDiagonal = true;
mathSolverOptions.numericPerturbation = 1e-10;
mathSolverOptions.numericPerturbationAffine = 0;
mathSolverOptions.linearConstraintTolerance = 1e-8;
mathSolverOptions.frankWolfeMethod = 'vanilla';
mathSolverOptions.stopNegGap = true;
mathSolverMod = QKDMathSolverModule(@renyiFrankWolfeSolver,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

%% global options
qkdInput.setGlobalOptions(struct("errorHandling",ErrorHandling.DontCatch,"verboseLevel",1,"cvxSolver","mosek", "cvxPrecision", "high"));