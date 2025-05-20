function qkdInput = RenyiDecoyBB84PassivePreset()
% BasicBB84LossyPreset a preset for a simple BB84 protocol with a lossy
% channel. Loss is included as a third dimension orthogonal to Here,
% Schmidt decomposition was used to shrink Alice from a 4d space to a 2d
% space.

qkdInput = QKDSolverInput();

%% Parameters

%Add loss
lossdB = linspace(0,60,16);
% transmittance = 10.^(-lossdB/10);
% transmittance = transmittance(1);
% transmittance = 10.^(-2.8);
% qkdInput.addScanParameter("transmittance", num2cell(transmittance));

%Add birefringence
% qkdInput.addFixedParameter("birefringenceAngle", 0);
% qkdInput.addScanParameter("birefringenceAngle", num2cell(linspace(0, pi/4, 11)));

%Add misalignment
qkdInput.addFixedParameter("misalignmentAngle",0.03);
%qkdInput.addScanParameter("misalignmentAngle",num2cell(linspace(0,pi/4,11)));

%Add dark counts
% qkdInput.addFixedParameter("darkCountRate", 0);

%Add depolarization
% qkdInput.addFixedParameter("depolarization", 0);

%Add detector efficiencies as vector
qkdInput.addFixedParameter("detectorEfficiency", ones(4,1));

%Decoy settings
decoy1 = 0.8;
decoy2 = 0.2;
decoy3 = 0.01;

delta = 0;

%Decoys 
qkdInput.addFixedParameter("GROUP_decoys_1", decoy1); %signal intensity
% qkdInput.addOptimizeParameter("GROUP_decoys_1", struct("lowerBound",0,"initVal",0.85,"upperBound",1)); %signal intensity

qkdInput.addFixedParameter("GROUP_decoys_2", decoy2); % decoy intensity 1
% qkdInput.addOptimizeParameter("GROUP_decoys_2", struct("lowerBound",0,"initVal",0.1,"upperBound",1)); % decoy intensity 1

% qkdInput.addFixedParameter("GROUP_decoys_3", decoy3); % decoy intensity 2 (something slightly above 0 is usually optimal.) %0.001 for EAT

%Tolerances in decoy intensities
qkdInput.addFixedParameter("GROUP_deltaDecoys_1", decoy1*(1/(1-delta) - 1)); %signal intensity
qkdInput.addFixedParameter("GROUP_deltaDecoys_2", decoy2*(1/(1-delta) - 1)); %decoy intensity 1
% qkdInput.addFixedParameter("GROUP_deltaDecoys_3", decoy3*(1/(1-delta) - 1)); %decoy intensity 2

%Probabilities of sending decoys conditioned on test rounds p(mu|test)
% qkdInput.addFixedParameter("GROUP_decoyProbs_1",1/3);
% qkdInput.addFixedParameter("GROUP_decoyProbs_2",1/3);
% qkdInput.addFixedParameter("GROUP_decoyProbs_3",1/3);

qkdInput.addFixedParameter("GROUP_decoyProbs_1",1/2);
qkdInput.addFixedParameter("GROUP_decoyProbs_2",1/2);

%Probability of sending intensities in generation rounds 
% qkdInput.addFixedParameter("probDecoyConGen", [1, 0, 0]);
qkdInput.addFixedParameter("probDecoyConGen", [1, 0]);

%Bob's probabilities of choosing each basis
qkdInput.addFixedParameter("probsB", [0.8, 0.2]);

%Bob's photon number cutoff
qkdInput.addFixedParameter("NB", 1);


%Define basis choice probability of Alice (X basis)
qkdInput.addFixedParameter("probTest",0.2);

% probTest.lowerBound = 0.0001;
% probTest.upperBound = 0.9;
% probTest.initVal = 0.05;
% qkdInput.addOptimizeParameter("probTest", probTest);

%Error correction efficiency f >= 1, f=1 is at Shanon limit
qkdInput.addFixedParameter("fEC",1.1);

%Alpha parameter of Renyi entropy
%fixed alpha
% qkdInput.addFixedParameter("logrenyiAlpha",-2);

% %optimzed alpha
% logrenyiAlpha.lowerBound = -5;
% logrenyiAlpha.upperBound = -0.5;
% logrenyiAlpha.initVal = -3;
% qkdInput.addOptimizeParameter("logrenyiAlpha", logrenyiAlpha);

% Finite-Size Correction parameters

% qkdInput.addFixedParameter("Ntot", 1e10); %total number of signals sent

epsilon.EC = (1/2)*1e-80; %failure probability for error-correction
epsilon.PA = (1/2)*1e-80; %failure probability for privacy amplification
qkdInput.addFixedParameter("epsilon", epsilon);


%% Modules
% description 
%Single Photon Description
% descriptionModule = QKDDescriptionModule(@RenyiDecoyBB84PassiveSinglePhotonDescriptionFunc);

%Vacuum and Single Photon Description
descriptionModule = QKDDescriptionModule(@RenyiDecoyBB84PassiveDescriptionFunc);

%Vacuum and Single Photon Description with arbitrary photon number cutoff
%on Bob's side
% descriptionModule = QKDDescriptionModule(@RenyiDecoy4BB84PassiveArbCutoffDescriptionFunc);

qkdInput.setDescriptionModule(descriptionModule);

% channel model
channelModule = QKDChannelModule(@RenyiDecoyBB84PassiveChannelFunc);
qkdInput.setChannelModule(channelModule);

% key rate function
keyModule = QKDKeyRateModule(@RenyiDecoyBB84PassiveKeyRateFunc,struct("photonCutOff",8));
qkdInput.setKeyRateModule(keyModule);

% optimization
%coordinate descent
optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),...
    struct("verboseLevel",0));

% direct search
% optimizerMod = QKDOptimizerModule(@directSearchOptimization, ...
%     struct("verboseLevel",0,"linearResolution",30,"meshSize",0.5),struct("verboseLevel",0));

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
mathSolverOptions.numericPerturbationAffine = 1e-12;
mathSolverOptions.linearConstraintTolerance = 1e-8;
mathSolverOptions.stopNegGap = true;
mathSolverOptions.linearSearchMinStep = 0;
mathSolverOptions.frankWolfeMethod = "vanilla";
mathSolverMod = QKDMathSolverModule(@renyiDecoyIntImpFrankWolfeSolver,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

%% global options
qkdInput.setGlobalOptions(struct("errorHandling",ErrorHandling.DontCatch,"verboseLevel",1,"cvxSolver","mosek", "cvxPrecision", "high"));