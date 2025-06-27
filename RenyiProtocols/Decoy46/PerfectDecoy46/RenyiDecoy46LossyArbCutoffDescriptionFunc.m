function [newParams,modParser] = RenyiDecoy46LossyArbCutoffDescriptionFunc(params, options, debugInfo)

arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
%Parsing technical options for the module
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser
%Parsing parameters for the module
modParser = moduleParser(mfilename);
modParser.addRequiredParam("probTest",@(x) mustBeInRange(x,0,1));
modParser.addRequiredParam("probsB",@(x) mustBeProbDist(x));
modParser.addRequiredParam("NB",@(x) mustBeInteger(x));
modParser.addAdditionalConstraint(@(x) mustBeGreaterThanOrEqual(x,1),"NB");
modParser.parse(params)
params = modParser.Results;

%% simple setup
%Setup parameters for later use
newParams = struct();

%z-basis choice
pz = 1-params.probTest;

%Define + and - states for later use
ketH = [1;0];
ketV = [0;1];
ketD = [1;1]/sqrt(2);
ketA = [1;-1]/sqrt(2);
ketR = [1;1i]/sqrt(2);
ketL = [1;-1i]/sqrt(2);

%Bob's basis choice probabilities
probsB = params.probsB;
pzB = probsB(1);
pxB = probsB(2);
pyB = probsB(3); 

%Dimensions of Alice and her sent state
dimA = 4; % |1,2,3,4>
dimAPrime = 2;

newParams.dimA = dimA;
newParams.dimAPrime = [1,dimAPrime];


%% Test and Gen states
% Defines states Alice send in test and generation rounds

% n=1 photon (Alice)
vecSentGenSinglePhoton = kron(zket(4,1),ketH) + kron(zket(4,2),ketV); %unnormalized!!
vecSentTestSinglePhoton = kron(zket(4,3),ketD) + kron(zket(4,4),ketA); %unnormalized!!

% n=0 photon (Alice)
rhoSentGenZeroPhoton = 1/2*blkdiag(ones(2),zeros(2));
rhoSentTestZeroPhoton = 1/2*blkdiag(zeros(2),ones(2));

newParams.stateGen = {rhoSentGenZeroPhoton,1/2*(vecSentGenSinglePhoton*vecSentGenSinglePhoton')};
newParams.stateTest = {rhoSentTestZeroPhoton,1/2*(vecSentTestSinglePhoton*vecSentTestSinglePhoton')};

%Probabilities of sending each signal conditioned on test and generation round
newParams.probSignalConTest = [1/2;1/2];
newParams.probSignalConGen = [1/2;1/2];

%% joint obserables
%ordered H,V and CONDITIONED ON GEN
POVMsAGen = {diag([1,0,0,0]), diag([0,1,0,0])};

%ordered D,A and CONDITIONED ON TEST
POVMsATest = {diag([0,0,1,0]), diag([0,0,0,1])};

newParams.POVMAGen = POVMsAGen;
newParams.POVMATest = POVMsATest;


% block diagonal structure and order: 1x1 vac, 2x2 qubit, 3x3 2-photon ... , flags
%order of POVM elements S_( H V D A R L ) Multi-click NON
blockDimsBob = [(1:params.NB+1),ones(1,7)];

%Bob's Povms generated in separate function
POVMsB = Create_POVM_6State_Flag(params.NB,probsB);

newParams.POVMB = POVMsB;
newParams.BobMultiPOVM = POVMsB{end-1};

%Dimensions of Bob's side
subspacedim = (params.NB+1)*(params.NB+2)/2;
dimB = subspacedim + 7;
newParams.dimB = dimB;

%projector onto Bob's <=N_B photon subspace
newParams.BobSubspaceProj = blkdiag(diag(ones(subspacedim,1)), zeros(7));

% each POVM element is assigned an announcement made by their respective
% party
newParams.announcementsA = ["Z","Z"];
newParams.announcementsB = ["Z","Z","X","X","Y","Y","Mult","vac"];

% For each pair of announcements Alice and Bob keep, we assign Alice's POVM
% elements a corresponding key dit value.
newParams.keyMap = KeyMapElement("Z","Z",[1,2]);

% add all joint observables conditioned on test rounds
observablesJointTest = cell(numel(POVMsATest),numel(POVMsB));
for indexA = 1:numel(POVMsATest)
    for indexB = 1:numel(POVMsB)
        observablesJointTest{indexA,indexB} = kron(POVMsATest{indexA},POVMsB{indexB});
    end
end
newParams.observablesJointTest = observablesJointTest;

% add all joint observables conditioned on gen rounds
observablesJointGen = cell(numel(POVMsAGen),numel(POVMsB));
for indexA = 1:numel(POVMsAGen)
    for indexB = 1:numel(POVMsB)
        observablesJointGen{indexA,indexB} = kron(POVMsAGen{indexA},POVMsB{indexB});
    end
end
newParams.observablesJointGen = observablesJointGen;

%% Kraus Ops (for G map)
% A: Alice's system, B: Bob's System, C: announcement register, R: key
% register. The Kraus operators are matrices from ABC \rightarrow RBC. We
% use results from https://arxiv.org/abs/1905.10896 (Appendix A) to shrink
% the Kraus operators from outputing RABC to just RBC.


%Bob is structered as 1x1 vac, 2x2 qubit, 3x3 2-photon, ... , flags
sqrtHV = sqrtm(POVMsB{1} + POVMsB{2});
sqrtHV = 1/2*(sqrtHV + sqrtHV');

%combine to Kraus op
krausOpZ = kron(zket(2,1)*zket(4,1)'+zket(2,2)*zket(4,2)',sqrtHV);
krausOps = {krausOpZ};

krausSum = 0;
for index = 1:numel(krausOps)
    krausSum = krausOps{index}'*krausOps{index};
end
debugInfo.storeInfo("krausSum", krausSum);

%% Pinching map 
proj0 = kron(diag([1,0]),eye(dimB));
proj1 = kron(diag([0,1]),eye(dimB));
keyProj = {proj0,proj1};

%% set key map, kraus ops, and block diagonal structure in new parameters
newParams.krausOps = {krausOps,krausOps};
newParams.keyProj = {keyProj,keyProj};
newParams.blockDimsAPrime = {1,dimAPrime};
newParams.blockDimsA = [dimA];
newParams.blockDimsB = blockDimsBob;
newParams.blockPhotonNum = [0,1];
end

