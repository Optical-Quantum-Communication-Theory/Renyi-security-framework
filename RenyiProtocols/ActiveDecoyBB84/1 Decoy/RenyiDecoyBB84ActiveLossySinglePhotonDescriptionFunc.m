function [newParams,modParser] = RenyiDecoyBB84ActiveLossySinglePhotonDescriptionFunc(params, options, debugInfo)

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
modParser.parse(params)
params = modParser.Results;

%% simple setup
%Setup parameters for later use
newParams = struct();

%z-basis choice
pz = params.probTest;

%Define + and - states for later use
ketP = [1;1]/sqrt(2);
ketM = [1;-1]/sqrt(2);

dimA = 2;
dimAPrime = 2;
dimB = 3;

newParams.dimA = dimA;
newParams.dimAPrime = dimAPrime;
newParams.dimB = dimB;

%% Test and Gen states
% Defines states Alice send in test and generation rounds
vecSent = MaxEntangled(dimAPrime,0,1);
vecSent = (vecSent*vecSent');

newParams.stateTest = {vecSent/trace(vecSent)};
newParams.stateGen = {vecSent/trace(vecSent)};

%Probabilities of sending each signal conditioned on test and generation round
newParams.probSignalConTest = [1/2;1/2];
newParams.probSignalConGen = [1/2;1/2];

%% joint obserables
%ordered H,V and CONDITIONED ON GEN
POVMsAGen = {diag([1,0]),diag([0,1])};

%ordered D,A and CONDITIONED ON TEST
POVMsATest = {ketP*ketP',ketM*ketM'};

newParams.POVMATest = POVMsATest;
newParams.POVMAGen = POVMsAGen;

% include vacuum
% block diagonal structure and order: 2x2 qubit, 1x1 vac
% ordered H,V,D,A, vac (no detection)
POVMsB = {pz*diag([1,0,0]),pz*diag([0,1,0]),(1-pz)*([ketP;0]*[ketP;0]'),(1-pz)*([ketM;0]*[ketM;0]'), diag([0,0,1])};

newParams.POVMB = POVMsB;

% each POVM element is assigned an announcement made by their respective
% party
newParams.announcementsA = ["Z","Z"];
newParams.announcementsB = ["Z","Z","X","X","vac"];

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

krausOpZ = sqrt(pz)*kron(diag([1,0])+diag([0,1]),eye(dimA,dimB));
krausOps = {krausOpZ};

krausSum = 0;
for index = 1:numel(krausOps)
    krausSum = krausOps{index}'*krausOps{index};
end
debugInfo.storeInfo("krausSum", krausSum);

%% Pinching map 
proj0 = kron(diag([1,0]),eye(dimB-1));
proj1 = kron(diag([0,1]),eye(dimB-1));
keyProj = {proj0,proj1};

%% set key map, kraus ops, and block diagonal structure in new parameters
newParams.krausOps = {krausOps};
newParams.keyProj = {keyProj};
newParams.blockDimsAPrime = {dimAPrime};
newParams.blockDimsA = [2];
newParams.blockDimsB = [2,1];
newParams.blockPhotonNum = 1;
end

