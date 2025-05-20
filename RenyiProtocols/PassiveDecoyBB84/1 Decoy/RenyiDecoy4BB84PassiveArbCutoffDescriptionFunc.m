function [newParams,modParser] = RenyiDecoy4BB84PassiveArbCutoffDescriptionFunc(params, options, debugInfo)
% BasicBB84LossyDescriptionFunc A simple description function for a qubit
% BB84 protocol with loss. Here, Schmidt decomposition was used to shrink
% Alice from a 4d space to a 2d space. This is used with BasicKeyRateFunc.
%
% Input parameters:
% * pz: The probability that Alice and Bob measure in the Z-basis. 
% Output parameters:
% * observablesJoint: The joint observables for Alice and Bob's measurement
%   of the signals.
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * probSignalsA: Probability of Alice selecting a signal to send to Bob.
%   In this protocol, it is half the probability of the basis choice.
% * rhoA: Alice's reduced density matrix for prepare-and-measure based
%   protocols.
% * POVMA: Alice's set of POVM operators which she measures her state with
%   in the source replacement scheme.
% * POVMB: Bob's set of POVM operators which he measures his state with.
% * announcementsA: Alice's announcements for each of her POVM operators.
%   Can be integers or strings.
% * announcementsB: Bob's announcements for each of his POVM operators.
%   Can be integers or strings.
% * keyMap: An array of KeyMap objects that contain pairs of accepted
%   announcements and an array dictating the mapping of Alice's measurement
%   outcome to key bits (May be written with Strings).
% * krausOps: A cell array of matrices. The Kraus operators that form the G
%   map on Alice and Bob's joint system. These should form a completely
%   postive trace non-increasing linear map. Each Kraus operator must be
%   the same size.
% * keyProj:  A cell array of projection operators that extract the key
%   from G(\rho). These projection operators should sum to identity. This
%   map is often called Z.
% Options:
% * none
% DebugInfo:
% * krausSum: sum_i K^\dagger_i*K_i which should be <= identity for
%   a CPTNI map.
%
% Reviewed by Devashish Tupkary 2023/09/18
% See also QKDDescriptionModule, BasicKeyRateFunc
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

%Bob's basis choice probabilities
probsB = params.probsB;
pzB = probsB(1);
pxB = probsB(2);

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
blockDimsBob = [(1:params.NB+1),ones(1,5)];

%Bob's Povms generated in separate function
POVMsB = Create_POVM_4State_Flag(params.NB,probsB);

newParams.POVMB = POVMsB;
newParams.BobMultiPOVM = POVMsB{end-1};

%Dimensions of Bob's side
subspacedim = (params.NB+1)*(params.NB+2)/2;
dimB = subspacedim + 5;
newParams.dimB = dimB;

%projector onto Bob's <=N_B photon subspace
newParams.BobSubspaceProj = blkdiag(diag(ones(subspacedim,1)), zeros(5));

% each POVM element is assigned an announcement made by their respective
% party
newParams.announcementsA = ["Z","Z"];
newParams.announcementsB = ["Z","Z","X","X","Mult","vac"];

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

