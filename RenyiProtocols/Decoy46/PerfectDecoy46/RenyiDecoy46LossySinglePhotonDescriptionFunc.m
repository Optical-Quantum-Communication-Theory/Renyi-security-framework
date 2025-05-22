function [newParams,modParser] = RenyiDecoy46LossySinglePhotonDescriptionFunc(params, options, debugInfo)
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
modParser.parse(params)
params = modParser.Results;

%% simple setup
%Setup parameters for later use
newParams = struct();

%z-basis choice
pz = 1 - params.probTest;

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

%Dimensions
dimA = 4; % |1,2,3,4>
dimAPrime = 2;
dimB = 10; % comb basis + Mult flag

newParams.dimA = dimA;
newParams.dimAPrime = [dimAPrime];
newParams.dimB = dimB;

%% Test and Gen states
% Defines states Alice send in test and generation rounds

% n=1 photon (Alice)
vecSentGenSinglePhoton = kron(zket(dimA,1),ketH) + kron(zket(dimA,2),ketV); %unnormalized!!
vecSentTestSinglePhoton = kron(zket(dimA,3),ketD) + kron(zket(dimA,4),ketA); %unnormalized!!

% n=0 photon (Alice)
rhoSentGenZeroPhoton = 1/2*blkdiag(ones(2),zeros(2));
rhoSentTestZeroPhoton = 1/2*blkdiag(zeros(2),ones(2));

newParams.stateGen = {1/2*(vecSentGenSinglePhoton*vecSentGenSinglePhoton')};
newParams.stateTest = {1/2*(vecSentTestSinglePhoton*vecSentTestSinglePhoton')};

%Probabilities of sending each signal conditioned on test and generation round
newParams.probSignalConTest = [1/2;1/2];
newParams.probSignalConGen = [1/2;1/2];

%% joint obserables
%ordered H,V and CONDITIONED ON GEN
POVMsAGen = {diag([1,0,0,0]), diag([0,1,0,0])}; % this is not a complete set of POVMs

%ordered D,A and CONDITIONED ON TEST
POVMsATest = {diag([0,0,1,0]), diag([0,0,0,1])}; % this is not a complete set of POVMs

newParams.POVMAGen = POVMsAGen;
newParams.POVMATest = POVMsATest;

% include vacuum
% block diagonal structure and order: 2x2 qubit, 1x1 vac, flags
%order of POVM elements S_( H V D A R L ) Multi-click NON
ketD = [ketD; 0]; 
ketA = [ketA; 0];
ketR = [ketR; 0]; 
ketL = [ketL; 0];
GammaH =    blkdiag(pzB*diag([1,0,0]),  diag(zket(7,1)));
GammaV =    blkdiag(pzB*diag([0,1,0]),  diag(zket(7,2)));
GammaD =    blkdiag(pxB*(ketD*(ketD')), diag(zket(7,3)));
GammaA =    blkdiag(pxB*(ketA*(ketA')), diag(zket(7,4)));
GammaR =    blkdiag(pyB*(ketR*(ketR')), diag(zket(7,5)));
GammaL =    blkdiag(pyB*(ketL*(ketL')), diag(zket(7,6)));
GammaMult = blkdiag(zeros(3),           diag(zket(7,7)));
GammaVAC =  blkdiag(diag([0,0,1]),      zeros(7));

POVMsB = {GammaH, GammaV, GammaD, GammaA, GammaR, GammaL, GammaMult, GammaVAC};

newParams.POVMB = POVMsB;
newParams.BobMultiPOVM = GammaMult;
newParams.BobSubspaceProj = blkdiag(diag([1,1,1]), zeros(7));

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
sqrtHV = diag([sqrt(pzB),sqrt(pzB),0,1,1,0,0,0,0,0]);
krausOpZ = kron(zket(2,1)*zket(dimA,1)'+zket(2,2)*zket(dimA,2)',sqrtHV);
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
newParams.krausOps = {krausOps};
newParams.keyProj = {keyProj};
newParams.blockDimsAPrime = {dimAPrime};
newParams.blockDimsA = dimA;
newParams.blockDimsB = [2,1,ones(1,7)];
newParams.blockPhotonNum = 1;
end

