function [newParams,modParser] = RenyiPhaseImpDecoy46LossyDescriptionFunc(params, options, debugInfo)
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

%testing probability
modParser.addRequiredParam("probTest",@(x) mustBeInRange(x,0,1));

%Bob's basis choices
modParser.addRequiredParam("probsB",@(x) mustBeProbDist(x));

%Bob's photon number cutoff
modParser.addRequiredParam("NB",@(x) mustBeInteger(x));
modParser.addAdditionalConstraint(@(x) mustBeGreaterThanOrEqual(x,1),"NB");

%Alice's photon number cutoff
modParser.addRequiredParam("NA",@(x) mustBeInteger(x));
modParser.addAdditionalConstraint(@(x) mustBeGreaterThanOrEqual(x,1),"NA");

%Decoy Intensities
modParser.addRequiredParam("decoys", @(x) mustBeCellOf(x, 'numeric'));
modParser.addAdditionalConstraint(@(x) allCells(x,@isscalar),"decoys");
modParser.addAdditionalConstraint(@(x) allCells(x,@(y) y>=0),"decoys");

%probs of decoy in test rounds
modParser.addRequiredParam("decoyProbs",@mustBeProbDistCell); %,must sum to 1. 

%probs of decoy in gen rounds
modParser.addRequiredParam("probDecoyConGen",@(x)mustBeProbDist(x));

%Phase randomization quality factor
modParser.addRequiredParam("phaseRandomQuality",@(x) mustBeInRange(x,0,1));

%Block cutoff used for key rate
modParser.addRequiredParam("blockCutoff",@(x) mustBeInteger(x));
modParser.addAdditionalConstraint(@(x) mustBeGreaterThanOrEqual(x,1),"blockCutoff");

modParser.parse(params)

params = modParser.Results;

%% simple setup
%Setup parameters for later use
newParams = struct();

%z-basis choice
pz = 1-params.probTest;

%Convert decoy intensities to matrix
decoys = cell2mat(params.decoys);

%decoy probabilities conditioned on test and gen
probDecoysCondTest = cell2mat(params.decoyProbs);
probDecoysCondGen = params.probDecoyConGen;

%Probabilities of sending each signal conditioned on test and generation round
probsSignalConTest = [1/2;1/2];
probsSignalConGen = [1/2;1/2];

%number of signals in test/gen rounds
numTestSigs = size(probsSignalConTest,1);
numGenSigs = size(probsSignalConGen,1);

%number of decoy intensities
numDecoy = numel(params.decoys);

%Probabilities of sending each signal conditioned on test and generation round
probSignalAndMuConTest = kron(probDecoysCondTest',probsSignalConTest);
newParams.probSignalConTest = probSignalAndMuConTest;
newParams.probSignalConGen = probsSignalConGen;

%tolerance used during calculation of Alice's POVM elements and states
tol = 1e-14;

%% Alice's POVM elements and states
% [rhoTest,rhoGen,POVMsATest,POVMsAGen,dimsA,dimsAPrime,probBlockCondTest,probBlockCondGen,eigenValueMat] = ...
%     create46PhaseImperfectStatesandPOVMsAlice(params.probTest,params.phaseRandomQuality,decoys,probDecoysCondTest,probDecoysCondGen,...
%     probsSignalConTest,probsSignalConGen,params.blockCutoff,params.NA,tol);
% deltaMat = ones((numTestSigs+numTestSigs)*numDecoy,params.blockCutoff+1);

try
    disp("Loading data")
    matData = load("phaseImperfect46DescriptionData.mat");
    [rhoTest,rhoGen,POVMsATest,POVMsAGen,dimsA,dimsAPrime,probBlockCondTest,probBlockCondGen,eigenValueMat,deltaMat] = parseData(matData);
catch
    disp("Loading data failed. Calculating states, POVMs, etc.")
    [rhoTest,rhoGen,POVMsATest,POVMsAGen,dimsA,dimsAPrime,probBlockCondTest,probBlockCondGen,eigenValueMat,deltaMat] = ...
        create46PhaseImperfectStatesandPOVMsAliceWithDelta(params.probTest,params.phaseRandomQuality,decoys,probDecoysCondTest,probDecoysCondGen,...
        probsSignalConTest,probsSignalConGen,params.blockCutoff,params.NA,tol);
    save("phaseImperfect46DescriptionData.mat","rhoTest","rhoGen","POVMsATest","POVMsAGen",...
        "dimsA","dimsAPrime","probBlockCondTest","probBlockCondGen", ...
        "eigenValueMat","deltaMat")
end

%Probabilities of sending each signal conditioned on test/gen-round and block
probSignalAndMuConTestAndBlock = zeros(size(probSignalAndMuConTest,1),params.blockCutoff+1);
probSignalAndMuConGenAndBlock = zeros(size(probsSignalConGen,1),params.blockCutoff+1);

%eigenvalues in test/gen rounds
eigenValueMatTest = eigenValueMat([3:4,5:6],:);
eigenValueMatGen = eigenValueMat(1:2,:);

%test
for indexBlock = 1:params.blockCutoff+1
    probSignalAndMuConTestAndBlock(:,indexBlock) = eigenValueMatTest(:,indexBlock).*probSignalAndMuConTest./probBlockCondTest(indexBlock)';
end

%gen
for indexBlock = 1:params.blockCutoff+1
    probSignalAndMuConGenAndBlock(:,indexBlock) = eigenValueMatGen(:,indexBlock).*probsSignalConGen./probBlockCondGen(indexBlock)';
end

%renormalize probSignalAndMuConTestAndBlock because numerical errors occur at 1e-11, no
%issue since our linear constraint tolerance is 1e-8, thus much bigger
renormprobSignalAndMuConTestAndBlock = probSignalAndMuConTestAndBlock;
for index = 1:size(renormprobSignalAndMuConTestAndBlock,2)
    tolalprob = sum(probSignalAndMuConTestAndBlock(:,index));
    renormprobSignalAndMuConTestAndBlock(:,index) = probSignalAndMuConTestAndBlock(:,index)/tolalprob;
end

%warn if abs(renormprobSignalAndMuConTestAndBlock -
%probSignalAndMuConTestAndBlock) > 1e-8
if max(abs(renormprobSignalAndMuConTestAndBlock - ...
        probSignalAndMuConTestAndBlock),[],"all") > 1e-8
    throwAsCaller(MException("renormProbs:DiftoRenormedProbsTooBig",...
        "Difference to renormalized probabilities too big."))
end

%p(a,mu|test,block)
% newParams.probSignalConTestAndBlock = probSignalAndMuConTestAndBlock;
newParams.probSignalConTestAndBlock = renormprobSignalAndMuConTestAndBlock;

%p(a,mu|gen,block)
newParams.probSignalConGenAndBlock = probSignalAndMuConGenAndBlock;

%Alice's POVM elements odered by each block
newParams.POVMAGen = POVMsAGen;
newParams.POVMATest = POVMsATest;

% dimensions of A and Aprime odered by each block
newParams.dimA = dimsA;
newParams.dimAPrime = dimsAPrime;

%Probabilities of each block conditioned on test and gen
%p(block|test)
newParams.probBlockCondTest = probBlockCondTest;

%p(block|gen)
altprobBlockCondGen = probBlockCondGen;
altprobBlockCondGen(3:end) = 0;

% newParams.probBlockCondGen = probBlockCondGen;
newParams.probBlockCondGen = altprobBlockCondGen;

%states Alice sends in test and generation rounds ordered by each block
newParams.stateTest = rhoTest;
newParams.stateGen = rhoGen;

%% Construction of correction terms

%Construct matrix which copies the weight appropriately for signal states
%in test rounds
Ar = repmat(ones(numTestSigs,1), 1, numDecoy);
Ac = mat2cell(Ar, numTestSigs, ones(1,numDecoy));
MatTest = blkdiag(Ac{:});

%Weight outside conditioned on (a,mu)
[decoysGrid,countGrid] = meshgrid(decoys,0:params.NA);
weightOutcondSignalAndMuAndTest = MatTest*(1-sum(poisson(sym(decoysGrid),sym(countGrid)),1)).';

%spectral gaps in test rounds
deltaMatTest = deltaMat([3:4,5:6],:) - (1-params.phaseRandomQuality)*repmat(sqrt(weightOutcondSignalAndMuAndTest),1,9);

%p(a,mu|test,block)/delta_m(a,mu)
weightedDeltaMatInv = probSignalAndMuConTestAndBlock./deltaMatTest;

%correction term for each block in test rounds 2*(1-q)*sum_{a,mu}
%p(a,mu|test,block)/delta_m(a,mu)sqrt(W_{a,mu})
epsilonBlockCondTest = double(2*(1-params.phaseRandomQuality)*(sqrt(weightOutcondSignalAndMuAndTest)'*weightedDeltaMatInv));


%weight outside conditioned on gen
weigthOutCondSignalAndMuAndGen = ones(numGenSigs,1)*(1-sum(poisson(sym(decoys(1)),sym(0:params.NA))));

%correction term for generation rounds
epsilonFull = double((1-params.phaseRandomQuality)*probsSignalConGen'*sqrt(weigthOutCondSignalAndMuAndGen));


%correction term for probabilities in test rounds
epsilonProbTest = double((1-params.phaseRandomQuality)*probSignalAndMuConTest'*sqrt(weightOutcondSignalAndMuAndTest));


%pass correction terms on to key rate function
newParams.epsilonFull = epsilonFull;
newParams.epsilonBlockCondTest = epsilonBlockCondTest;
newParams.epsilonProbTest = epsilonProbTest;

%% Bob's POVM elements
%Bob's basis choice probabilities
probsB = params.probsB;

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

%% Joint observables
% each POVM element is assigned an announcement made by their respective
% party
newParams.announcementsA = ["Z","Z"];
newParams.announcementsB = ["Z","Z","X","X","Y","Y","Mult","vac"];

% For each pair of announcements Alice and Bob keep, we assign Alice's POVM
% elements a corresponding key dit value.
newParams.keyMap = KeyMapElement("Z","Z",[1,2]);

% add all joint observables conditioned on test rounds
observablesJointTest = cell(size(POVMsATest,1),numel(POVMsB),params.blockCutoff+1);
for indexBlock = 1:params.blockCutoff+1
    for indexA = 1:size(POVMsATest,1)
        for indexB = 1:numel(POVMsB)
            observablesJointTest{indexA,indexB,indexBlock} = kron(POVMsATest{indexA,indexBlock},POVMsB{indexB});
        end
    end
end
newParams.observablesJointTest = observablesJointTest;

% add all joint observables conditioned on gen rounds
observablesJointGen = cell(size(POVMsAGen,1),numel(POVMsB),params.blockCutoff+1);
for indexBlock = 1:params.blockCutoff+1
    for indexA = 1:size(POVMsAGen,1)
        for indexB = 1:numel(POVMsB)
            observablesJointGen{indexA,indexB,indexBlock} = kron(POVMsAGen{indexA,indexBlock},POVMsB{indexB});
        end
    end
end
newParams.observablesJointGen = observablesJointGen;

%% Kraus Ops (for G map)
% A: Alice's system, B: Bob's System, C: announcement register, R: key
% register. The Kraus operators are matrices from ABC \rightarrow RBC. We
% use results from https://arxiv.org/abs/1905.10896 (Appendix A) to shrink
% the Kraus operators from outputing RABC to just RBC.


%Bob is structered as 1x1 vac, 2x2 qubit, 3x3 2-photon, ... , flags
%find sqrt of Bob's POVM
sqrtHVBob = removeNumImprecHermMat(sqrtm(POVMsB{1} + POVMsB{2}));

%find sqrt of Alice's POVMs
% H POVMs
PovmAliceHs = arrayfun(@(x) POVMsAGen{1,x}, 1:1:params.blockCutoff+1,'UniformOutput',false);
sqrtPovmAliceHs = cellfun(@(x) removeNumImprecHermMat(sqrtm(x),tol), PovmAliceHs,'UniformOutput',false);

% V POVMs
PovmAliceVs = arrayfun(@(x) POVMsAGen{2,x}, 1:1:params.blockCutoff+1,'UniformOutput',false);
sqrtPovmAliceVs = cellfun(@(x) removeNumImprecHermMat(sqrtm(x),tol), PovmAliceVs,'UniformOutput',false);

%combine to Kraus op
krausOpZ = arrayfun(@(x) {kron(kron(zket(2,1),sqrtPovmAliceHs{x}) + kron(zket(2,2),sqrtPovmAliceVs{x}),sqrtHVBob)},1:1:params.blockCutoff+1,'UniformOutput',false);
krausOps = krausOpZ;

krausSum = cell(1,size(krausOps,2));

for indexBlock = 1:size(krausOps,2)
    krausSumTemp = 0;
    for index = 1:size(krausOps,1)
        krausSumTemp = krausSumTemp + krausOps{index,indexBlock}{1,1}'*krausOps{index,indexBlock}{1,1};
    end
    krausSum{indexBlock} = krausSumTemp;
end
debugInfo.storeInfo("krausSum", krausSum);

%% Pinching map
keyProj = cell(1,params.blockCutoff+1);
for indexBlock = 1:params.blockCutoff+1
    proj0 = kron(diag([1,0]),eye(dimsA(indexBlock)*dimB));
    proj1 = kron(diag([0,1]),eye(dimsA(indexBlock)*dimB));
    keyProj{indexBlock} = {proj0,proj1};
end

%% set key map, kraus ops, and block diagonal structure in new parameters
newParams.krausOps = krausOps;
newParams.keyProj = keyProj;
newParams.blockDimsAPrime = num2cell(dimsAPrime);
newParams.blockDimsA = dimsA;
newParams.blockDimsB = blockDimsBob;
newParams.blockPhotonNum = 0:1:params.blockCutoff;
end

%% validation functions
function mustBeProbDistCell(input)
    mustBeProbDist([input{:}])
end

%helper function for numerical imprecision of hermitian matrices
function [MatOut,MatIn] = removeNumImprecHermMat(MatIn,tol)
    arguments
        MatIn (:,:) double 
        tol (1,1) double {mustBePositive} = 1e-15;
    end
    %Calculate maximum nonhermitian part
    maxNonHerm = max(abs(MatIn - MatIn'),[],"all");
    
    if maxNonHerm > 1e-6
        throwAsCaller(MException("removeNumImprec:MatrixNotHermitian",...
            "Matrix is not Hermitian."))
    end

    %Set output equal input
    MatOut = 1/2*(MatIn+MatIn');

    %set small values = 0
    MatOut(abs(MatOut) <= tol ) = 0;
    
    %set small imaginary values = 0
    imMat = imag(MatOut);
    imMat(abs(imMat) <= tol ) = 0;

    %set small real values = 0
    reMat = real(MatOut);
    reMat(abs(reMat) <= tol ) = 0;
    
    %put matrix back together
    MatOut = reMat + 1i*imMat;

    %guarantee MatOut is hermitian
    MatOut = 1/2*(MatOut + MatOut');
end

function observablesAndDimensionsMustBeTheSameForAllPages(observables,dimsA,dimsB)
    %number of pages
    numPages = size(observables,3);

    %if dimsA = scalar, repeat same number
    if (size(dimsA,1) == size(dimsA,2)) == 1
        dimsA = dimsA*ones(1,numPages);
    end
    %if dimsB = scalar, repeat same number
    if (size(dimsB,1) == size(dimsB,2)) == 1
        dimsB = dimsB*ones(1,numPages);
    end

    %perform checks
    for index = 1:numPages
        if ~allCells(observables(:,:,index),@(x) size(x,1) == dimsA(index)*dimsB(index))
            throwAsCaller(MException("ValidationFunction:ObservablesHaveMismatchedSize",...
                "The Observables must have the same dimensions as Alice and Bob multiplied together."));
        end
    end
end

function allPages(cellMat,func)
    numPages = size(cellMat,3);
    for index = 1:numPages
        allCells(cellMat(:,:,index),func);
    end
end

%Poisson distribution
function prob = poisson(intensity, count)
    prob = exp(-intensity).*intensity.^count./factorial(count);
end


%helper function parse data
function [rhoTest,rhoGen,POVMsATest,POVMsAGen,dimsA,dimsAPrime,probBlockCondTest,probBlockCondGen,eigenValueMat,probBlockCondTestNPlus] = parseData(matData)
    rhoTest = matData.rhoTest;
    rhoGen = matData.rhoGen;
    POVMsATest = matData.POVMsATest;
    POVMsAGen = matData.POVMsAGen;
    dimsA = matData.dimsA;
    dimsAPrime = matData.dimsAPrime;
    probBlockCondTest = matData.probBlockCondTest;
    probBlockCondGen = matData.probBlockCondGen;
    eigenValueMat = matData.eigenValueMat;
    probBlockCondTestNPlus = matData.deltaMat;
end