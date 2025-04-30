function [newParams, modParser]= RenyiDecoyBB84PassiveChannelFunc(params,options, debugInfo)
% BasicBB84_WCPActiveChannel A channel function for BB84 using decoy state
% analysis. Given a collection of decoy intensities, this channel produces
% a group of 4x5 tables of squashed expectations, one for each decoy intensity,
% which are the conditional probability for each of Bob's 5 detector
% patterns given Alice's signal sent.
%
% Input parameters:
% * decoys: a cell of the intensities used in decoy analysis. These are the
%   mean photon numbers that Alice can choose from when performing the
%   decoy protocol. The first element in the cell array is treated as the
%   intensity used for key generation.
% * transmittance (1): the transmissivity of the quantum channel; equivalent to 1 -
%   loss. Must be between 0 and 1 inclusive.
% * detectorEfficiency (1): the efficiency of Bob's detectors. Must be
%   between 0 and 1 inclusive.
% * ed (0): Angle Alice and Bob's bases are misaligned by
%   around the Y-axis. For example, Bob's detectors could be slightly
%   rotated away from the incoming signals. Angles must be real numbers.
%   The angle is the physical rotation angle (2pi period) and not an angle
%   on the bloch sphere (4pi period).
% * darkCountRate (0): The probability that a detector that recieves no
%   photons will still randomly click anyway. Must be between 0 and 1.
% Output parameters:
% * expectationsConditional: The conditional expectations (as a 3D array)
%   from Alice and Bob's measurements. This should be organized as a 4 x 16
%   x n array, where n = the number of intensities used in the decoy
%   protocol. In each table, these should line up with the corresponding
%   observables at each entry.
% Options:
% * None.
% DebugInfo:
% * None.
%
% See also QKDChannelModule,  makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser
modParser = moduleParser(mfilename);

%Decoy Intensities
modParser.addRequiredParam("decoys", @(x) mustBeCellOf(x, 'numeric'));
modParser.addAdditionalConstraint(@(x) allCells(x,@isscalar),"decoys");
modParser.addAdditionalConstraint(@(x) allCells(x,@(y) y>=0),"decoys");

%Sending probabilities of decoy intensities
modParser.addRequiredParam("decoyProbs",@mustBeProbDistCell); %,must sum to 1. 

%Bob's basis choices
modParser.addRequiredParam("probsB", @(x) mustBeProbDist(x));

%Transmittance
modParser.addOptionalParam("transmittance", 1, @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"transmittance");

%Detector efficiencies as vector
modParser.addOptionalParam("detectorEfficiency", ones(4,1), @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint( @(x) mustBeEqualSize(x, ones(4,1)) ,"detectorEfficiency");

%Misalignment
modParser.addOptionalParam("misalignmentAngle",0,@mustBeReal);
modParser.addAdditionalConstraint(@isscalar,"misalignmentAngle");

%Birefringence angle
modParser.addOptionalParam("birefringenceAngle",0,@mustBeReal);
modParser.addAdditionalConstraint(@isscalar,"birefringenceAngle");

%Dark count rate
modParser.addOptionalParam("darkCountRate", 0, @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"darkCountRate");

%probability of sending decoy states in generation rounds
modParser.addRequiredParam("probDecoyConGen",@(x)mustBeProbDist(x));

%Bob's photon cutoff
modParser.addOptionalParam("NB",1,@(x) mustBeInteger(x));
modParser.addAdditionalConstraint(@(x) mustBeGreaterThanOrEqual(x,1),"NB");

modParser.parse(params);

params = modParser.Results;


% construct the signals

% we will use an intensity of 1 for now as we can scale that up after the
% fact.
signals = {Coherent.pauliCoherentState(1,1,1);... %H
    Coherent.pauliCoherentState(1,1,2);... %V
    Coherent.pauliCoherentState(1,2,1);... %D
    Coherent.pauliCoherentState(1,2,2)}; %A

% build the (sub) isometry transition matrix that represents the channel
% and Bob's measurement except for the dark counts which must be handled
% later.
[transMat,detectorMat] = fourstateLinearOpticsSetup(params.transmittance,params.misalignmentAngle,params.birefringenceAngle,params.detectorEfficiency,params.probsB);
debugInfo.storeInfo("transMat",transMat);
debugInfo.storeInfo("detectorMat",detectorMat);

newParams.detectorMat = detectorMat;
%Calcualte subspace estimation from detector matrix for each n<= NB
vecSubspaceBnds = zeros(params.NB,1);
for indexn = 1:params.NB
    [~,cn,~,~] = detectorDecompFlag(detectorMat,indexn,true);
    vecSubspaceBnds(indexn) = cn;
end
newParams.subspaceBnd = vecSubspaceBnds;

%Calculate the conditional probabilities for the click patterns
probEachDetectorClickCon = zeros(numel(signals),size(transMat,1),numel(params.decoys));
expectationsCon = zeros(numel(signals),2^size(transMat,1),numel(params.decoys));

for index = 1:numel(params.decoys)
    %scale the signal states for the intensity
    signalsDecoy = cellfun(@(x) sqrt(params.decoys{index})*x,signals,"UniformOutput",false);
    [expectationsCon(:,:,index), probEachDetectorClickCon(:,:,index)] = simulateChannel(signalsDecoy,transMat,params.darkCountRate);
end
debugInfo.storeInfo("probEachDetectorClickCon",probEachDetectorClickCon);

%Expectations conditioned on test and decoy
newParams.expectationsConTestConDecoy = expectationsCon(3:4,:,:);

%Expectations conditioned on gen
expectationsConGen =0;
for index = 1:numel(params.probDecoyConGen)
    expectationsConGen = expectationsConGen + params.probDecoyConGen(index)*expectationsCon(1:2,:,index);
end
newParams.expectationsConGen = expectationsConGen;

end

function [transMat,detectorMat] = fourstateLinearOpticsSetup(eta,misalignmentAngle,birefringenceAngle,detectorEfficiency,probsB)

    % probDist = (pzB, pxB)
    % detectorEfficiency = (deteffH, ..., deteffA), 4x1 vector
    
    %% construct channel transition marix
    %loss/transmittance
    channelMat = Coherent.copyChannel(Coherent.transmittanceChannel(eta),2);
    
    %misalignment rotation (rotation around Y)
    channelMat = Coherent.rotateStateZXY(misalignmentAngle,[0,0,1],"angleOnBlochSphere",false)*channelMat;
    
    %birefringence rotation (rotation around Z)
    channelMat = Coherent.rotateStateZXY(birefringenceAngle,[1,0,0],"angleOnBlochSphere",false)*channelMat;
    
    %% Build up Bob's detector transition matrix
    
    % Bob applies a beam splitter to send signals to each detector basis setup
    detectorMat = Coherent.copyChannel(Coherent.singleInputMultiBeamSpliter(probsB),...
        2,"weaveCopies",true);
    
    % Bob applies a rotation to convert A and D back to H and V for easier
    % measurement
    detectorMat = blkdiag(pauliBasis(1,false)',pauliBasis(2,false)')*detectorMat;
    
    % We have to handle dark counts after we get the click probabilities for
    % each detector, so no changes here.
    
    % Each detector has a different efficiency so we can't pull it right to the
    % start.
    detectorMat = Coherent.transmittanceChannel(detectorEfficiency)*detectorMat;
    
    transMat = detectorMat*channelMat;
end

function probDetectorClickCon = applyDarkCounts(probDetectorClickCon,darkCountRate)
    probDetectorClickCon = 1-(1-probDetectorClickCon)*(1-darkCountRate);
end
    
function [probDetectorClickPatternCon, probEachDetectorClicksCon] = simulateChannel(signals,transMat,darkCountRate)
        %Construct the independent detector click probabilities for each signal
        probEachDetectorClicksCon = detectorClickProbabilities(signals,transMat);
        %simulate the effects of dark counts
        probEachDetectorClicksCon  = applyDarkCounts(probEachDetectorClicksCon,darkCountRate);
        %Construct all combinations of detector firing patterns from the
        %independent detectors.
        probDetectorClickPatternCon = detectorClickPatterns(probEachDetectorClicksCon);
end



function probDetectorClickCon = detectorClickProbabilities(signals,transMat)

    probDetectorClickCon = zeros(numel(signals),size(transMat,1));
    
    for index = 1:numel(signals)
        bobsSignal = transMat*signals{index};
        probDetectorClickCon(index,:) = 1-Coherent.fockCoherentProb(zeros(size(transMat,1),1),bobsSignal,"combineModes",false);
    end
end


function probDetectorClickPatternCon = detectorClickPatterns(probClickCon)
    %Because we have coherent states, each detector acts independently. This
    %function takes the independent results from each detector and computes the
    %probabilities of each click pattern outcome.
    numSignals = size(probClickCon,1);
    numDetectors = size(probClickCon,2);
    
    probDetectorClickPatternCon = zeros(numSignals,2^numDetectors);
    
    sizeDetectorPatterns = 2*ones(1,numDetectors);
    clickProbSwitch = @(click, clickProb) (click==0).*(1-clickProb) + (click~=0).*clickProb;
    for signalIndex = 1:numSignals
        for indexPat = 1:2^numDetectors
            patternVec = ind2subPlus(sizeDetectorPatterns,indexPat)-1;
            probDetectorClickPatternCon(signalIndex,indexPat) = prod(clickProbSwitch(patternVec,probClickCon(signalIndex,:)));
        end
    end
end

%% validation function helpers

function mustBeProbDistCell(input)
    mustBeProbDist([input{:}])
end