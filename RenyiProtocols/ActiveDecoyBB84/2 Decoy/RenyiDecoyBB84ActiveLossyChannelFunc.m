function [newParams, modParser]= RenyiDecoyBB84ActiveLossyChannelFunc(params,options, debugInfo)
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
modParser.addRequiredParam("decoys", @(x) mustBeCellOf(x, 'numeric'));
modParser.addAdditionalConstraint(@(x) allCells(x,@isscalar),"decoys");
modParser.addAdditionalConstraint(@(x) allCells(x,@(y) y>=0),"decoys");

modParser.addRequiredParam("decoyProbs",@mustBeProbDistCell); %,must sum to 1. 

modParser.addRequiredParam("probTest", @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"probTest");

modParser.addOptionalParam("transmittance", 1, @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"transmittance");

modParser.addOptionalParam("detectorEfficiency", 1, @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"detectorEfficiency");

modParser.addOptionalParam("misalignmentAngle",0,@mustBeReal);
modParser.addAdditionalConstraint(@isscalar,"misalignmentAngle");

modParser.addOptionalParam("darkCountRate", 0, @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"darkCountRate");

modParser.parse(params);

params = modParser.Results;

%z basis only used for key generation
pz = 1 - params.probTest;


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

[transMatHV,transMatDA] = simpleBB84ActiveLinearOpticsSetup(params.transmittance,params.misalignmentAngle,params.detectorEfficiency);
debugInfo.storeInfo("transMat",{transMatHV,transMatDA});

%Calculate the conditional probabilities for the click patterns
expectationsCon = zeros(numel(signals),2*2^size(transMatHV,1),numel(params.decoys));

for index = 1:numel(params.decoys)
    %scale the signal states for the intensity
    signalsDecoy = cellfun(@(x) sqrt(params.decoys{index})*x,signals,"UniformOutput",false);
    expectationsCon(:,1:4,index) = pz *     simulateChannel(signalsDecoy,transMatHV,params.darkCountRate);
    expectationsCon(:,5:8,index) = (1-pz) * simulateChannel(signalsDecoy,transMatDA,params.darkCountRate);
end

newParams.expectationsConTestConDecoy = expectationsCon(3:4,:,:); 
newParams.expectationsConGen = expectationsCon(1:2,:,1);

end

function [transMatHV,transMatDA] = simpleBB84ActiveLinearOpticsSetup(transmittance,misalignmentAngle,detectorEfficiency)

%% construct channel transition marix
%loss/transmittance
channelMat = Coherent.copyChannel(Coherent.transmittanceChannel(transmittance),2);

%misalignment rotation
channelMat = Coherent.rotateStateZXY(misalignmentAngle,[0,0,1],"angleOnBlochSphere",false)*channelMat;


%% Build up Bob's detector transition matrix
% Each detector has the same efficiency so we can pull it right to the
% start.
detectorMatHV = Coherent.transmittanceChannel(detectorEfficiency);
detectorMatDA = detectorMatHV;

% Bob applies a rotation to convert A and D back to H and V for easier
% measurement
detectorMatDA = pauliBasis(2,false).'*detectorMatDA;

% We have to handle dark counts after we get the click probabilities for
% each detector, so no changes here.

transMatHV = detectorMatHV*channelMat;
transMatDA = detectorMatDA*channelMat;
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



function mustBeProbDistCell(input)
mustBeProbDist([input{:}])
end