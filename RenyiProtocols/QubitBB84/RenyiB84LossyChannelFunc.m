function [newParams, modParser]= RenyiB84LossyChannelFunc(params,options, debugInfo)
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
modParser.addRequiredParam("observablesJointTest",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("observablesJointGen",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("dimA",@(x) x==2);
modParser.addRequiredParam("dimAPrime",@(x) x==2);
modParser.addRequiredParam("dimB", @(x) x ==3); 

modParser.addOptionalParam("transmittance", 1, @(x) mustBeInRange(x,0,1));
modParser.addOptionalParam("depolarization",0, @(x) mustBeInRange(x,0,1));
modParser.addOptionalParam("misalignmentAngle",0,@(x) mustBeReal(x));

modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJointTest","dimA","dimB"]);
modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJointGen","dimA","dimB"]);

modParser.addRequiredParam("stateTest",@isDensityOperator);
modParser.addRequiredParam("stateGen",@isDensityOperator);

modParser.parse(params);

params = modParser.Results;

%% simple setup

newParams = struct();

%Observables conditioned on test
observablesJointTest = params.observablesJointTest;

%Observables conditioned on gen
observablesJointGen = params.observablesJointGen;

dimA = params.dimA;
dimAPrime = params.dimAPrime;


%% generate rho and apply quantum channel

% depolarization
depolChoiMat =  Qudit.depolarizationChoiMat(dimA,params.depolarization);

% rotation
rotationKrausOps = {Qudit.rotateStateZXY(params.misalignmentAngle,[0,0,1],"angleOnBlochSphere",false)};

% transmittance/loss
transmittanceChoiMat = Qudit.transmittanceChoiMat(params.transmittance,dimA); %last dimension is loss

channelChoi = PartialMap(depolChoiMat,rotationKrausOps,2,[dimA,dimAPrime]);
channelChoi = PartialMap(channelChoi,transmittanceChoiMat,2,[dimA,dimAPrime]);

rhoABTest = PartialMap(params.stateTest,channelChoi,2,[dimA,dimAPrime]);
rhoABGen = PartialMap(params.stateGen,channelChoi,2,[dimA,dimAPrime]);

debugInfo.storeInfo("rhoABTest",rhoABTest);
debugInfo.storeInfo("rhoABGen",rhoABGen);

%% compute expectations
% load observables from description to compute expectations
newParams.expectationsJointTest = calcExpectationsJoint(rhoABTest,observablesJointTest);
newParams.expectationsJointGen = calcExpectationsJoint(rhoABGen,observablesJointGen);
end

function expectationsJoint = calcExpectationsJoint(rhoAB,observablesJoint)
expectationsJoint = zeros(size(observablesJoint));
for index = 1:numel(observablesJoint)
    expectationsJoint(index) = real(trace(observablesJoint{index}'*rhoAB));
end
end