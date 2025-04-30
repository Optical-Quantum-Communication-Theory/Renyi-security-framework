function [newParams, modParser]= RenyiB84LossyChannelFunc(params,options, debugInfo)
% BasicBB84LossyChannelFunc A channel function for qubit BB84 with loss.
% Here, Schmidt decomposition was used to shrink Alice from a 4d space to a
% 2d space.
% 
% Input parameters:
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * observablesJoint: The joint observables from Alice and Bob's
%   measurements. The observables must be hermitian and each must be the
%   size dimA*dimB by dimA*dimB. The observables assume the spaces are
%   ordered A \otimes B. They also should be positive semi-definite and
%   should sum to identity, but this is hard to check because of machine
%   precision issues.
% * transmittance (1): the transmissivity of the quantum channel. Must be
%   between 0 and 1 inclusive.
% * depolarization (0): The amount of depolarization applied to the signal
%   Alice sends to Bob. At maximum depolarization (depolariztion =1) a pure
%   qubit state is converted to a maximally mixed state. Depolarization
%   should be between 0 and 1.
% * misalignmentAngle (0): Physical angle of misalignment between Alice and
%   Bob's measurements around Y axix. This angle is measured as the
%   physical rotation of the device (period 2pi). Although calculations are
%   done on the Bloch sphere, angles should not be given in that form
%   (period 4pi).
% Output parameters:
% * expectationsJoint: The joint expectations for Alice and Bob's
%   measurement of the signals. Computed by taking the
%   observablesJoint and applying them to a simulated rhoAB.
% Options:
% * none
% DebugInfo:
% * rhoAB: Alice and Bob's shared density matrix after the channel has
%   acted on it. Usefull for checking the channel has been applied
%   correctly.
%
% Reviewed by Devashish Tupkary 2023/09/18
% See also QKDChannelModule
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