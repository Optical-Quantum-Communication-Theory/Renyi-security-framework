function [keyRate, modParser] = RenyiDecoyBB84PassiveKeyRateFunc(params,options,mathSolverFunc,debugInfo)

arguments
    params (1,1) struct
    options (1,1) struct
    mathSolverFunc (1,1) function_handle
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.addOptionalParam("photonCutOff",10, @(x) mustBePositive(x));
optionsParser.addAdditionalConstraint(@(x) mustBeInteger(x), "photonCutOff");
optionsParser.addAdditionalConstraint(@(x) mustBeFinite(x), "photonCutOff");

optionsParser.parse(options);
options = optionsParser.Results;


%% modParser
modParser = moduleParser(mfilename);

modParser.addRequiredParam("observablesJointTest",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("expectationsConTestConDecoy",@(x) eachRowMustBeAProbDist(x));

modParser.addRequiredParam("observablesJointGen",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("expectationsConGen",@(x) eachRowMustBeAProbDist(x));
modParser.addRequiredParam("decoys",@(x) allCells(x,@(y) y>=0));
modParser.addRequiredParam("deltaDecoys",@(x) allCells(x,@(y) y>=0));
modParser.addRequiredParam("decoyProbs",@mustBeProbDistCell); %,must sum to 1. 

modParser.addRequiredParam("probSignalConTest",@(x)mustBeProbDist(x));
modParser.addRequiredParam("probSignalConGen",@(x)mustBeProbDist(x));

modParser.addRequiredParam("probDecoyConGen",@(x)mustBeProbDist(x));

modParser.addRequiredParam("BobMultiPOVM",@(x) ishermitian(x));
modParser.addRequiredParam("BobSubspaceProj",@(x) ishermitian(x));
modParser.addRequiredParam("subspaceBnd", @(x) mustBeInRange(x, 0, 1));
modParser.addOptionalParam("NB",1,@(x) mustBeInteger(x));
modParser.addAdditionalConstraint(@(x) mustBeGreaterThanOrEqual(x,1),"NB");

modParser.addRequiredParam("POVMAGen",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("POVMATest",@(x) allCells(x,@ishermitian));

modParser.addRequiredParam("krausOps", @(x) allCells(x, @isCPTNIKrausOps));
modParser.addRequiredParam("keyProj", @(x) allCellsMust(x, @mustBeAKeyProj));

modParser.addRequiredParam("dimA",@mustBeInteger);
modParser.addRequiredParam("dimAPrime",@mustBeInteger);
modParser.addRequiredParam("dimB", @mustBeInteger);
modParser.addAdditionalConstraint(@mustBePositive,"dimA");
modParser.addAdditionalConstraint(@mustBePositive,"dimAPrime");
modParser.addAdditionalConstraint(@mustBePositive,"dimB");
modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJointTest","dimA","dimB"]);
modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJointGen","dimA","dimB"]);

modParser.addRequiredParam("announcementsA")
modParser.addRequiredParam("announcementsB")
modParser.addRequiredParam("keyMap",@(x)mustBeA(x,"KeyMapElement"))
% modParser.addAdditionalConstraint(@mustBeSizedLikeAnnouncements,["expectationsJointTestDecoy","announcementsA","announcementsB"])
% Change to add check that each page is the correct size
% modParser.addAdditionalConstraint(@mustBeSizedLikeAnnouncements,["expectationsJointGen","announcementsA","announcementsB"])

modParser.addRequiredParam("fEC", @(x) mustBeGreaterThanOrEqual(x,1));

modParser.addRequiredParam("probTest",@(x) x>=0);
% modParser.addRequiredParam("renyiAlpha",@(x)mustBeInRange(x,1,2,"exclude-lower"));
modParser.addRequiredParam("logrenyiAlpha",@(x)mustBeInRange(x,-10,0,"exclude-lower"));
modParser.addOptionalParam("epsilon", @(x) mustBeInRange(x, 0, 1));
modParser.addRequiredParam("Ntot",@(x) x>0);

modParser.addRequiredParam("stateTest",@(x) allCells(x,@isDensityOperator));
modParser.addRequiredParam("stateGen",@(x) allCells(x,@isDensityOperator));

modParser.addRequiredParam("blockPhotonNum", @(x) mustBeInteger(x));
modParser.addAdditionalConstraint(@(x) mustBeNonnegative(x),"blockPhotonNum");

modParser.addOptionalParam("blockDimsAPrime", nan, @(x) allCells(x, @isBlockDimsWellFormated));
modParser.addOptionalParam("blockDimsB", nan, @isBlockDimsWellFormated);
modParser.addAdditionalConstraint(@(x,y) eachCellBlockDimsMustMatch(x,y), ["blockDimsAPrime","dimAPrime"]);
modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsB","dimB"]);
modParser.addAdditionalConstraint(@(blockDimsAPrime,blockDimsB) ~xor(isequaln(blockDimsAPrime,nan),isequaln(blockDimsB,nan)),["blockDimsAPrime","blockDimsB"]);


modParser.parse(params);

params = modParser.Results;

%% simple setup
debugMathSolver = debugInfo.addLeaves("mathSolver");
mathSolverInput = struct();

%% Apply squashing
squashingMap = BB84PassivePostProcessingSquashingMap();

%test
squashedConExpTest = pagemtimes(params.expectationsConTestConDecoy,squashingMap);
expectationsJointTestConDecoy = pagemtimes(diag(params.probSignalConTest),squashedConExpTest);

%gen
squashedConExpGen = params.expectationsConGen*squashingMap;
expectationsJointGen = diag(params.probSignalConGen)*squashedConExpGen;

%create operators for squashing constraints in Gen and Test rounds (both
%can be used)
[squashingOpsTest,squashingOpsGen,lowerBndSubspaceTest,lowerBndSubspaceGen,upperBndSubspaceTest,upperBndSubspaceGen] = ...
    createSquashingCons(params.POVMATest,params.POVMAGen,params.probSignalConTest, ...
    params.probSignalConGen,params.subspaceBnd,params.BobMultiPOVM,params.NB);

%Store squashing constraints in Test rounds
mathSolverInput.squashingConstraintsTest = arrayfun(@(x)...
    InequalityConstraint(squashingOpsTest{x},lowerBndSubspaceTest(x), upperBndSubspaceTest(x)),1:numel(squashingOpsTest));

%Store squashing constraints in Gen rounds
mathSolverInput.squashingConstraintsGen = arrayfun(@(x)...
    InequalityConstraint(squashingOpsGen{x},lowerBndSubspaceGen(x), upperBndSubspaceGen(x)),1:numel(squashingOpsGen));

%% Error correction

deltaLeak = errorCorrectionCost(params.announcementsA,params.announcementsB,...
    expectationsJointGen,params.keyMap,params.fEC);
debugInfo.storeInfo("deltaLeak",deltaLeak);


%% translate for the math solver
%now we have all the parts we need to get a key rate from the a math
%solver, but we need to put it into a form it can understand.


%Convert decoys and delta decoys into array
decoysMat = cell2mat(params.decoys);
deltaDecoysMat = cell2mat(params.deltaDecoys);

%Bound photon number distribution in test rounds
probDistPhotonConMuTestLower = zeros(options.photonCutOff+1,numel(params.decoys));
probDistPhotonConMuTestUpper = zeros(options.photonCutOff+1,numel(params.decoys));

for indexDecoy = 1:numel(params.decoys)
    [probDistPhotonConMuTestLower(:,indexDecoy),probDistPhotonConMuTestUpper(:,indexDecoy)] = ...
        minmaxPoisVec(decoysMat(indexDecoy),deltaDecoysMat(indexDecoy),0:options.photonCutOff);
end
mathSolverInput.probDistPhotonConMuTestLower = probDistPhotonConMuTestLower;
mathSolverInput.probDistPhotonConMuTestUpper = probDistPhotonConMuTestUpper;

%first we give it the kraus operators for the G map and the projection
%operators for the key map (Z).

mathSolverInput.krausOps = params.krausOps;
mathSolverInput.keyProj = params.keyProj;

%Recast logAlpha to alpha
params.renyiAlpha = 1 + 10^(params.logrenyiAlpha);

mathSolverInput.renyiAlpha = params.renyiAlpha;
mathSolverInput.dimAPrime = params.dimAPrime;

mathSolverInput.stateTest = params.stateTest;
mathSolverInput.stateGen = params.stateGen;

%Bounds on photon distribution in generation round
photonDistributionGenUpper = 0;
photonDistributionGenLower = 0;

for index = 1:numel(params.probDecoyConGen)
    photonDistributionGenLower = photonDistributionGenLower + ...
        params.probDecoyConGen(index)*poisson(params.decoys{index} - params.deltaDecoys{index},params.blockPhotonNum);
    photonDistributionGenUpper = photonDistributionGenUpper + ...
        params.probDecoyConGen(index)*poisson(params.decoys{index} + params.deltaDecoys{index},params.blockPhotonNum);
end

mathSolverInput.photonDistributionGenLower = photonDistributionGenLower;
mathSolverInput.photonDistributionGenUpper = photonDistributionGenUpper;

%lower bound on sum_n p_mu(n)
[decoysGrid,countGrid] = meshgrid(cell2mat(params.decoys) + cell2mat(params.deltaDecoys),0:options.photonCutOff);
probDistPhotonConMuLowerTest = poisson(decoysGrid,countGrid);
probRemaining = sum(probDistPhotonConMuLowerTest,1);
mathSolverInput.probRemaining = min(probRemaining,1);

mathSolverInput.decoyTest = params.decoys;
mathSolverInput.probDecoyConTest = cell2mat(params.decoyProbs);
mathSolverInput.probSignalConTest = params.probSignalConTest;

mathSolverInput.blockPhotonNum = params.blockPhotonNum;

%reshape expectations
numObs = numel(params.observablesJointTest); %number of observations
numInt = size(params.decoys,2); %number of intensities in test rounds

reshapedExpectationsJointTestConDecoy = zeros(numObs,numInt);
for indexDecoy = 1:numInt
    reshapedExpectationsJointTestConDecoy(:,indexDecoy) = reshape(expectationsJointTestConDecoy(:,:,indexDecoy),[],1);
end

mathSolverInput.testConstraints = arrayfun(@(x)...
    EqualityConstraintDecoy(params.observablesJointTest{x},reshapedExpectationsJointTestConDecoy(x,:)),1:numObs);

mathSolverInput.probTest = params.probTest;

% if block diag information was given, then pass it to the solver.
if ~isequaln(params.blockDimsAPrime,nan)
    mathSolverInput.blockDimsAPrime = params.blockDimsAPrime;
    mathSolverInput.blockDimsB = params.blockDimsB;
end

% now we call the math solver function on the formulated inputs, with
% it's options.
[relEnt,~] = mathSolverFunc(mathSolverInput,debugMathSolver);

%store the key rate (even if negative)
keyRateFixed = finiteKeyRate(relEnt, deltaLeak, params.renyiAlpha, params.epsilon, params.Ntot, 1-params.probTest, debugInfo);
debugInfo.storeInfo("keyRateFixed",keyRateFixed);

%extract QES from debugInfo if toggled
if isfield(debugInfo.leaves.mathSolver.info,"qesCell")
    %extract qes
    qesCell = debugInfo.leaves.mathSolver.info.qesCell;

    %test constraints
    testCons = mathSolverInput.testConstraints;

    %testCons*p(mu|test)
    pMuTestcons = reshape(cat(1,testCons.vector)*diag(cell2mat(params.decoyProbs)),[],1);

    %Lower bound on relative entropy from QES
    QESRelEnt = qesCell{1}.'*([params.probTest*pMuTestcons;1-params.probTest]) + qesCell{2};

    %Key rate from QES
    QESRate = finiteKeyRate(QESRelEnt, deltaLeak, params.renyiAlpha, params.epsilon, params.Ntot, 1-params.probTest, debugInfo);

    %Store QES Rate in debuginfo
    debugInfo.storeInfo("QESRate",QESRate);
    debugInfo.storeInfo("qesCell",qesCell);

    %set key rate to QES rate
    keyRate = QESRate;
else
    keyRate = keyRateFixed;
end

if options.verboseLevel>=1
    %ensure that we cut off at 0 when we display this for the user.
    fprintf("Key rate: %e\n",max(keyRateFixed,0));
    if isfield(debugInfo.leaves.mathSolver.info,"qesCell")
        fprintf("QES rate = %e and difference from fixed length = %e \n", max(QESRate,0), keyRateFixed - QESRate);
    end
end

%set the linearization estimate key rate as well for debuging
if isfield(debugMathSolver.info,"relEntStep2Linearization")
    keyRateStep2Linearization = debugMathSolver.info.relEntStep2Linearization - deltaLeak; 
    debugInfo.storeInfo("keyRateRelEntStep2Linearization",keyRateStep2Linearization)

    if options.verboseLevel>=2
        fprintf("Key rate using step 2 linearization intial value: %e\n",max(keyRateStep2Linearization,0))
    end
end

%%%%%%%%%%%  FINITE Key Rate Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finite-size key rate function combining results from 2 works
% 
% Insert eq. (20) into Theorem 1 from: 
% Finite-size analysis of prepare-and-measure and decoy-state QKD via entropy accumulation 
% http://arxiv.org/abs/2406.10198
%
% and then use the lower bound on the Renyi entropy from Lemma 10 of:
% Generalized R\'enyi entropy accumulation theorem and generalized quantum
% probability estimation
% http://arxiv.org/abs/2405.05912

    function [keyRate, debugInfo] = finiteKeyRate(relEnt, deltaLeak, renyiAlpha, epsilon, Ntotal, pGen, debugInfo)
    %computes the finite size keyrate.

    %Error correction including error verification
    ECLeakage = pGen*deltaLeak + ceil(log2(1/epsilon.EC))/Ntotal;
    
    %Privacy Amplification
    privacyAmplification = renyiAlpha/(renyiAlpha-1)*log2(1/epsilon.PA)/Ntotal;

    keyRate = relEnt - ECLeakage - privacyAmplification + 2/Ntotal;
end

end

%% helper functions
%Poisson disctribution
function prob = poisson(mu,m)
    prob = mu.^m.*exp(-mu)./factorial(m);
end

%creates operators and bounds for squashing constraints
function [squashingOpsTest,squashingOpsGen,lowerBndSubspaceTest,lowerBndSubspaceGen,upperBndSubspaceTest,upperBndSubspaceGen] = ...
    createSquashingCons(POVMATest,POVMAGen,probSignalConTest,probSignalConGen,subspaceBnd,BobMultiPOVM,NB)
    %create operators for squashing constraints in Gen and Test rounds (both
    %can be used)
    squashingOpsGen = cell(1,numel(POVMAGen)*NB);
    squashingOpsTest = cell(1,numel(POVMATest)*NB);
    
    %lower bounds in Gen and test
    lowerBndSubspaceGen = zeros(1,numel(POVMAGen)*NB);
    lowerBndSubspaceTest = zeros(1,numel(POVMATest)*NB);
    
    %lower bounds in Gen and test
    upperBndSubspaceGen = zeros(1,numel(POVMAGen)*NB);
    upperBndSubspaceTest = zeros(1,numel(POVMATest)*NB);
    
    %dimension of Bob's system
    dimB = size(BobMultiPOVM,1);
    
    %Number of Alice's POVM elements used for bounds in gen and test rounds    
    numPovmGen = numel(POVMAGen);
    numPovmTest = numel(POVMATest);
    
    %we iterate over all subspace with <= NB photons, i.e. n = (0,1),
    %n = (0,1,2), n = (0,1,2,3), ...
    for indexn = 1:NB
        %calculate dimension of current subspace <= indexn
        subspacedim = (indexn+1)*(indexn+2)/2;
        
        %vector with all ones in current subspace and 0s else 
        vector = zeros(1,dimB);
        vector(1:subspacedim) = 1;

        %projector onto subspace defined via vector
        BobSubspaceProj = diag(vector);
        
        %Bob's operator involved in squashing is sum of multiclick POVM and
        %subspace projector
        BobSquashingOp = BobMultiPOVM + subspaceBnd(indexn)*BobSubspaceProj;
        
        %iterate over Alice's POVM elements in gen rounds
        for indexA = 1:numPovmGen
            %Define joint operator for squashing
            squashingOpsGen{(indexn-1)*numPovmGen+indexA} = kron(POVMAGen{indexA},BobSquashingOp);
            %lower bound on subspace weight given by 
            % subspacebnd* prob(Alice sends state)
            lowerBndSubspaceGen((indexn-1)*numPovmGen+indexA) = subspaceBnd(indexn) *probSignalConGen(indexA);
            %upper bound on subspace weight given by prob(Alice sends state)
            upperBndSubspaceGen((indexn-1)*numPovmGen+indexA) = probSignalConGen(indexA);
        end
        
        for indexA = 1:numPovmTest
            %Define joint operator for squashing
            squashingOpsTest{(indexn-1)*numPovmTest+indexA} = kron(POVMATest{indexA},BobSquashingOp);
            %lowerbound on subspace weight given by 
            % subspacebnd* prob(Alice sends state)
            lowerBndSubspaceTest((indexn-1)*numPovmTest+indexA) = subspaceBnd(indexn) *probSignalConTest(indexA);
            %upper bound on subspace weight given by prob(Alice sends state)
            upperBndSubspaceTest((indexn-1)*numPovmTest+indexA) = probSignalConTest(indexA);
        end
    end
end


%% validation functions
function mustBeProbDistCell(input)
    mustBeProbDist([input{:}])
end


function eachRowMustBeAProbDist(expectationsConditional)

    % get the dimensions of the conditional expectations. Then based on that
    % pick a strategy to handle it
    dimExpCon = size(expectationsConditional);
    
    errorID ="BasicBB84WCPDecoyKeyRateFunc:InvalidRowsAreNotProbDists";
    errorTXT = "A row in the conditional distribution is not a valid probability distribution.";
    
    if numel(dimExpCon) == 2 % Matlab's minimum number of dimensions is 2.
        % The array is 2d and the slicing is easy
        for index = 1:dimExpCon(1)
            if~isProbDist(expectationsConditional(index,:))
               throwAsCaller(MException(errorID,errorTXT));
            end
        end
    else
        % We have some tricky slicing to do for 3 plus dimensions.
        % We need to index the first dimension and the combination of
        % dimensions 3 and up. The second dimension will just use :.
        maskedDims = [dimExpCon(1),prod(dimExpCon(3:end))];
    
        for index = 1:prod(maskedDims)
            vecIndex = ind2subPlus(maskedDims,index);
            if ~isProbDist(expectationsConditional(vecIndex(1),:,vecIndex(2)))
                throwAsCaller(MException(errorID,errorTXT));
            end
        end
    end
end

function eachCellBlockDimsMustMatch(blockdimsCell,dimsCell)
    if numel(blockdimsCell) ~= numel(dimsCell)
        errorID ="RenyiDecoyBB84KeyRateFunc:InvalidNumberOfBlocksDims";
        errorTXT = "The number of blockdimensions does not match the number of dimensions.";
        throwAsCaller(MException(errorID,errorTXT));
    else    
        for index = 1:numel(dimsCell)
            blockDimsMustMatch(blockdimsCell{index},dimsCell(index));
        end
    end
end