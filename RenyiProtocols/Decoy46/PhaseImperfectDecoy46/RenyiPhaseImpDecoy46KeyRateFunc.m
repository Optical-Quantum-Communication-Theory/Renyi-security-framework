function [keyRate, modParser] = RenyiPhaseImpDecoy46KeyRateFunc(params,options,mathSolverFunc,debugInfo)

arguments
    params (1,1) struct
    options (1,1) struct
    mathSolverFunc (1,1) function_handle
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
% optionsParser.addOptionalParam("photonCutOff",10, @(x) mustBePositive(x));
% optionsParser.addAdditionalConstraint(@(x) mustBeInteger(x), "photonCutOff");
% optionsParser.addAdditionalConstraint(@(x) mustBeFinite(x), "photonCutOff");

optionsParser.parse(options);
options = optionsParser.Results;


%% modParser
modParser = moduleParser(mfilename);

modParser.addRequiredParam("observablesJointTest",@(x) allPages(x,@ishermitian));
modParser.addRequiredParam("expectationsConTestConDecoy",@(x) eachRowMustBeAProbDist(x));

modParser.addRequiredParam("observablesJointGen",@(x) allPages(x,@ishermitian));
modParser.addRequiredParam("expectationsConGen",@(x) eachRowMustBeAProbDist(x));
modParser.addRequiredParam("decoys",@(x) allCells(x,@(y) y>=0));
modParser.addRequiredParam("decoyProbs",@mustBeProbDistCell); %,must sum to 1. 

modParser.addRequiredParam("probSignalConTest",@(x)mustBeProbDist(x));
modParser.addRequiredParam("probSignalConGen",@(x)mustBeProbDist(x));

modParser.addRequiredParam("probSignalConTestAndBlock",@(x) eachRowMustBeAProbDist(x.'));
modParser.addRequiredParam("probSignalConGenAndBlock",@(x) eachRowMustBeAProbDist(x.'));

modParser.addRequiredParam("probBlockCondTest",@(x) sum(x,2)<=1);
modParser.addRequiredParam("probBlockCondGen",@(x) sum(x,2)<=1); 

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
modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSameForAllPages,...
    ["observablesJointTest","dimA","dimB"]);
modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSameForAllPages,...
    ["observablesJointGen","dimA","dimB"]);

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

% correction terms due to phase imperfection
modParser.addRequiredParam("epsilonFull", @(x) mustBeInRange(x, 0, 1));
modParser.addRequiredParam("epsilonBlockCondTest", @(x) mustBeInRange(x, 0, 1));
modParser.addRequiredParam("epsilonProbTest", @(x) mustBeInRange(x, 0, 1));

%Block cutoff used for key rate
modParser.addRequiredParam("blockCutoff",@(x) mustBeInteger(x));
modParser.addAdditionalConstraint(@(x) mustBeGreaterThanOrEqual(x,1),"blockCutoff");

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
squashingMap = fourSixPostProcessingSquashingMap();

%test
squashedConExpTest = pagemtimes(params.expectationsConTestConDecoy,squashingMap);

% %keep page structure of expectations the same
% probSignalOnly = [1,1,0,0;0,0,1,1]*params.probSignalConTest ; %only works here!
% expectationsJointTestConDecoy = pagemtimes(diag(probSignalOnly),squashedConExpTest);

splitSquashedConExpTest = num2cell(squashedConExpTest, [1 2]);
expectationsJointTestConDecoy = diag(params.probSignalConTest)*vertcat(splitSquashedConExpTest{:});

%gen
squashedConExpGen = params.expectationsConGen*squashingMap;
expectationsJointGen = diag(params.probSignalConGen)*squashedConExpGen;

%create operators for squashing constraints in Gen and Test rounds (both
%can be used)
[squashingOpsTest,squashingOpsGen,lowerBndSubspaceTest,lowerBndSubspaceGen,upperBndSubspaceTest,upperBndSubspaceGen] = ...
    createSquashingCons(params.POVMATest,params.POVMAGen,params.probSignalConTestAndBlock, ...
    params.probSignalConGenAndBlock,params.subspaceBnd,params.BobMultiPOVM,params.NB);


%Store squashing constraints in Test rounds
for indexBlock = 1:params.blockCutoff+1
    squashingConstraintsTest(:,indexBlock) = arrayfun(@(x)...
    InequalityConstraint(squashingOpsTest{x,indexBlock},lowerBndSubspaceTest(x,indexBlock), upperBndSubspaceTest(x,indexBlock)),1:size(squashingOpsTest,1));
end
mathSolverInput.squashingConstraintsTest = squashingConstraintsTest;

%Store squashing constraints in Gen rounds
for indexBlock = 1:params.blockCutoff+1
    squashingConstraintsGen(:,indexBlock) = arrayfun(@(x)...
    InequalityConstraint(squashingOpsGen{x,indexBlock},lowerBndSubspaceGen(x,indexBlock), upperBndSubspaceGen(x,indexBlock)),1:size(squashingOpsGen,1));
end
mathSolverInput.squashingConstraintsGen = squashingConstraintsGen;

%% Error correction

deltaLeak = errorCorrectionCost(params.announcementsA,params.announcementsB,...
    expectationsJointGen,params.keyMap,params.fEC);
debugInfo.storeInfo("deltaLeak",deltaLeak);

%% Continuity bound correction from phase imperfection
continuityBnd = continuityCorr(1 + 10^(params.logrenyiAlpha),params.epsilonFull,2);

%% translate for the math solver
%now we have all the parts we need to get a key rate from the a math
%solver, but we need to put it into a form it can understand.
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

%Block distribution in generation and test round
mathSolverInput.probBlockCondGen = params.probBlockCondGen;
mathSolverInput.probBlockCondTest = params.probBlockCondTest;
mathSolverInput.probSignalConTest = params.probSignalConTest;
mathSolverInput.probSignalConTestAndBlock = params.probSignalConTestAndBlock;
mathSolverInput.blockPhotonNum = params.blockPhotonNum;

mathSolverInput.probTest = params.probTest;

%total probability of all blocks in test rounds
pTotTest = sum(params.probBlockCondTest,2);

%corrections due to phase imperfections
mathSolverInput.continuityBnd = continuityBnd;
mathSolverInput.epsilonBlockCondTest = params.epsilonBlockCondTest;
mathSolverInput.epsilonProbTest = params.epsilonProbTest;

%Construct constraints
numObsAlice = size(expectationsJointTestConDecoy,1); %number of observations from Alice
numObsBob = size(expectationsJointTestConDecoy,2); %number of observations from Bob

for indexA = 1:numObsAlice
    for indexB = 1:numObsBob
        testConstraints(indexA, indexB) = EqualityConstraintChoiBlock(params.observablesJointTest(indexA,indexB,:),...
            expectationsJointTestConDecoy(indexA,indexB));
        probRemaining(indexA,indexB) = params.probSignalConTest(indexA)*(1-pTotTest); 
    end
end
mathSolverInput.testConstraints = testConstraints(:);
mathSolverInput.probRemaining = probRemaining(:);

% if block diag information was given, then pass it to the solver.
if ~isequaln(params.blockDimsAPrime,nan)
    mathSolverInput.blockDimsAPrime = params.blockDimsAPrime;
    mathSolverInput.blockDimsB = params.blockDimsB;
end

% now we call the math solver function on the formulated inputs, with
% it's options.
[relEnt,~] = mathSolverFunc(mathSolverInput,debugMathSolver);

%subtract continuity bnd from relEnt
relEnt = relEnt - continuityBnd;

%store the key rate (even if negative)
keyRateFixed = finiteKeyRate(relEnt, deltaLeak, params.renyiAlpha, params.epsilon, params.Ntot, 1-params.probTest, debugInfo);
debugInfo.storeInfo("keyRateFixed",keyRateFixed);

%extract QES from debugInfo if toggled
if isfield(debugInfo.leaves.mathSolver.info,"qesCell")
    %extract qes
    qesCell = debugInfo.leaves.mathSolver.info.qesCell;

    %test constraints
    testCons = mathSolverInput.testConstraints;

    %testCons
    testConsProbs = reshape(cat(1,testCons.scalar),[],1);

    %Lower bound on relative entropy from QES
    QESRelEnt = qesCell{1}.'*([params.probTest*testConsProbs;1-params.probTest]) + qesCell{2};

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

%creates operators and bounds for squashing constraints
function [squashingOpsTest,squashingOpsGen,lowerBndSubspaceTest,lowerBndSubspaceGen,upperBndSubspaceTest,upperBndSubspaceGen] = ...
    createSquashingCons(POVMATest,POVMAGen,probSignalConTestAndBlock,probSignalConGenAndBlock,subspaceBnd,BobMultiPOVM,NB)
    %create operators for squashing constraints in Gen and Test rounds (both
    %can be used)
    squashingOpsGen = cell(size(POVMAGen,1)*NB,size(POVMAGen,2));
    squashingOpsTest = cell(size(POVMATest,1)*NB,size(POVMATest,2));
    
    %lower bounds in Gen and test
    lowerBndSubspaceGen = zeros(size(POVMAGen,1)*NB,size(POVMAGen,2));
    lowerBndSubspaceTest = zeros(size(POVMATest,1)*NB,size(POVMATest,2));
    
    %lower bounds in Gen and test
    upperBndSubspaceGen = zeros(size(POVMAGen,1)*NB,size(POVMAGen,2));
    upperBndSubspaceTest = zeros(size(POVMATest,1)*NB,size(POVMATest,2));
    
    %dimension of Bob's system
    dimB = size(BobMultiPOVM,1);
    
    %Number of Alice's POVM elements used for bounds in gen and test rounds    
    numPovmGen = size(POVMAGen,1);
    numPovmTest = size(POVMATest,1);
    numBlocks = size(POVMATest,2);
        
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
        
        %iterate over Alice's POVM elements in gen rounds and blocks
        for indexBlock = 1:numBlocks
            for indexA = 1:numPovmGen
                %Define joint operator for squashing
                squashingOpsGen{(indexn-1)*numPovmGen+indexA,indexBlock} = kron(POVMAGen{indexA,indexBlock},BobSquashingOp);
                %lower bound on subspace weight given by 
                % subspacebnd* prob(Alice sends state)
                lowerBndSubspaceGen((indexn-1)*numPovmGen+indexA,indexBlock) = subspaceBnd(indexn) *probSignalConGenAndBlock(indexA,indexBlock);
                %upper bound on subspace weight given by prob(Alice sends state)
                upperBndSubspaceGen((indexn-1)*numPovmGen+indexA,indexBlock) = probSignalConGenAndBlock(indexA,indexBlock);
            end
        end
        
        for indexBlock = 1:numBlocks
            for indexA = 1:numPovmTest
                %Define joint operator for squashing
                squashingOpsTest{(indexn-1)*numPovmTest+indexA,indexBlock} = kron(POVMATest{indexA,indexBlock},BobSquashingOp);
                %lowerbound on subspace weight given by 
                % subspacebnd* prob(Alice sends state)
                lowerBndSubspaceTest((indexn-1)*numPovmTest+indexA,indexBlock) = subspaceBnd(indexn)*probSignalConTestAndBlock(indexA,indexBlock);
                %upper bound on subspace weight given by prob(Alice sends state)
                upperBndSubspaceTest((indexn-1)*numPovmTest+indexA,indexBlock) = probSignalConTestAndBlock(indexA,indexBlock);
            end
        end
    end
end

%helper function for continuity bound
function corr = continuityCorr(alpha,epsilonFull,dS)
    %Different cases for correction dpending on proof technique see Bluhm
    %et. al.
    %Convert alpha and epsilon to symbolic values
    alphaSym = sym(alpha);
    epsilonFullSym = sym(epsilonFull);

    % First term
    term1 = log2(1 + epsilonFullSym) + (1 / (alphaSym - 1)) * ...
        log2(1 + epsilonFullSym * dS^(alphaSym - 1) - (epsilonFullSym^alphaSym) / ((1 + epsilonFullSym)^(alphaSym - 1)));
    
    % Second term
    term2 = (alphaSym / (alphaSym - 1)) * log2(1 + epsilonFullSym * dS^((alphaSym-1) / alphaSym));
    
    % Third term
    term3 = log2(1 + epsilonFullSym) + (alphaSym / (alphaSym - 1)) * ...
        log2(1 + epsilonFullSym * dS^((alphaSym - 1)/alphaSym) - (epsilonFullSym^(2 - 1/alphaSym)) / ((1 + epsilonFullSym)^((alphaSym - 1) / alphaSym)));
    
    % Compute the minimum of the terms
    corr = min(double([term1, term2, term3]),[],"all");
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

%validation function applying checks on each page of cell array
function allPages(cellMat,func)
    numPages = size(cellMat,3);
    for index = 1:numPages
        allCells(cellMat(:,:,index),func);
    end
end

%validation function checking dimensions across pages
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