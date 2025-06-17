classdef RenyiTildeDownDecoyPA

    methods (Static)

        %% Frank Wolfe functions
        function val = funcFW(vecFW,testProb,alpha,perturbation,perturbationAff,dimAPrime,...
                genStateAAPrime,keyProj,krausOps,probBlocks)
            arguments
                vecFW (:,1) double
                testProb (1,1) double
                alpha (1,1) double
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
                perturbationAff (1,1) double {mustBeInRange(perturbationAff,0,1)}
                dimAPrime (1,:) double {mustBeInteger,mustBePositive}
                genStateAAPrime (1,:) cell {mustBeNonempty,mustBeCellOf(genStateAAPrime,"double")}% cell of (:,:) double
                keyProj (1,:) cell {mustBeNonempty,mustBeCellOf(keyProj,"cell")}% cell (:,1) cell of (:,:) double
                krausOps (1,:) cell {mustBeNonempty,mustBeCellOf(krausOps,"cell")}% cell of (:,1) cell of (:,:) double
                probBlocks (1,:) double {mustBeNonnegative,mustBeNonempty} = 1
            end

            %Dimensions of the Choi matrices
            dimsAB = cellfun(@(x) size(x{1},2), krausOps);
            dimAAPrime = cellfun(@(x) size(x,1), genStateAAPrime);
            dimA = dimAAPrime./dimAPrime;
            dimsAPrimeB = dimsAB./dimA .* dimAPrime;


            % parse classical relative entropy upper bound and choi
            % matrices
            [phi,choiEveAPrimeB] = unpackVecFW(vecFW, dimsAPrimeB);

            val = privacyAmpSplit(phi,choiEveAPrimeB,genStateAAPrime, ...
                keyProj,krausOps,testProb,alpha,dimAPrime,perturbation,perturbationAff,probBlocks);
        end


        function gradValVec = gradFuncFW(vecFW,testProb,alpha,perturbation,perturbationAff,dimAPrime,genStateAAPrime,keyProj,krausOps,probBlocks)
            arguments
                vecFW (:,1) double
                testProb (1,1) double
                alpha (1,1) double
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
                perturbationAff (1,1) double {mustBeInRange(perturbationAff,0,1)}
                dimAPrime (1,:) double {mustBeInteger,mustBePositive}
                genStateAAPrime (1,:) cell {mustBeNonempty,mustBeCellOf(genStateAAPrime,"double")}% cell of (:,:) double
                keyProj (1,:) cell {mustBeNonempty,mustBeCellOf(keyProj,"cell")}% cell (:,1) cell of (:,:) double
                krausOps (1,:) cell {mustBeNonempty,mustBeCellOf(krausOps,"cell")}% cell of (:,1) cell of (:,:) double
                probBlocks (1,:) double {mustBeNonnegative,mustBeNonempty} = 1
            end

            %Dimensions of the Choi matrices
            dimsAB = cellfun(@(x) size(x{1},2), krausOps);
            dimAAPrime = cellfun(@(x) size(x,1), genStateAAPrime);
            dimA = dimAAPrime./dimAPrime;
            dimsAPrimeB = dimsAB./dimA .* dimAPrime;


            % parse classical relative entropy upper bound and choi
            % matrices
            [phi,choiEveAPrimeB] = unpackVecFW(vecFW, dimsAPrimeB);

            [gradByChoi,gradByPhi] = gradPrivacyAmpSplit(phi, ...
                choiEveAPrimeB,genStateAAPrime,keyProj,krausOps,testProb,...
                alpha,dimAPrime,perturbation,perturbationAff,probBlocks);

            %Reshape each gradient in the cell gradbyChoi to be linearized
            gradByChoi = cell2mat(cellfun(@(x) x(:),gradByChoi, "uniformoutput", false));

            %stack all gradients together in the same shape as vecFW
            gradValVec = [gradByChoi; gradByPhi];
        end


        function [deltaVecFW,exitFlag,statusMessage] = subproblemUniqueDecoy(vecFW,vecFWGrad,...
                testStatesAAPrime,genStatesAAprime,dimsAPrime,testProb,testCons,squashingConsTest,squashingConsGen,probDistPhotonConMuTest,...
                decoyProbs,probSignalConTest,blockPhotonNum,options,feasibleMode)
            arguments
                vecFW (:,1) double
                vecFWGrad (:,1) double
                testStatesAAPrime (1,:) cell {mustBeNonempty,mustBeCellOf(testStatesAAPrime,"double")}% cell of (:,:) double
                genStatesAAprime (1,:) cell {mustBeNonempty,mustBeCellOf(genStatesAAprime,"double")}
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
                testProb (1,1) double {mustBeInRange(testProb, 0, 1)}
                testCons (:,1) EqualityConstraintDecoy
                squashingConsTest (:,1) InequalityConstraint
                squashingConsGen (:,1) InequalityConstraint
                probDistPhotonConMuTest (:,:) double {mustBePositive}
                decoyProbs (:,1) double {mustBeProbDist}
                probSignalConTest (:,1) double {mustBeProbDist}
                blockPhotonNum (:,1) uint64
                options (1,1) struct
                feasibleMode (1,1) logical = false;
            end

            %%%% !!!!!Testcons need to reflect block structure!!!!!
            
            %Linear constraint tolerance
            linConTol = options.linearConstraintTolerance;

            %Dimensions of the Choi matrices
            dimsAB = size(testCons(1).operator,1);
            dimsAAPrime = cellfun(@(x) size(x,1), testStatesAAPrime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;


            %read of various sizes
            photonCutoff = size(probDistPhotonConMuTest,1)-1; %Photon Cutoff
            numObs = numel(testCons); %Number of observation
            numTestInt = numel(decoyProbs); %Number of intensities in test rounds
            numBlocks = numel(dimsAPrimeB); %Number of Blocks

            % unpack all the values from vecFW and vecFWGrad
            [phi,choiEveAPrimeB]         = unpackVecFW(vecFW,    dimsAPrimeB);
            [gradPhi,gradChoiEveAPrimeB] = unpackVecFW(vecFWGrad,dimsAPrimeB);

            %% CVX
            cvx_begin
            cvx_solver(convertStringsToChars(options.cvxSolver));
            cvx_precision(convertStringsToChars(options.cvxPrecision));
            cvx_quiet(options.verboseLevel<2);


            %% variables
            deltaPhi = variable('deltaPhi(1,1)');
            YieldMat = variable(sprintf('YieldMat(%d,%d)',numObs,photonCutoff+1), 'nonnegative');
            rhoCTestMat = variable(sprintf('rhoCTestMat(%d,%d)',numObs,numTestInt), 'nonnegative');

            if options.blockDiagonal
                %define each block and its corresponding string as elements
                %of cell array
                deltaChoiBlocks = cell(1,numBlocks);
                deltaChoiBlocksString = cell(1,numBlocks);

                %% block diagonal set up which cant be placed in a separate function
                for indexBlock = 1:numBlocks
                    permMat = options.blockP{indexBlock};
                    newDims = options.newDims{indexBlock}.';
    
                    deltaChoiStrings = compose("deltaChoi%1d%1d(%d,%d)",indexBlock,(1:numel(newDims)).',newDims,newDims);
                    % CVX MUST generates variables in this work space. We have to do each
                    % separately.
                    for index = 1:numel(deltaChoiStrings)
                        variable(convertStringsToChars(deltaChoiStrings(index)),'hermitian')
                    end
                    % Get a single cell array to capture those variables.
                    deltaChoiBlocksString{indexBlock} = "{"+strjoin(compose("deltaChoi%d%d",indexBlock,1:numel(newDims)),",")+"}";
                    deltaChoiCurrentBlocks = eval(deltaChoiBlocksString{indexBlock});    
    
                    % there is a bug in cvx's blkdiag override that messes up some times.
                    % We can get around it by doing block diagonal manually.
                    deltaChoi = expression(compose("deltaChoi(%1d,%1d)",dimsAPrimeB(indexBlock),dimsAPrimeB(indexBlock)));
                    deltaChoi = directSumWorkArround(deltaChoi, deltaChoiCurrentBlocks);
    
                    % unscramble from block diagonal form.
                    deltaChoi = permMat'*deltaChoi*permMat;
                    deltaChoi = (deltaChoi+deltaChoi')/2;
                    deltaChoiBlocks{indexBlock} = deltaChoi;
                end

            else
                deltaChoiStrings = compose("deltaChoi%1d(%d,%d)",(1:numBlocks).',dimsAPrimeB.',dimsAPrimeB.');
                % CVX MUST generates variables in this work space. We have to do each
                % separately.
                for index = 1:numBlocks
                    variable(convertStringsToChars(deltaChoiStrings(index)),'hermitian');
                end
                % Get a single cell array to capture those variables.
                deltaChoiBlocksString = "{"+strjoin(compose("deltaChoi%d",1:numel(dimsAPrimeB)),",")+"}";
                deltaChoiBlocks = eval(deltaChoiBlocksString);
            end


            %% objective
            if feasibleMode
                minimize(deltaPhi);
            else
                temp = cellfun(@(x,y) real(trace(x*y)),gradChoiEveAPrimeB,deltaChoiBlocks,...
                    "UniformOutput",false);
                minimize(gradPhi*deltaPhi + sum([temp{:}]));
            end

            %% constraints

            % unique acceptance penalty
            phi+deltaPhi >= penaltyTermUniqueDecoy(rhoCTestMat,testCons,testProb,decoyProbs);
            phi+deltaPhi >= 0;

            % constraints from test states and channel
            rhoABTests = cell(size(testStatesAAPrime));

            %For squashing we also need states from generation rounds
            rhoABGens = cell(size(genStatesAAprime));

            for index = 1:numBlocks
                
                %Dimensions of current block
                dimA = dimsA(index);
                dimAPrime = dimsAPrime(index);
                dimAPrimeB = dimsAPrimeB(index);
                dimB = dimsB(index);
                choiTemp = deltaChoiBlocks{index} + choiEveAPrimeB{index};
                
                %test state of current block
                testStateAAPrime = testStatesAAPrime{index};
                
                %gen state of current block (only needed for squashing)
                genStateAAPrime = genStatesAAprime{index};

                % Eve's channels are CPTP
                choiTemp == hermitian_semidefinite(double(dimAPrimeB));
                PartialTrace(choiTemp, 2,[dimAPrime,dimB]) == eye(dimAPrime);

                % Calculate rhoAB conditioned on test for each block.
                % This is needed for the decoy constraints.
                rhoABTests{index} = PartialMap(testStateAAPrime,...
                    choiTemp,2,[dimA,dimAPrime]);
                
                %if squashing constraints are supplied calculate rho
                %conditioned on gen rounds
                if ~isempty(squashingConsGen)
                    rhoABGens{index} = PartialMap(genStateAAPrime,...
                    choiTemp,2,[dimA,dimAPrime]);
                end
            end

            % Decoy constraints
            decoyConstraints(cvx_problem, YieldMat, rhoCTestMat, rhoABTests,...
                probDistPhotonConMuTest, probSignalConTest, testCons, blockPhotonNum, linConTol)
            
            % Squashing constraints (apply only if they exist)
            if ~isempty(squashingConsTest) || ~isempty(squashingConsGen)
                squashingConstraints(cvx_problem, rhoABTests, rhoABGens, ...
                    squashingConsTest, squashingConsGen, numBlocks, linConTol);
            end

            cvx_end

            % reconstruct rho from blocks incase it didn't get filled in by CVX
            if options.blockDiagonal                
                deltaChoiBlocks = cell(1,numBlocks);
                for indexBlock = 1:numBlocks
                    %pick permutation matrix of current block
                    permMat = options.blockP{indexBlock};
                    %delta choi of current block
                    deltaChoi = zeros(dimsAPrimeB(indexBlock));

                    %evaluate cvx variables of current block (required
                    %because deltaChoiBlocks is only a derived variable)
                    deltaChoiCurrentBlocks = eval(deltaChoiBlocksString{indexBlock});

                    %stack together in block diagonal form and undo
                    %permutation
                    deltaChoi = directSumWorkArround(deltaChoi, deltaChoiCurrentBlocks);
                    deltaChoi = permMat'*deltaChoi*permMat;
                    deltaChoiBlocks{indexBlock} = (deltaChoi+deltaChoi')/2;
                end
            else
                %Reconstruct Blocks
                deltaChoiBlocks = eval(deltaChoiBlocksString);
            end            
            
            %Reformat into vector
            deltaChoiBlocks = cell2mat(cellfun(@(x) x(:),deltaChoiBlocks, "uniformoutput", false).');
            deltaVecFW = [deltaChoiBlocks;deltaPhi];

            statusMessage = string(cvx_status);

            switch statusMessage
                case "Solved"
                    exitFlag = SubProblemExitFlag.solved;
                case "Inaccurate/Solved"
                    exitFlag = SubProblemExitFlag.inaccurate;
                case "Failed"
                    exitFlag = SubProblemExitFlag.failed;
                otherwise
                    exitFlag = SubProblemExitFlag.failed;
            end
        end

        function [deltaVecFW,exitFlag,statusMessage] = subproblemUniqueDecoyIntImperfect(vecFW,vecFWGrad,...
                testStatesAAPrime,genStatesAAprime,dimsAPrime,testProb,testCons,squashingConsTest,squashingConsGen,....
                probDistPhotonConMuTestLower, probDistPhotonConMuTestUpper, probRemaining,...
                decoyProbs,probSignalConTest,blockPhotonNum,options,feasibleMode)
            arguments
                vecFW (:,1) double
                vecFWGrad (:,1) double
                testStatesAAPrime (1,:) cell {mustBeNonempty,mustBeCellOf(testStatesAAPrime,"double")}% cell of (:,:) double
                genStatesAAprime (1,:) cell {mustBeNonempty,mustBeCellOf(genStatesAAprime,"double")}
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
                testProb (1,1) double {mustBeInRange(testProb, 0, 1)}
                testCons (:,1) EqualityConstraintDecoy
                squashingConsTest (:,1) InequalityConstraint
                squashingConsGen (:,1) InequalityConstraint
                probDistPhotonConMuTestLower (:,:) double {mustBeNonnegative}
                probDistPhotonConMuTestUpper (:,:) double {mustBeNonnegative}
                probRemaining (1,:) double {mustBeNonnegative}
                decoyProbs (:,1) double {mustBeProbDist}
                probSignalConTest (:,1) double {mustBeProbDist}
                blockPhotonNum (:,1) uint64
                options (1,1) struct
                feasibleMode (1,1) logical = false;
            end

            %%%% !!!!!Testcons need to reflect block structure!!!!!
            
            %Linear constraint tolerance
            linConTol = options.linearConstraintTolerance;

            %Dimensions of the Choi matrices
            dimsAB = size(testCons(1).operator,1);
            dimsAAPrime = cellfun(@(x) size(x,1), testStatesAAPrime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;


            %read of various sizes
            photonCutoff = size(probDistPhotonConMuTestLower,1)-1; %Photon Cutoff
            numObs = numel(testCons); %Number of observation
            numTestInt = numel(decoyProbs); %Number of intensities in test rounds
            numBlocks = numel(dimsAPrimeB); %Number of Blocks

            % unpack all the values from vecFW and vecFWGrad
            [phi,choiEveAPrimeB]         = unpackVecFW(vecFW,    dimsAPrimeB);
            [gradPhi,gradChoiEveAPrimeB] = unpackVecFW(vecFWGrad,dimsAPrimeB);

            %% CVX
            cvx_begin
            cvx_solver(convertStringsToChars(options.cvxSolver));
            cvx_precision(convertStringsToChars(options.cvxPrecision));
            cvx_quiet(options.verboseLevel<2);


            %% variables
            deltaPhi = variable('deltaPhi(1,1)');
            YieldMat = variable(sprintf('YieldMat(%d,%d)',numObs,photonCutoff+1), 'nonnegative');
            rhoCTestMat = variable(sprintf('rhoCTestMat(%d,%d)',numObs,numTestInt), 'nonnegative');

            if options.blockDiagonal
                %define each block and its corresponding string as elements
                %of cell array
                deltaChoiBlocks = cell(1,numBlocks);
                deltaChoiBlocksString = cell(1,numBlocks);

                %% block diagonal set up which cant be placed in a separate function
                for indexBlock = 1:numBlocks
                    permMat = options.blockP{indexBlock};
                    newDims = options.newDims{indexBlock}.';
    
                    deltaChoiStrings = compose("deltaChoi%1d%1d(%d,%d)",indexBlock,(1:numel(newDims)).',newDims,newDims);
                    % CVX MUST generates variables in this work space. We have to do each
                    % separately.
                    for index = 1:numel(deltaChoiStrings)
                        variable(convertStringsToChars(deltaChoiStrings(index)),'hermitian')
                    end
                    % Get a single cell array to capture those variables.
                    deltaChoiBlocksString{indexBlock} = "{"+strjoin(compose("deltaChoi%d%d",indexBlock,1:numel(newDims)),",")+"}";
                    deltaChoiCurrentBlocks = eval(deltaChoiBlocksString{indexBlock});    
    
                    % there is a bug in cvx's blkdiag override that messes up some times.
                    % We can get around it by doing block diagonal manually.
                    deltaChoi = expression(compose("deltaChoi(%1d,%1d)",dimsAPrimeB(indexBlock),dimsAPrimeB(indexBlock)));
                    deltaChoi = directSumWorkArround(deltaChoi, deltaChoiCurrentBlocks);
    
                    % unscramble from block diagonal form.
                    deltaChoi = permMat'*deltaChoi*permMat;
                    deltaChoi = (deltaChoi+deltaChoi')/2;
                    deltaChoiBlocks{indexBlock} = deltaChoi;
                end

            else
                deltaChoiStrings = compose("deltaChoi%1d(%d,%d)",(1:numBlocks).',dimsAPrimeB.',dimsAPrimeB.');
                % CVX MUST generates variables in this work space. We have to do each
                % separately.
                for index = 1:numBlocks
                    variable(convertStringsToChars(deltaChoiStrings(index)),'hermitian');
                end
                % Get a single cell array to capture those variables.
                deltaChoiBlocksString = "{"+strjoin(compose("deltaChoi%d",1:numel(dimsAPrimeB)),",")+"}";
                deltaChoiBlocks = eval(deltaChoiBlocksString);
            end


            %% objective
            if feasibleMode
                minimize(deltaPhi);
            else
                temp = cellfun(@(x,y) real(trace(x*y)),gradChoiEveAPrimeB,deltaChoiBlocks,...
                    "UniformOutput",false);
                minimize(gradPhi*deltaPhi + sum([temp{:}]));
            end

            %% constraints

            % unique acceptance penalty
            phi+deltaPhi >= penaltyTermUniqueDecoy(rhoCTestMat,testCons,testProb,decoyProbs);
            phi+deltaPhi >= 0;

            % constraints from test states and channel
            rhoABTests = cell(size(testStatesAAPrime));

            %For squashing we also need states from generation rounds
            rhoABGens = cell(size(genStatesAAprime));

            for index = 1:numBlocks
                
                %Dimensions of current block
                dimA = dimsA(index);
                dimAPrime = dimsAPrime(index);
                dimAPrimeB = dimsAPrimeB(index);
                dimB = dimsB(index);
                choiTemp = deltaChoiBlocks{index} + choiEveAPrimeB{index};
                
                %test state of current block
                testStateAAPrime = testStatesAAPrime{index};
                
                %gen state of current block (only needed for squashing)
                genStateAAPrime = genStatesAAprime{index};

                % Eve's channels are CPTP
                choiTemp == hermitian_semidefinite(double(dimAPrimeB));
                PartialTrace(choiTemp, 2,[dimAPrime,dimB]) == eye(dimAPrime);

                % Calculate rhoAB conditioned on test for each block.
                % This is needed for the decoy constraints.
                rhoABTests{index} = PartialMap(testStateAAPrime,...
                    choiTemp,2,[dimA,dimAPrime]);
                
                %if squashing constraints are supplied calculate rho
                %conditioned on gen rounds
                if ~isempty(squashingConsGen)
                    rhoABGens{index} = PartialMap(genStateAAPrime,...
                    choiTemp,2,[dimA,dimAPrime]);
                end
            end

            % Decoy constraints
            decoyConstraintsIntImperfect(cvx_problem, YieldMat, rhoCTestMat, rhoABTests,...
                probDistPhotonConMuTestLower, probDistPhotonConMuTestUpper,probRemaining, probSignalConTest, testCons, blockPhotonNum, linConTol)
            
            % Squashing constraints (apply only if they exist)
            if ~isempty(squashingConsTest) || ~isempty(squashingConsGen)
                squashingConstraints(cvx_problem, rhoABTests, rhoABGens, ...
                    squashingConsTest, squashingConsGen, numBlocks, linConTol);
            end

            cvx_end

            % reconstruct rho from blocks incase it didn't get filled in by CVX
            if options.blockDiagonal                
                deltaChoiBlocks = cell(1,numBlocks);
                for indexBlock = 1:numBlocks
                    %pick permutation matrix of current block
                    permMat = options.blockP{indexBlock};
                    %delta choi of current block
                    deltaChoi = zeros(dimsAPrimeB(indexBlock));

                    %evaluate cvx variables of current block (required
                    %because deltaChoiBlocks is only a derived variable)
                    deltaChoiCurrentBlocks = eval(deltaChoiBlocksString{indexBlock});

                    %stack together in block diagonal form and undo
                    %permutation
                    deltaChoi = directSumWorkArround(deltaChoi, deltaChoiCurrentBlocks);
                    deltaChoi = permMat'*deltaChoi*permMat;
                    deltaChoiBlocks{indexBlock} = (deltaChoi+deltaChoi')/2;
                end
            else
                %Reconstruct Blocks
                deltaChoiBlocks = eval(deltaChoiBlocksString);
            end            
            
            %Reformat into vector
            deltaChoiBlocks = cell2mat(cellfun(@(x) x(:),deltaChoiBlocks, "uniformoutput", false).');
            deltaVecFW = [deltaChoiBlocks;deltaPhi];

            statusMessage = string(cvx_status);

            switch statusMessage
                case "Solved"
                    exitFlag = SubProblemExitFlag.solved;
                case "Inaccurate/Solved"
                    exitFlag = SubProblemExitFlag.inaccurate;
                case "Failed"
                    exitFlag = SubProblemExitFlag.failed;
                otherwise
                    exitFlag = SubProblemExitFlag.failed;
            end
        end
        
    end
end


function cost = privacyAmpSplit(phi,choiEveAPrimeB,genStateAAPrime,keyProj,krausOps,testProb,alpha,dimAPrime,perturbation,perturbationAff,probBlocks)
arguments
    % minimal checks for shape and type
    phi (1,1) double
    choiEveAPrimeB (1,:) cell {mustBeCellOf(choiEveAPrimeB,"double")}% cell of (:,:) double
    genStateAAPrime (1,:) cell {mustBeCellOf(genStateAAPrime,"double")}% cell of (:,:) double
    keyProj (1,:) cell {mustBeCellOf(keyProj,"cell")}% cell (:,1) cell of (:,:) double
    krausOps (1,:) cell {mustBeCellOf(krausOps,"cell")}% cell of (:,1) cell of (:,:) double
    testProb (1,1) double
    alpha (1,1) double
    dimAPrime (1,:) uint64
    perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
    perturbationAff (1,1) double {mustBeInRange(perturbationAff,0,1)}
    probBlocks (1,:) double {mustBeNonnegative,mustBeNonempty} = 1
end

rhoAB = cell(1,numel(probBlocks));
for index=1:numel(probBlocks)
    dimA = size(genStateAAPrime{index},1)/dimAPrime(index);
    rhoAB{index} = PartialMap(genStateAAPrime{index},choiEveAPrimeB{index},2,[dimA,dimAPrime(index)]);
end

%Concatenate rhoAB and keyProj and krausOps as repeating arguments,
%i.e. (a1,b1,a2,b2,...)
repeatingArgs = [rhoAB;keyProj;krausOps];

%version with alphaHat
% cost = phi/(alpha-1) + (1-testProb)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,probBlocks,repeatingArgs{:});

%version with alpha
cost = alpha/(alpha-1)*phi + (1-testProb)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,probBlocks,repeatingArgs{:});
end

function [gradByChoi,gradByPhi] = gradPrivacyAmpSplit(~,choiEvesAPrimeB,genStatesAAPrime,keyProjs,krausOpss,testProb,alpha,dimsAPrime,perturbation,perturbationAff,probBlocks)
    arguments
        % minimal checks for shape and type
        ~ %phi (1,1) double Not needed
        choiEvesAPrimeB (1,:) cell {mustBeCellOf(choiEvesAPrimeB,"double")}% cell of (:,:) double
        genStatesAAPrime (1,:) cell {mustBeCellOf(genStatesAAPrime,"double")}% cell of (:,:) double
        keyProjs (1,:) cell {mustBeCellOf(keyProjs,"cell")}% cell (:,1) cell of (:,:) double
        krausOpss (1,:) cell {mustBeCellOf(krausOpss,"cell")}% cell of (:,1) cell of (:,:) double
        testProb (1,1) double
        alpha (1,1) double
        dimsAPrime (1,:) uint64
        perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
        perturbationAff (1,1) double {mustBeInRange(perturbationAff,0,1)}
        probBlocks (1,:) double {mustBeNonnegative,mustBeNonempty} = 1
    end
    %gradient of phi term
    %version with alphaHat
    % gradByPhi = 1/(alpha-1);

    %version with alpha
    gradByPhi = alpha/(alpha-1);
    
    %gradient by Choi
    rhoAB = cell(1,numel(probBlocks));
    dimB = zeros(1,numel(probBlocks));
    for index=1:numel(probBlocks)
        dimA = size(genStatesAAPrime{index},1)/dimsAPrime(index);

        rhoAB{index} = PartialMap(genStatesAAPrime{index},choiEvesAPrimeB{index},2,[dimA,dimsAPrime(index)]);
        dimB(index) = size(rhoAB{index},1)/dimA;
    end
    
    %Concatenate rhoAB and keyProj and krausOps as repeating arguments,
    %i.e. (a1,b1,a2,b2,...)
    repeatingArgs = [rhoAB;keyProjs;krausOpss];
    
    gradByChoi = RenyiTildeDown.gradConEnt(alpha,perturbation,perturbationAff,probBlocks,repeatingArgs{:});
    for index=1:numel(probBlocks)
        dimA = size(genStatesAAPrime{index},1)/dimsAPrime(index);

        % The choi matrix must be ordered input system tensor output
        % system.
        genChoiAPrimeA = Swap(genStatesAAPrime{index},[1,2],[dimA,dimsAPrime(index)]);

        dualGenChoi = DualMap(genChoiAPrimeA,[dimsAPrime(index),dimA]);
        gradByChoi{index} = (1-testProb)*PartialMap(gradByChoi{index},dualGenChoi,1,[dimA,dimB(index)]);
    end
end

function penaltyVal = penaltyTermUniqueDecoy(rhoCTestMat,testCons,testProb,decoyProbs)
    arguments
        rhoCTestMat (:,:) % double or CVX
        testCons (:,1) EqualityConstraintDecoy
        testProb (1,1) double
        decoyProbs (:,1) double
    end
    
    % q = (gamma* p(mu_i|test) * p(a,b|test,mu)
    acceptPoint = reshape(testProb*cat(1,testCons.vector)*diag(decoyProbs),[],1);
    
    % rho_C = (gamma* p(mu_i|test) * rho_C|Mu,Test
    rhoCTestMatJoint = reshape(testProb*rhoCTestMat*diag(decoyProbs),[],1);
    
    penaltyVal = sum(rel_entr(acceptPoint,rhoCTestMatJoint))/log(2)...
        ...% this part should always be zero for unique acceptance
        + (1-testProb)*log2((1-testProb)/(1-testProb));
    end
    
function decoyConstraints(cvx_problem, YieldMat, rhoCTestMat, rhoABTests,...
        probDistPhotonConMuTest, probSignalConTest, testCons, blockPhotonNum, linConTol)
    arguments
        cvx_problem (1,1)
        YieldMat (:,:)
        rhoCTestMat (:,:)
        rhoABTests
        probDistPhotonConMuTest
        probSignalConTest (:,1) double
        testCons (:,1) EqualityConstraintDecoy
        blockPhotonNum (:,1) uint64
        linConTol (1,1) double {mustBeNonnegative}
    end
    
    
    photonCutoff = size(probDistPhotonConMuTest,1)-1;
    numSignalsAlice = numel(probSignalConTest);
    numObs = numel(testCons);
    
    %sum of yields is close (<= 1- ptot(mu)) to statistics generated from state
    pTotMu = sum(probDistPhotonConMuTest,1);
    probRemaining = ones(numObs,1)*(1-pTotMu);
    
    % equivalent to 0 <= rhoCTestMat - YieldMat*probDistPhotonConMuTest <= ones(numObs,1)*(1-pTotMu);
    YieldMat*probDistPhotonConMuTest <= rhoCTestMat                 + linConTol;
    YieldMat*probDistPhotonConMuTest >= rhoCTestMat - probRemaining - linConTol;
    
    %sum_b Y_n(a,b) = p(a|test,n)
    YieldTensor = reshape(YieldMat,numSignalsAlice,[],photonCutoff+1);
    % squeeze(sum(YieldTensor,2)) == probSignalConTest*ones(1,photonCutoff+1);
    abs(squeeze(sum(YieldTensor,2)) - probSignalConTest*ones(1,photonCutoff+1)) <= linConTol;
    
    %Y_n(a,b) <= p(a|test,n) Not needed, but sometimes it works better. No clue
    %why.
    YieldTensor <= repmat(probSignalConTest,1,size(YieldTensor,2),size(YieldTensor,3));
    
    
    %Trace constraints on e.g. single photons
    for index=1:numel(blockPhotonNum)
        for indexCons = 1:numObs
            YieldMat(indexCons,double(blockPhotonNum(index)+1)) == ...
                real(trace(testCons(indexCons).operator*rhoABTests{index}));
        end
    end
end


function decoyConstraintsIntImperfect(cvx_problem, YieldMat, rhoCTestMat, rhoABTests,...
        probDistPhotonConMuTestLower, probDistPhotonConMuTestUpper, probRemaining, probSignalConTest, testCons, blockPhotonNum, linConTol)
    arguments
        cvx_problem (1,1)
        YieldMat (:,:)
        rhoCTestMat (:,:)
        rhoABTests
        probDistPhotonConMuTestLower
        probDistPhotonConMuTestUpper
        probRemaining (1,:) double
        probSignalConTest (:,1) double        
        testCons (:,1) EqualityConstraintDecoy
        blockPhotonNum (:,1) uint64
        linConTol (1,1) double {mustBeNonnegative}
    end
    
    
    photonCutoff = size(probDistPhotonConMuTestLower,1)-1;
    numSignalsAlice = numel(probSignalConTest);
    numObs = numel(testCons);
    
    %sum of yields is close (<= 1- ptot(mu)) to statistics generated from state
    probRemainingReshaped = ones(numObs,1)*(1-probRemaining);
    
    % equivalent to 0 <= rhoCTestMat - YieldMat*probDistPhotonConMuTest <= ones(numObs,1)*(1-pTotMu);
    YieldMat*probDistPhotonConMuTestLower <= rhoCTestMat                         + linConTol;
    YieldMat*probDistPhotonConMuTestUpper >= rhoCTestMat - probRemainingReshaped - linConTol;
    
    %sum_b Y_n(a,b) = p(a|test,n)
    YieldTensor = reshape(YieldMat,numSignalsAlice,[],photonCutoff+1);
    % squeeze(sum(YieldTensor,2)) == probSignalConTest*ones(1,photonCutoff+1);
    abs(squeeze(sum(YieldTensor,2)) - probSignalConTest*ones(1,photonCutoff+1)) <= linConTol;
    
    %Y_n(a,b) <= p(a|test,n) Not needed, but sometimes it works better. No clue
    %why.
    YieldTensor <= repmat(probSignalConTest,1,size(YieldTensor,2),size(YieldTensor,3));
    
    
    %Trace constraints on e.g. single photons
    for index=1:numel(blockPhotonNum)
        for indexCons = 1:numObs
            YieldMat(indexCons,double(blockPhotonNum(index)+1)) == ...
                real(trace(testCons(indexCons).operator*rhoABTests{index}));
        end
    end
end


function squashingConstraints(cvx_problem, rhoABTests, rhoABGens, ...
        squashingConsTest, squashingConsGen, numBlocks, linConTol)
    arguments
        cvx_problem (1,1)            
        rhoABTests
        rhoABGens
        squashingConsTest (:,1) InequalityConstraint
        squashingConsGen (:,1) InequalityConstraint
        numBlocks (1,1) uint64
        linConTol (1,1) double {mustBeNonnegative}
    end
    % Squashing constraints

    % Create squashing constraints for test rounds
    for indexBlock = 1:numBlocks
        %iterate over constraints of Gen/Test block
        for indexCons = 1:numel(squashingConsTest)
            if ~isinf(squashingConsTest(indexCons).upperBound) %upper bound
                real(trace(squashingConsTest(indexCons).operator'*(rhoABTests{indexBlock}))) <= ...
                    squashingConsTest(indexCons).upperBound + linConTol;
            end
            if ~isinf(squashingConsTest(indexCons).lowerBound) %lower bound
                real(trace(squashingConsTest(indexCons).operator'*(rhoABTests{indexBlock}))) >= ... 
                    squashingConsTest(indexCons).lowerBound - linConTol;
            end
        end
    end

    % Create squashing constraints for gen rounds
    for indexBlock = 1:numBlocks
        %iterate over constraints of Gen/Test block
        for indexCons = 1:numel(squashingConsGen)
            if ~isinf(squashingConsGen(indexCons).upperBound) %upper bound
                real(trace(squashingConsGen(indexCons).operator'*(rhoABGens{indexBlock}))) <= ...
                    squashingConsGen(indexCons).upperBound + linConTol;
            end
            if ~isinf(squashingConsGen(indexCons).lowerBound) %lower bound
                real(trace(squashingConsGen(indexCons).operator'*(rhoABGens{indexBlock}))) >= ... 
                    squashingConsGen(indexCons).lowerBound - linConTol;
            end
        end
    end
end


function [phi,choiMatrices] = unpackVecFW(vecFW,dimsAPrimeB)

phi = vecFW(end);

numBlocks = numel(dimsAPrimeB);
choiMatrices = cell(1,numBlocks);

totalIndex = 1;
for index = 1:numBlocks
    %Get dimensions of A'B of current block
    currentDim = dimsAPrimeB(index);
    lastIndex = totalIndex+currentDim^2-1;

    %Pick current choi matrix and reshape it.
    choiMatrices{index} = reshape(vecFW(totalIndex:lastIndex), ...
        [currentDim,currentDim]);

    totalIndex = lastIndex+1;
end
end