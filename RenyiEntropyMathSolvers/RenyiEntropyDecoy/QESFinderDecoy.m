classdef QESFinderDecoy
    methods (Static)
        function [gradQES,tol,exitFlag] = getGradQES(vecFW,vecFWGrad,...
                testStatesAAPrime,genStatesAAprime,dimsAPrime,testProb,testCons, ...
                squashingConsTest,squashingConsGen,probDistPhotonConMuTest, ...
                decoyProbs,probSignalConTest,blockPhotonNum,options)
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
            end

            %Linear constraint tolerance
            linConTol = 1/2*options.linearConstraintTolerance;

            %Dimensions of the Choi matrices
            dimsAB = size(testCons(1).operator,1);
            dimsAAPrime = cellfun(@(x) size(x,1), testStatesAAPrime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;

            %Photon Cutoff
            photonCutoff = size(probDistPhotonConMuTest,1)-1;

            %Number of observation
            numObs = numel(testCons);

            %Number of intensities in test rounds
            numTestInt = numel(decoyProbs);

            %Number of Blocks
            numBlocks = numel(dimsAPrimeB);

            %Number of Signals sent by Alice in test rounds
            numSignalsAlice = numel(probSignalConTest);


            %Define accept point
            % q = (gamma* p(mu_i|test) * p(a,b|test,mu), 1-gamma)^T
            acceptPoint = [reshape(testProb*cat(1,testCons.vector)*diag(decoyProbs),[],1); 1-testProb];
            
            % unpack all the values from vecFW and vecFWGrad
            [phi,choiEveAPrimeB]         = unpackVecFW(vecFW,    dimsAPrimeB);
            [gradPhi,gradChoiEveAPrimeB] = unpackVecFW(vecFWGrad,dimsAPrimeB);

            %% CVX
            cvx_begin
            cvx_solver(convertStringsToChars(options.cvxSolver));
            cvx_precision(convertStringsToChars(options.cvxPrecision));
            cvx_quiet(options.verboseLevel<2);
            
            lambda = variable('lambda(numObs*numTestInt+1,1)');
            dualLambda = dual('variable','dualLambda'); %declares dual variable
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
            temp = cellfun(@(x,y) real(trace(x*y)),gradChoiEveAPrimeB,deltaChoiBlocks,"UniformOutput",false);
            minimize(gradPhi*deltaPhi + sum([temp{:}]));

            %% constraints

            % unique acceptance penalty
            phi+deltaPhi >= penaltyTermVariableLengthDecoy(cvx_problem,rhoCTestMat,testProb,decoyProbs,lambda);
            phi+deltaPhi >= 0;           
                        
            % constraints from test states and channel
            rhoABTests = cell(size(testStatesAAPrime));

            %For squashing we also need states from generation rounds
            rhoABGens = cell(size(genStatesAAprime));

            for index = 1:numBlocks
                %Find dimensions of current block
                currentdimA = dimsA(index);
                currentdimAPrime = dimsAPrime(index);
                currentdimAPrimeB = dimsAPrimeB(index);
                currentdimB = dimsB(index);
                choiTemp = deltaChoiBlocks{index} + choiEveAPrimeB{index};
                
                %test state of current block
                currentTestStateAAPrime = testStatesAAPrime{index};

                %gen state of current block (only needed for squashing)
                currentGenStateAAPrime = genStatesAAprime{index};

                % Eve's channels are CPTP
                choiTemp == hermitian_semidefinite(double(currentdimAPrimeB));
                PartialTrace(choiTemp, 2,[currentdimAPrime,currentdimB]) == eye(currentdimAPrime);

                % Calculate rhoAB conditioned on test for each block.
                % This is needed for the decoy constraints.
                rhoABTests{index} = PartialMap(currentTestStateAAPrime,...
                    choiTemp,2,[currentdimA,currentdimAPrime]);

                %if squashing constraints are supplied calculate rho
                %conditioned on gen rounds
                if ~isempty(squashingConsGen)
                    rhoABGens{index} = PartialMap(currentGenStateAAPrime,...
                    choiTemp,2,[currentdimA,currentdimAPrime]);
                end
            end
            
            % Decoy constraints
            decoyConstraints(cvx_problem, YieldMat, rhoCTestMat, rhoABTests,...
                probDistPhotonConMuTest, probSignalConTest, testCons, blockPhotonNum, linConTol);

            % Squashing constraints (apply only if they exist)
            if ~isempty(squashingConsTest) || ~isempty(squashingConsGen)
                squashingConstraints(cvx_problem, rhoABTests, rhoABGens, ...
                    squashingConsTest, squashingConsGen, numBlocks, linConTol);
            end

            %% Get gradQES from dual variables
            acceptPoint - lambda == 0 : dualLambda;
            
            cvx_end
            
            %store status message
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
            
            %Gradient of QES given by dual variables to q - lambda == 0
            %(cvx assings the sign the other way around....)

            gradQES = -full(dualLambda);

            %store tolerance of solution
            tol = cvx_slvtol;
        end

        function [gradQES,tol,exitFlag] = getGradQESIntImperfect(vecFW,vecFWGrad,...
                testStatesAAPrime,genStatesAAprime,dimsAPrime,testProb,testCons, ...
                squashingConsTest,squashingConsGen, ... 
                probDistPhotonConMuTestLower, probDistPhotonConMuTestUpper, probRemaining, ...
                decoyProbs,probSignalConTest,blockPhotonNum,options)
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
            end

            %Linear constraint tolerance
            linConTol = 1/2*options.linearConstraintTolerance;

            %Dimensions of the Choi matrices
            dimsAB = size(testCons(1).operator,1);
            dimsAAPrime = cellfun(@(x) size(x,1), testStatesAAPrime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;

            %Photon Cutoff
            photonCutoff = size(probDistPhotonConMuTestLower,1)-1;

            %Number of observation
            numObs = numel(testCons);

            %Number of intensities in test rounds
            numTestInt = numel(decoyProbs);

            %Number of Blocks
            numBlocks = numel(dimsAPrimeB);

            %Number of Signals sent by Alice in test rounds
            numSignalsAlice = numel(probSignalConTest);


            %Define accept point
            % q = (gamma* p(mu_i|test) * p(a,b|test,mu), 1-gamma)^T
            acceptPoint = [reshape(testProb*cat(1,testCons.vector)*diag(decoyProbs),[],1); 1-testProb];
            
            % unpack all the values from vecFW and vecFWGrad
            [phi,choiEveAPrimeB]         = unpackVecFW(vecFW,    dimsAPrimeB);
            [gradPhi,gradChoiEveAPrimeB] = unpackVecFW(vecFWGrad,dimsAPrimeB);

            %% CVX
            cvx_begin
            cvx_solver(convertStringsToChars(options.cvxSolver));
            cvx_precision(convertStringsToChars(options.cvxPrecision));
            cvx_quiet(options.verboseLevel<2);
            
            lambda = variable('lambda(numObs*numTestInt+1,1)');
            dualLambda = dual('variable','dualLambda'); %declares dual variable
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
            temp = cellfun(@(x,y) real(trace(x*y)),gradChoiEveAPrimeB,deltaChoiBlocks,"UniformOutput",false);
            minimize(gradPhi*deltaPhi + sum([temp{:}]));

            %% constraints

            % unique acceptance penalty
            phi+deltaPhi >= penaltyTermVariableLengthDecoy(cvx_problem,rhoCTestMat,testProb,decoyProbs,lambda);
            phi+deltaPhi >= 0;           
                        
            % constraints from test states and channel
            rhoABTests = cell(size(testStatesAAPrime));

            %For squashing we also need states from generation rounds
            rhoABGens = cell(size(genStatesAAprime));

            for index = 1:numBlocks
                %Find dimensions of current block
                currentdimA = dimsA(index);
                currentdimAPrime = dimsAPrime(index);
                currentdimAPrimeB = dimsAPrimeB(index);
                currentdimB = dimsB(index);
                choiTemp = deltaChoiBlocks{index} + choiEveAPrimeB{index};
                
                %test state of current block
                currentTestStateAAPrime = testStatesAAPrime{index};

                %gen state of current block (only needed for squashing)
                currentGenStateAAPrime = genStatesAAprime{index};

                % Eve's channels are CPTP
                choiTemp == hermitian_semidefinite(double(currentdimAPrimeB));
                PartialTrace(choiTemp, 2,[currentdimAPrime,currentdimB]) == eye(currentdimAPrime);

                % Calculate rhoAB conditioned on test for each block.
                % This is needed for the decoy constraints.
                rhoABTests{index} = PartialMap(currentTestStateAAPrime,...
                    choiTemp,2,[currentdimA,currentdimAPrime]);

                %if squashing constraints are supplied calculate rho
                %conditioned on gen rounds
                if ~isempty(squashingConsGen)
                    rhoABGens{index} = PartialMap(currentGenStateAAPrime,...
                    choiTemp,2,[currentdimA,currentdimAPrime]);
                end
            end
            
            % Decoy constraints
            decoyConstraintsIntImperfect(cvx_problem, YieldMat, rhoCTestMat, rhoABTests,...
                probDistPhotonConMuTestLower, probDistPhotonConMuTestUpper,probRemaining, probSignalConTest, testCons, blockPhotonNum, linConTol)
            
            % rhoCTestMat is a probability distribution in each column
            abs(sum(rhoCTestMat,1) - ones(1,numTestInt)) <= linConTol

            % Squashing constraints (apply only if they exist)
            if ~isempty(squashingConsTest) || ~isempty(squashingConsGen)
                squashingConstraints(cvx_problem, rhoABTests, rhoABGens, ...
                    squashingConsTest, squashingConsGen, numBlocks, linConTol);
            end

            %% Get gradQES from dual variables
            acceptPoint - lambda == 0 : dualLambda;
            
            cvx_end
            
            %store status message
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
            
            %Gradient of QES given by dual variables to q - lambda == 0
            %(cvx assings the sign the other way around....)

            gradQES = -full(dualLambda);

            %store tolerance of solution
            tol = cvx_slvtol;
        end
    
        function val = funcFW(vecFWQes,genStateAAPrime,keyProj,krausOps,testCons, ...
                testProb,gradQES,alpha,dimsAPrime,decoyProbs,probBlocks,perturbation,perturbationAff)
            arguments
                vecFWQes (:,1) double
                genStateAAPrime (1,:) cell {mustBeCellOf(genStateAAPrime,"double")}
                keyProj (1,:) cell {mustBeCellOf(keyProj,"cell")}
                krausOps (1,:) cell {mustBeCellOf(krausOps,"cell")}
                testCons (:,1) EqualityConstraintDecoy
                testProb (1,1) double
                gradQES (:,1) double
                alpha (1,1) double
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
                decoyProbs (:,1) double {mustBeProbDist}
                probBlocks (:,1) double {mustBePositive}
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
                perturbationAff(1,1) double {mustBeInRange(perturbationAff,0,1)}                
            end

            %Dimensions of the Choi matrices
            dimsAB = size(testCons(1).operator,1);
            dimsAAPrime = cellfun(@(x) size(x,1), genStateAAPrime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;

            %Number of observation
            numObs = numel(testCons);

            %Number of intensities in test rounds
            numTestInt = numel(decoyProbs);

            %Number of Blocks
            numBlocks = numel(dimsAPrimeB);

            %unpack choi states and outputs in test rounds
            [choiEveAPrimeB,rhoCTestResults] = unpackVecFWQes(vecFWQes,dimsAPrimeB,numObs,numTestInt);
            
            %Find corresponding states rhoAB
            rhoAB = cell(1,numel(probBlocks));
            for index=1:numel(probBlocks)
                rhoAB{index} = PartialMap(genStateAAPrime{index},choiEveAPrimeB{index},...
                    2,[dimsA(index),dimsAPrime(index)]);
            end
            
            %Concatenate rhoAB and keyProj and krausOps as repeating arguments,
            %i.e. (a1,b1,a2,b2,...)
            repeatingArgs = [rhoAB;keyProj;krausOps];                        

            %rhoCtest*p(mu|test)            
            pMuRhoCTest = reshape(reshape(rhoCTestResults,[],numTestInt)*diag(decoyProbs),[],1).';
            sumTraceConTest = pMuRhoCTest*2.^((alpha-1)*gradQES(1:end-1));
            
            %Value of QES entropy
            val = -1/(alpha-1)*log2((1-testProb)*2^(-(alpha-1)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,probBlocks,repeatingArgs{:}))*2^((alpha-1)*gradQES(end)) ...
                + testProb*sumTraceConTest);
        end

        function val = funcFW2(vecFWQes,genStateAAPrime,keyProj,krausOps,testCons, ...
                testProb,gradQES,alpha,dimsAPrime,decoyProbs,probBlocks,perturbation,perturbationAff)
            arguments
                vecFWQes (:,1) double
                genStateAAPrime (1,:) cell {mustBeCellOf(genStateAAPrime,"double")}
                keyProj (1,:) cell {mustBeCellOf(keyProj,"cell")}
                krausOps (1,:) cell {mustBeCellOf(krausOps,"cell")}
                testCons (:,1) EqualityConstraintDecoy
                testProb (1,1) double
                gradQES (:,1) double
                alpha (1,1) double
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
                decoyProbs (:,1) double {mustBeProbDist}
                probBlocks (:,1) double {mustBePositive}
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
                perturbationAff(1,1) double {mustBeInRange(perturbationAff,0,1)}                
            end

            %Dimensions of the Choi matrices
            dimsAB = size(testCons(1).operator,1);
            dimsAAPrime = cellfun(@(x) size(x,1), genStateAAPrime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;

            %Number of observation
            numObs = numel(testCons);

            %Number of intensities in test rounds
            numTestInt = numel(decoyProbs);

            %Number of Blocks
            numBlocks = numel(dimsAPrimeB);

            %unpack choi states and outputs in test rounds
            [choiEveAPrimeB,rhoCTestResults] = unpackVecFWQes(vecFWQes,dimsAPrimeB,numObs,numTestInt);
            
            %Find corresponding states rhoAB
            rhoAB = cell(1,numel(probBlocks));
            for index=1:numel(probBlocks)
                rhoAB{index} = PartialMap(genStateAAPrime{index},choiEveAPrimeB{index},...
                    2,[dimsA(index),dimsAPrime(index)]);
            end
            
            %Concatenate rhoAB and keyProj and krausOps as repeating arguments,
            %i.e. (a1,b1,a2,b2,...)
            repeatingArgs = [rhoAB;keyProj;krausOps];                        

            %rhoCtest*p(mu|test)            
            pMuRhoCTest = reshape(reshape(rhoCTestResults,[],numTestInt)*diag(decoyProbs),[],1).';
            sumTraceConTest = pMuRhoCTest*2.^((alpha-1)*gradQES(1:end-1));
            
            %Value of QES entropy
            val = -((1-testProb)*2^(-(alpha-1)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,probBlocks,repeatingArgs{:}))*2^((alpha-1)*gradQES(end)) ...
                + testProb*sumTraceConTest);
        end

        function grad = gradFuncFW(vecFWQes,genStatesAAPrime,keyProj,krausOps,testCons, ...
                testProb,gradQES,alpha,dimsAPrime,decoyProbs,probBlocks,perturbation,perturbationAff)
            arguments
                vecFWQes (:,1) double
                genStatesAAPrime (1,:) cell {mustBeCellOf(genStatesAAPrime,"double")}
                keyProj (1,:) cell {mustBeCellOf(keyProj,"cell")}
                krausOps (1,:) cell {mustBeCellOf(krausOps,"cell")}
                testCons (:,1) EqualityConstraintDecoy
                testProb (1,1) double
                gradQES (:,1) double
                alpha (1,1) double
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
                decoyProbs (:,1) double {mustBeProbDist}
                probBlocks (:,1) double {mustBePositive}
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
                perturbationAff(1,1) double {mustBeInRange(perturbationAff,0,1)} 
            end

            %Dimensions of the Choi matrices
            dimsAB = size(testCons(1).operator,1);
            dimsAAPrime = cellfun(@(x) size(x,1), genStatesAAPrime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;

            %Number of observation
            numObs = numel(testCons);

            %Number of intensities in test rounds
            numTestInt = numel(decoyProbs);

            %unpack choi states and outputs in test rounds
            [choiEveAPrimeB,rhoCTestResults] = unpackVecFWQes(vecFWQes,dimsAPrimeB,numObs,numTestInt);
            
            %Find corresponding states rhoAB
            rhoAB = cell(1,numel(probBlocks));
            for index=1:numel(probBlocks)
                currentdimA = size(genStatesAAPrime{index},1)/dimsAPrime(index);
                rhoAB{index} = PartialMap(genStatesAAPrime{index},choiEveAPrimeB{index},2,[currentdimA,dimsAPrime(index)]);
            end
            %Concatenate rhoAB and keyProj and krausOps as repeating arguments,
            %i.e. (a1,b1,a2,b2,...)
            repeatingArgs = [rhoAB;keyProj;krausOps];            

            %rhoCtest*p(mu|test)            
            pMuRhoCTest = reshape(reshape(rhoCTestResults,[],numTestInt)*diag(decoyProbs),[],1).';
            sumTraceConTest = pMuRhoCTest*2.^((alpha-1)*gradQES(1:end-1));


            %Calculate inner term of log
            innerTerm = (1-testProb)*2^(-(alpha-1)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,probBlocks,repeatingArgs{:}))*2^((alpha-1)*gradQES(end)) ...
                + testProb*sumTraceConTest;
            
            %Gradient by RhoCTest
            % Here we want the gradient and the function is morally
            % speaking, f(x) = const*x. So we need to construct the constant,
            % p(mu|test), in the same shape.
            pmuVec = reshape(ones(numObs,numTestInt)*diag(decoyProbs),[],1);
            gradByRhoCTest = -testProb/(alpha-1)/log(2)*1/innerTerm*2.^((alpha-1)*gradQES(1:end-1)).*pmuVec;

            %Gradient by Choi
            gradByChoi = RenyiTildeDown.gradConEnt(alpha,perturbation,perturbationAff,probBlocks,repeatingArgs{:});
            for index=1:numel(probBlocks)
                currentdimA = size(genStatesAAPrime{index},1)/dimsAPrime(index);
                
                % The choi matrix must be ordered input system tensor output
                % system.
                genChoiAPrimeA = Swap(genStatesAAPrime{index},[1,2],[currentdimA,dimsAPrime(index)]);        
                tempDualMap = DualMap(genChoiAPrimeA,[dimsAPrime(index),currentdimA]);

                gradConEntByChoi = PartialMap(gradByChoi{index},tempDualMap,1,[currentdimA,dimsB(index)]);
                
                gradByChoi{index} = (1-testProb)/innerTerm*2^(-(alpha-1)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,probBlocks,repeatingArgs{:}))...
                *2^((alpha-1)*gradQES(end))*gradConEntByChoi;
            end
            %Reshape each gradient in the cell gradbyChoi to be linearized
            gradByChoi = cell2mat(cellfun(@(x) x(:),gradByChoi, "uniformoutput", false));

            %Final gradient
            grad = [gradByChoi(:);gradByRhoCTest];
        end

        function grad = gradFuncFW2(vecFWQes,genStatesAAPrime,keyProj,krausOps,testCons, ...
                testProb,gradQES,alpha,dimsAPrime,decoyProbs,probBlocks,perturbation,perturbationAff)
            arguments
                vecFWQes (:,1) double
                genStatesAAPrime (1,:) cell {mustBeCellOf(genStatesAAPrime,"double")}
                keyProj (1,:) cell {mustBeCellOf(keyProj,"cell")}
                krausOps (1,:) cell {mustBeCellOf(krausOps,"cell")}
                testCons (:,1) EqualityConstraintDecoy
                testProb (1,1) double
                gradQES (:,1) double
                alpha (1,1) double
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
                decoyProbs (:,1) double {mustBeProbDist}
                probBlocks (:,1) double {mustBePositive}
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
                perturbationAff(1,1) double {mustBeInRange(perturbationAff,0,1)} 
            end

            %Dimensions of the Choi matrices
            dimsAB = size(testCons(1).operator,1);
            dimsAAPrime = cellfun(@(x) size(x,1), genStatesAAPrime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;

            %Number of observation
            numObs = numel(testCons);

            %Number of intensities in test rounds
            numTestInt = numel(decoyProbs);

            %unpack choi states and outputs in test rounds
            [choiEveAPrimeB,rhoCTestResults] = unpackVecFWQes(vecFWQes,dimsAPrimeB,numObs,numTestInt);
            
            %Find corresponding states rhoAB
            rhoAB = cell(1,numel(probBlocks));
            for index=1:numel(probBlocks)
                currentdimA = size(genStatesAAPrime{index},1)/dimsAPrime(index);
                rhoAB{index} = PartialMap(genStatesAAPrime{index},choiEveAPrimeB{index},2,[currentdimA,dimsAPrime(index)]);
            end
            %Concatenate rhoAB and keyProj and krausOps as repeating arguments,
            %i.e. (a1,b1,a2,b2,...)
            repeatingArgs = [rhoAB;keyProj;krausOps];            

            %Gradient by RhoCTest
            % Here we want the gradient and the function is morally
            % speaking, f(x) = const*x. so We need toff the constant,
            % p(mu|test), in the same shape.
            pmuVec = reshape(ones(numObs,numTestInt)*diag(decoyProbs),[],1);
            gradByRhoCTest = -testProb*2.^((alpha-1)*gradQES(1:end-1)).*pmuVec;

            %Gradient by Choi
            gradByChoi = RenyiTildeDown.gradConEnt(alpha,perturbation,perturbationAff,probBlocks,repeatingArgs{:});
            for index=1:numel(probBlocks)
                currentdimA = size(genStatesAAPrime{index},1)/dimsAPrime(index);
                
                % The choi matrix must be ordered input system tensor output
                % system.
                genChoiAPrimeA = Swap(genStatesAAPrime{index},[1,2],[currentdimA,dimsAPrime(index)]);        
                tempDualMap = DualMap(genChoiAPrimeA,[dimsAPrime(index),currentdimA]);

                gradConEntByChoi = PartialMap(gradByChoi{index},tempDualMap,1,[currentdimA,dimsB(index)]);
                
                gradByChoi{index} = (alpha-1)*log(2)*(1-testProb)*2^((alpha-1)*gradQES(end))*...
                                    2^(-(alpha-1)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,probBlocks,repeatingArgs{:}))*...
                                    gradConEntByChoi;
            end
            %Reshape each gradient in the cell gradbyChoi to be linearized
            gradByChoi = cell2mat(cellfun(@(x) x(:),gradByChoi, "uniformoutput", false));

            %Final gradient
            grad = [gradByChoi(:);gradByRhoCTest];
        end

        function [deltaVecFW,exitFlag,statusMessage,tol] = subproblemQes(vecFWQes,vecFWQesGrad,...
                testStatesAAPrime,genStatesAAprime,dimsAPrime,testCons,squashingConsTest,squashingConsGen, ...
                probDistPhotonConMuTest,decoyProbs,probSignalConTest,blockPhotonNum,options,feasibleMode)
            arguments
                vecFWQes (:,1) double
                vecFWQesGrad (:,1) double
                testStatesAAPrime (1,:) cell {mustBeNonempty,mustBeCellOf(testStatesAAPrime,"double")}% cell of (:,:) double 
                genStatesAAprime (1,:) cell {mustBeNonempty,mustBeCellOf(genStatesAAprime,"double")}
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
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

            %Linear constraint tolerance
            linConTol = options.linearConstraintTolerance;

            %Dimensions of the Choi matrices
            dimsAB = size(testCons(1).operator,1);
            dimsAAPrime = cellfun(@(x) size(x,1), testStatesAAPrime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;

            %Photon Cutoff
            photonCutoff = size(probDistPhotonConMuTest,1)-1;

            %Number of observation
            numObs = numel(testCons);

            %Number of intensities in test rounds
            numTestInt = numel(decoyProbs);

            %Number of Blocks
            numBlocks = numel(dimsAPrimeB);

            %Number of Signals sent by Alice in test rounds
            numSignalsAlice = numel(probSignalConTest);

            %unpack choi states and outputs in test rounds and 
            % gradient of choi states and gradient outputs in test rounds
            [choiEveAPrimeB,rhoCTestResults]    = unpackVecFWQes(vecFWQes,      dimsAPrimeB,numObs,numTestInt);
            [gradChoiEveAPrimeB,gradRhoC]       = unpackVecFWQes(vecFWQesGrad,  dimsAPrimeB,numObs,numTestInt);

            %RhoCTest reshaped into correct form for optimization
            %size(numobs,numint)
            rhoCTestResults = reshape(rhoCTestResults,[],numTestInt);

            %% CVX
            cvx_begin
            cvx_solver(convertStringsToChars(options.cvxSolver));
            cvx_precision(convertStringsToChars(options.cvxPrecision));
            cvx_quiet(options.verboseLevel<2);
            
            %% Define variables
            deltaRhoC = variable(sprintf('deltaRhoC(%d,%d)',numObs,numTestInt)); 
            YieldMat = variable(sprintf('YieldMat(%d,%d)',numObs,photonCutoff+1), 'nonnegative');
            
            %Sum of current rhoCTestResults and deltaRhoC is more
            %convenient still needs to be nonnegative
            rhoCTestMat = deltaRhoC + rhoCTestResults;            
            rhoCTestMat == nonnegative([numObs,numTestInt]);            

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
                minimize(1);
            else
                tempChoi = cellfun(@(x,y) real(trace(x*y)),gradChoiEveAPrimeB,deltaChoiBlocks,"UniformOutput",false);           
                minimize(gradRhoC.'*reshape(deltaRhoC,[],1) + sum([tempChoi{:}]));
            end
            
            %% constraints
            
            % constraints from test states and channel
            rhoABTests = cell(size(testStatesAAPrime));

            %For squashing we also need states from generation rounds
            rhoABGens = cell(size(genStatesAAprime));

            for index = 1:numBlocks
                %Find dimensions of current block
                currentdimA = dimsA(index);
                currentdimAPrime = dimsAPrime(index);
                currentdimAPrimeB = dimsAPrimeB(index);
                currentdimB = dimsB(index);
                choiTemp = deltaChoiBlocks{index} + choiEveAPrimeB{index};

                %test state of current block
                currentTestStateAAPrime = testStatesAAPrime{index};

                %gen state of current block (only needed for squashing)
                currentGenStateAAPrime = genStatesAAprime{index};

                % Eve's channels are CPTP
                choiTemp == hermitian_semidefinite(double(currentdimAPrimeB));
                PartialTrace(choiTemp, 2,[currentdimAPrime,currentdimB]) == eye(currentdimAPrime);

                % Calculate rhoAB conditioned on test for each block.
                % This is needed for the decoy constraints.
                rhoABTests{index} = PartialMap(currentTestStateAAPrime,...
                    choiTemp,2,[currentdimA,currentdimAPrime]);

                %if squashing constraints are supplied calculate rho
                %conditioned on gen rounds
                if ~isempty(squashingConsGen)
                    rhoABGens{index} = PartialMap(currentGenStateAAPrime,...
                    choiTemp,2,[currentdimA,currentdimAPrime]);
                end
            end

            %Decoy constraints
            decoyConstraints(cvx_problem, YieldMat, rhoCTestMat, rhoABTests,...
                probDistPhotonConMuTest, probSignalConTest, testCons, blockPhotonNum, linConTol);

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
            
            %Reshape deltaChoiBlocks into vectors
            deltaChoiBlocks = cell2mat(cellfun(@(x) x(:),deltaChoiBlocks, "uniformoutput", false).');

            %Output overall delta vector
            deltaVecFW = [deltaChoiBlocks;reshape(deltaRhoC,[],1)];

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

            %store tolerance of solution
            tol = cvx_slvtol;
        end

        function [deltaVecFW,exitFlag,statusMessage,tol] = subproblemQesIntImperfect(vecFWQes,vecFWQesGrad,...
                testStatesAAPrime,genStatesAAprime,dimsAPrime,testCons,squashingConsTest,squashingConsGen, ...
                probDistPhotonConMuTestLower,probDistPhotonConMuTestUpper,probRemaining,decoyProbs,probSignalConTest,blockPhotonNum,options,feasibleMode)
            arguments
                vecFWQes (:,1) double
                vecFWQesGrad (:,1) double
                testStatesAAPrime (1,:) cell {mustBeNonempty,mustBeCellOf(testStatesAAPrime,"double")}% cell of (:,:) double 
                genStatesAAprime (1,:) cell {mustBeNonempty,mustBeCellOf(genStatesAAprime,"double")}
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
                testCons (:,1) EqualityConstraintDecoy
                squashingConsTest (:,1) InequalityConstraint
                squashingConsGen (:,1) InequalityConstraint
                probDistPhotonConMuTestLower (:,:) double {mustBePositive}
                probDistPhotonConMuTestUpper (:,:) double {mustBePositive}
                probRemaining (1,:) double {mustBePositive}
                decoyProbs (:,1) double {mustBeProbDist}
                probSignalConTest (:,1) double {mustBeProbDist}                
                blockPhotonNum (:,1) uint64     
                options (1,1) struct
                feasibleMode (1,1) logical = false;
            end

            %Linear constraint tolerance
            linConTol = options.linearConstraintTolerance;

            %Dimensions of the Choi matrices
            dimsAB = size(testCons(1).operator,1);
            dimsAAPrime = cellfun(@(x) size(x,1), testStatesAAPrime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;

            %Photon Cutoff
            photonCutoff = size(probDistPhotonConMuTestLower,1)-1;

            %Number of observation
            numObs = numel(testCons);

            %Number of intensities in test rounds
            numTestInt = numel(decoyProbs);

            %Number of Blocks
            numBlocks = numel(dimsAPrimeB);

            %Number of Signals sent by Alice in test rounds
            numSignalsAlice = numel(probSignalConTest);

            %unpack choi states and outputs in test rounds and 
            % gradient of choi states and gradient outputs in test rounds
            [choiEveAPrimeB,rhoCTestResults]    = unpackVecFWQes(vecFWQes,      dimsAPrimeB,numObs,numTestInt);
            [gradChoiEveAPrimeB,gradRhoC]       = unpackVecFWQes(vecFWQesGrad,  dimsAPrimeB,numObs,numTestInt);

            %RhoCTest reshaped into correct form for optimization
            %size(numobs,numint)
            rhoCTestResults = reshape(rhoCTestResults,[],numTestInt);

            %% CVX
            cvx_begin
            cvx_solver(convertStringsToChars(options.cvxSolver));
            cvx_precision(convertStringsToChars(options.cvxPrecision));
            cvx_quiet(options.verboseLevel<2);
            
            %% Define variables
            deltaRhoC = variable(sprintf('deltaRhoC(%d,%d)',numObs,numTestInt)); 
            YieldMat = variable(sprintf('YieldMat(%d,%d)',numObs,photonCutoff+1), 'nonnegative');
            
            %Sum of current rhoCTestResults and deltaRhoC is more
            %convenient still needs to be nonnegative
            rhoCTestMat = deltaRhoC + rhoCTestResults;            
            rhoCTestMat == nonnegative([numObs,numTestInt]);            

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
                minimize(1);
            else
                tempChoi = cellfun(@(x,y) real(trace(x*y)),gradChoiEveAPrimeB,deltaChoiBlocks,"UniformOutput",false);           
                minimize(gradRhoC.'*reshape(deltaRhoC,[],1) + sum([tempChoi{:}]));
            end
            
            %% constraints
            
            % constraints from test states and channel
            rhoABTests = cell(size(testStatesAAPrime));

            %For squashing we also need states from generation rounds
            rhoABGens = cell(size(genStatesAAprime));

            for index = 1:numBlocks
                %Find dimensions of current block
                currentdimA = dimsA(index);
                currentdimAPrime = dimsAPrime(index);
                currentdimAPrimeB = dimsAPrimeB(index);
                currentdimB = dimsB(index);
                choiTemp = deltaChoiBlocks{index} + choiEveAPrimeB{index};

                %test state of current block
                currentTestStateAAPrime = testStatesAAPrime{index};

                %gen state of current block (only needed for squashing)
                currentGenStateAAPrime = genStatesAAprime{index};

                % Eve's channels are CPTP
                choiTemp == hermitian_semidefinite(double(currentdimAPrimeB));
                PartialTrace(choiTemp, 2,[currentdimAPrime,currentdimB]) == eye(currentdimAPrime);

                % Calculate rhoAB conditioned on test for each block.
                % This is needed for the decoy constraints.
                rhoABTests{index} = PartialMap(currentTestStateAAPrime,...
                    choiTemp,2,[currentdimA,currentdimAPrime]);

                %if squashing constraints are supplied calculate rho
                %conditioned on gen rounds
                if ~isempty(squashingConsGen)
                    rhoABGens{index} = PartialMap(currentGenStateAAPrime,...
                    choiTemp,2,[currentdimA,currentdimAPrime]);
                end
            end

            %Decoy constraints
            decoyConstraintsIntImperfect(cvx_problem, YieldMat, rhoCTestMat, rhoABTests,...
                probDistPhotonConMuTestLower, probDistPhotonConMuTestUpper,probRemaining, probSignalConTest, testCons, blockPhotonNum, linConTol)
            
            % rhoCTestMat is a probability distribution in each column
            abs(sum(rhoCTestMat,1) - ones(1,numTestInt)) <= linConTol

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
            
            %Reshape deltaChoiBlocks into vectors
            deltaChoiBlocks = cell2mat(cellfun(@(x) x(:),deltaChoiBlocks, "uniformoutput", false).');

            %Output overall delta vector
            deltaVecFW = [deltaChoiBlocks;reshape(deltaRhoC,[],1)];

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

            %store tolerance of solution
            tol = cvx_slvtol;
        end

        function [vecQesInit,exitFlag,statusMessage]  = initPoint(vecFW,testStatesAAPrime,genStatesAAprime,dimsAPrime,...
                testCons,squashingConsTest,squashingConsGen,probDistPhotonConMuTest,decoyProbs,probSignalConTest,blockPhotonNum,options)
            arguments
                vecFW (:,1) double
                testStatesAAPrime (1,:) cell {mustBeNonempty,mustBeCellOf(testStatesAAPrime,"double")}% cell of (:,:) double 
                genStatesAAprime (1,:) cell {mustBeNonempty,mustBeCellOf(genStatesAAprime,"double")}
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
                testCons (:,1) EqualityConstraintDecoy
                squashingConsTest (:,1) InequalityConstraint
                squashingConsGen (:,1) InequalityConstraint
                probDistPhotonConMuTest (:,:) double {mustBePositive}
                decoyProbs (:,1) double {mustBeProbDist}
                probSignalConTest (:,1) double {mustBeProbDist}                
                blockPhotonNum (:,1) uint64     
                options (1,1) struct
            end

            %Linear constraint tolerance
            linConTol = options.linearConstraintTolerance;

            %Dimensions of the Choi matrices
            dimsAB = size(testCons(1).operator,1);
            dimsAAPrime = cellfun(@(x) size(x,1), testStatesAAPrime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;

            %Photon Cutoff
            photonCutoff = size(probDistPhotonConMuTest,1)-1;

            %Number of observation
            numObs = numel(testCons);

            %Number of intensities in test rounds
            numTestInt = numel(decoyProbs);

            %Number of Blocks
            numBlocks = numel(dimsAPrimeB);

            %Number of Signals sent by Alice in test rounds
            numSignalsAlice = numel(probSignalConTest);

            % unpack Choi states
            [~,choiEveAPrimeB] = unpackVecFW(vecFW,dimsAPrimeB);

            % states in test rounds
            rhoABTests = cell(size(testStatesAAPrime));

            %For squashing we also need states from generation rounds
            rhoABGens = cell(size(genStatesAAprime));

            for index = 1:numBlocks
                %Find dimensions of current block
                currentdimA = dimsA(index);
                currentdimAPrime = dimsAPrime(index);
                choiTemp = choiEveAPrimeB{index};                

                %test state of current block
                currentTestStateAAPrime = testStatesAAPrime{index};

                %gen state of current block (only needed for squashing)
                currentGenStateAAPrime = genStatesAAprime{index};

                % Calculate rhoAB conditioned on test for each block.
                % This is needed for the decoy constraints.
                rhoABTests{index} = PartialMap(currentTestStateAAPrime,...
                    choiTemp,2,[currentdimA,currentdimAPrime]);

                %if squashing constraints are supplied calculate rho
                %conditioned on gen rounds
                if ~isempty(squashingConsGen)
                    rhoABGens{index} = PartialMap(currentGenStateAAPrime,...
                    choiTemp,2,[currentdimA,currentdimAPrime]);
                end
            end

            %% CVX
            cvx_begin
            cvx_solver(convertStringsToChars(options.cvxSolver));
            cvx_precision(convertStringsToChars(options.cvxPrecision));
            cvx_quiet(options.verboseLevel<2);
            
            %% Define variables
            rhoCTestMat = variable(sprintf('rhoCTestMat(%d,%d)',numObs,numTestInt), 'nonnegative');
            YieldMat = variable(sprintf('YieldMat(%d,%d)',numObs,photonCutoff+1), 'nonnegative');


            %% objective
            % minimize(1)
            minimize( norm(reshape(rhoCTestMat,[],1) - reshape(cat(1,testCons.vector),[],1), 1))
            
            %% constraints

            %Decoy constraints
            decoyConstraints(cvx_problem, YieldMat, rhoCTestMat, rhoABTests,...
                probDistPhotonConMuTest, probSignalConTest, testCons, blockPhotonNum, linConTol);

            % Squashing constraints (apply only if they exist)
            if ~isempty(squashingConsTest) || ~isempty(squashingConsGen)
                squashingConstraints(cvx_problem, rhoABTests, rhoABGens, ...
                    squashingConsTest, squashingConsGen, numBlocks, linConTol);
            end            

            cvx_end
            
            %Output vector containing vectorized Choi state and
            %corresponding rhoCtest
            vecQesInit = [vecFW;reshape(rhoCTestMat,[],1)];

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

        function [vecQesInit,exitFlag,statusMessage] = initPointIntImperfect(vecFW,testStatesAAPrime,genStatesAAprime,dimsAPrime,...
                testCons,squashingConsTest,squashingConsGen,probDistPhotonConMuTestLower,probDistPhotonConMuTestUpper,probRemaining,...
                decoyProbs,probSignalConTest,blockPhotonNum,options)
            arguments
                vecFW (:,1) double
                testStatesAAPrime (1,:) cell {mustBeNonempty,mustBeCellOf(testStatesAAPrime,"double")}% cell of (:,:) double 
                genStatesAAprime (1,:) cell {mustBeNonempty,mustBeCellOf(genStatesAAprime,"double")}
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
                testCons (:,1) EqualityConstraintDecoy
                squashingConsTest (:,1) InequalityConstraint
                squashingConsGen (:,1) InequalityConstraint
                probDistPhotonConMuTestLower (:,:) double {mustBePositive}
                probDistPhotonConMuTestUpper (:,:) double {mustBePositive}
                probRemaining (1,:) double {mustBeNonnegative}
                decoyProbs (:,1) double {mustBeProbDist}
                probSignalConTest (:,1) double {mustBeProbDist}                
                blockPhotonNum (:,1) uint64     
                options (1,1) struct
            end

            %Linear constraint tolerance
            linConTol = options.linearConstraintTolerance;

            %Dimensions of the Choi matrices
            dimsAB = size(testCons(1).operator,1);
            dimsAAPrime = cellfun(@(x) size(x,1), testStatesAAPrime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;

            %Photon Cutoff
            photonCutoff = size(probDistPhotonConMuTestLower,1)-1;

            %Number of observation
            numObs = numel(testCons);

            %Number of intensities in test rounds
            numTestInt = numel(decoyProbs);

            %Number of Blocks
            numBlocks = numel(dimsAPrimeB);

            %Number of Signals sent by Alice in test rounds
            numSignalsAlice = numel(probSignalConTest);

            % unpack Choi states
            [~,choiEveAPrimeB] = unpackVecFW(vecFW,dimsAPrimeB);

            % states in test rounds
            rhoABTests = cell(size(testStatesAAPrime));

            %For squashing we also need states from generation rounds
            rhoABGens = cell(size(genStatesAAprime));

            for index = 1:numBlocks
                %Find dimensions of current block
                currentdimA = dimsA(index);
                currentdimAPrime = dimsAPrime(index);
                choiTemp = choiEveAPrimeB{index};                

                %test state of current block
                currentTestStateAAPrime = testStatesAAPrime{index};

                %gen state of current block (only needed for squashing)
                currentGenStateAAPrime = genStatesAAprime{index};

                % Calculate rhoAB conditioned on test for each block.
                % This is needed for the decoy constraints.
                rhoABTests{index} = PartialMap(currentTestStateAAPrime,...
                    choiTemp,2,[currentdimA,currentdimAPrime]);

                %if squashing constraints are supplied calculate rho
                %conditioned on gen rounds
                if ~isempty(squashingConsGen)
                    rhoABGens{index} = PartialMap(currentGenStateAAPrime,...
                    choiTemp,2,[currentdimA,currentdimAPrime]);
                end
            end

            %% CVX
            cvx_begin
            cvx_solver(convertStringsToChars(options.cvxSolver));
            cvx_precision(convertStringsToChars(options.cvxPrecision));
            cvx_quiet(options.verboseLevel<2);
            
            %% Define variables
            rhoCTestMat = variable(sprintf('rhoCTestMat(%d,%d)',numObs,numTestInt), 'nonnegative');
            YieldMat = variable(sprintf('YieldMat(%d,%d)',numObs,photonCutoff+1), 'nonnegative');


            %% objective
            % minimize(1)
            minimize( norm(reshape(rhoCTestMat,[],1) - reshape(cat(1,testCons.vector),[],1), 1))
            
            %% constraints

            %Decoy constraints
            decoyConstraintsIntImperfect(cvx_problem, YieldMat, rhoCTestMat, rhoABTests,...
                probDistPhotonConMuTestLower, probDistPhotonConMuTestUpper, probRemaining, probSignalConTest, testCons, blockPhotonNum, linConTol)
            
            % rhoCTestMat is a probability distribution in each column
            abs(sum(rhoCTestMat,1) - ones(1,numTestInt)) <= linConTol

            % Squashing constraints (apply only if they exist)
            if ~isempty(squashingConsTest) || ~isempty(squashingConsGen)
                squashingConstraints(cvx_problem, rhoABTests, rhoABGens, ...
                    squashingConsTest, squashingConsGen, numBlocks, linConTol);
            end            

            cvx_end
            
            %Output vector containing vectorized Choi state and
            %corresponding rhoCtest
            vecQesInit = [vecFW;reshape(rhoCTestMat,[],1)];

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

        function perturbedVec = perturbEigvalsChoi(vecFWQES,dimsAPrimeB,includeRhoC,minEigenvalue)
            arguments
                vecFWQES (:,1) double
                dimsAPrimeB (1,:) double {mustBeInteger,mustBePositive}
                includeRhoC (1,1) logical
                minEigenvalue (1,1) double {mustBeGreaterThanOrEqual(minEigenvalue,0)} = 1e-12;
            end          

            %% parsing Choi matrices from vecQES to matrices
            totalindex = 1;
            choiEveAPrimeB = cell(1,numel(dimsAPrimeB));

            for index = 1:numel(dimsAPrimeB)
                %Get dimensions of A'B of current block
                currentDimsAPrimeB = dimsAPrimeB(index);
                
                %Get choi state
                choiEveAPrimeB{index} = reshape(vecFWQES(totalindex:totalindex+currentDimsAPrimeB^2-1), ...
                    [currentDimsAPrimeB,currentDimsAPrimeB]);

                totalindex = totalindex + currentDimsAPrimeB^2;
            end

            %% perturb Choi states
            %find necessary perturbations
            pertubation = cellfun(@(x) perturbationChannelEpsilon(x,minEigenvalue), choiEveAPrimeB);
            
            %perturb
            perturbedChoi = arrayfun(@(x) perturbationChannel(choiEveAPrimeB{x},pertubation(x)), 1:numel(choiEveAPrimeB),"UniformOutput",false);
            
            %Reshape back to vectors
            perturbedChoi = cell2mat(cellfun(@(x) x(:),perturbedChoi, "uniformoutput", false).');
            
            %return pertubed Choi states and remainder of vector
            if includeRhoC
                perturbedVec = [perturbedChoi;vecFWQES(numel(perturbedChoi)+1:end)];                
            else
                perturbedVec = [perturbedChoi;vecFWQES(end)];
            end
        end
    end
end

function penaltyVal = penaltyTermVariableLengthDecoy(cvx_problem,rhoCTestresults,testProb,decoyProbs,lambda)
    arguments
        cvx_problem
        rhoCTestresults (:,:) % double or CVX            
        testProb (1,1) double
        decoyProbs (:,1) double
        lambda (:,:) % double or CVX
    end
    
    % rho_C = (gamma* p(mu_i|test) * rho_C|Mu,Test
    rhoCTestMatJoint = reshape(testProb*rhoCTestresults*diag(decoyProbs),[],1);       
    
    %D(lambda||rhoC)
    penaltyVal = sum(rel_entr(lambda(1:end-1),rhoCTestMatJoint))/log(2)...
                + rel_entr(lambda(end),(1-testProb))/log(2);
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
        testCons
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
    YieldMat*probDistPhotonConMuTest <= rhoCTestMat + linConTol;
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
        squashingConsTest
        squashingConsGen
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
                    squashingConsTest(indexCons).upperBound;% + linConTol;
            end
            if ~isinf(squashingConsTest(indexCons).lowerBound) %lower bound
                real(trace(squashingConsTest(indexCons).operator'*(rhoABTests{indexBlock}))) >= ... 
                    squashingConsTest(indexCons).lowerBound;% - linConTol;
            end
        end
    end

    % Create squashing constraints for gen rounds
    for indexBlock = 1:numBlocks
        %iterate over constraints of Gen/Test block
        for indexCons = 1:numel(squashingConsGen)
            if ~isinf(squashingConsGen(indexCons).upperBound) %upper bound
                real(trace(squashingConsGen(indexCons).operator'*(rhoABGens{indexBlock}))) <= ...
                    squashingConsGen(indexCons).upperBound;% + linConTol;
            end
            if ~isinf(squashingConsGen(indexCons).lowerBound) %lower bound
                real(trace(squashingConsGen(indexCons).operator'*(rhoABGens{indexBlock}))) >= ... 
                    squashingConsGen(indexCons).lowerBound;% - linConTol;
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

function [choiMatsEveAPrimeB,rhoCTestResults] = unpackVecFWQes(vecFWQes,dimsAPrimeB,numObs,numTestInt)
    %number of blocks
    numBlocks = numel(dimsAPrimeB);    

    %parsing Choi matrices from vecFWQes to matrices
    totalindex = 1;
    choiMatsEveAPrimeB = cell(1,numBlocks);

    for index = 1:numBlocks
        %Get dimensions of A'B of current block
        currentDimsAPrimeB = dimsAPrimeB(index);
        
        %Get choi state
        choiMatsEveAPrimeB{index} = reshape(vecFWQes(totalindex:totalindex+currentDimsAPrimeB^2-1), ...
            [currentDimsAPrimeB,currentDimsAPrimeB]);

        totalindex = totalindex + currentDimsAPrimeB^2;
    end

    %Extract outputs in test rounds from vecFWQes
    rhoCTestResults = vecFWQes(end-numObs*numTestInt+1:end);
end