classdef QESFinderChoiBlock
    methods (Static)
        function [gradQES,tol,exitFlag] = getGradQES(vecFW, vecFWGrad, testStatesAAPrime,...
                genStatesAAprime, dimsAPrime, testProb, testCons,...
                squashingConsTest, squashingConsGen, probBlockConTest, ...
                probSignalConTestAndBlock, probRemaining, epsilonProb, epsilonBlock, options)
            arguments
                vecFW (:,1) double
                vecFWGrad (:,1) double
                testStatesAAPrime (1,:) cell {mustBeNonempty,mustBeCellOf(testStatesAAPrime,"double")}% cell of (:,:) double
                genStatesAAprime (1,:) cell {mustBeNonempty,mustBeCellOf(genStatesAAprime,"double")}
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
                testProb (1,1) double {mustBeInRange(testProb, 0, 1)}
                testCons (:,1) EqualityConstraintChoiBlock
                squashingConsTest (:,:) InequalityConstraint
                squashingConsGen (:,:) InequalityConstraint
                probBlockConTest (1,:) double {mustBeNonnegative}
                probSignalConTestAndBlock (:,:) double {mustBeNonnegative}
                probRemaining (:,1) double {mustBeNonnegative,mustBeInRange(probRemaining, 0, 1)}
                epsilonProb (1,1) double {mustBeNonnegative}
                epsilonBlock (1,:) double {mustBeNonnegative}
                options (1,1) struct
            end

            %Linear constraint tolerance
            linConTol = 1/2*options.linearConstraintTolerance;

            %Dimensions of the Choi matrices
            dimsAB = testCons(1).rhoDim;
            dimsAAPrime = cellfun(@(x) size(x,1), testStatesAAPrime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;

            %Number of observation
            numObs = numel(testCons);

            %Number of Blocks
            numBlocks = numel(dimsAPrimeB);


            %Define accept point
            % q = (gamma* p((a,mu),b|test), 1-gamma)^T
            acceptPoint = [reshape(testProb*cat(1,testCons.scalar),[],1); 1-testProb];
                       
            % unpack all the values from vecFW and vecFWGrad
            [phi,choiEveAPrimeB]         = unpackVecFW(vecFW,    dimsAPrimeB);
            [gradPhi,gradChoiEveAPrimeB] = unpackVecFW(vecFWGrad,dimsAPrimeB);

            %% CVX
            cvx_begin
            cvx_solver(convertStringsToChars(options.cvxSolver));
            cvx_precision(convertStringsToChars(options.cvxPrecision));
            cvx_quiet(options.verboseLevel<2);
            
            lambda = variable('lambda(numObs+1,1)');
            dualLambda = dual('variable','dualLambda'); %declares dual variable
            deltaPhi = variable('deltaPhi(1,1)');            
            rhoCTestMat = variable(sprintf('rhoCTestMat(%d)',numObs), 'nonnegative');

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
            phi+deltaPhi >= penaltyTermVariableLengthChoiBlock(cvx_problem,rhoCTestMat,testProb,lambda);
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
            
            % block constraints
            blockConstraints(cvx_problem, rhoCTestMat, rhoABTests,...
                    probBlockConTest, probSignalConTestAndBlock, probRemaining, testCons,...
                    epsilonProb, epsilonBlock, linConTol);
            
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
    
        function val = funcFW(vecFWQes,genStatesAAprime,keyProj,krausOps,testCons, ...
                testProb,gradQES,alpha,dimsAPrime,probBlockConGen,perturbation,perturbationAff)
            arguments
                vecFWQes (:,1) double
                genStatesAAprime (1,:) cell {mustBeNonempty,mustBeCellOf(genStatesAAprime,"double")}
                keyProj (1,:) cell
                krausOps (1,:) cell
                testCons (:,1) EqualityConstraintChoiBlock
                testProb (1,1) double
                gradQES (:,1) double
                alpha (1,1) double
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
                probBlockConGen (:,1) double {mustBeNonnegative}
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
                perturbationAff (1,1) double {mustBeInRange(perturbationAff,0,1)}
            end

            %Dimensions of the Choi matrices
            dimsAB = testCons(1).rhoDim;
            dimsAAPrime = cellfun(@(x) size(x,1), genStatesAAprime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;

            %Number of observation
            numObs = numel(testCons);

            %Number of Blocks
            numBlocks = numel(dimsAPrimeB);
            
            %unpack choi states and outputs in test rounds
            [choiEveAPrimeB,rhoCTestResults] = unpackVecFWQes(vecFWQes,dimsAPrimeB,numObs);         
                        
            %Find corresponding states rhoAB
            rhoAB = cell(1,numel(probBlockConGen));
            for index=1:numel(probBlockConGen)                
                rhoOp = PartialMap(genStatesAAprime{index},choiEveAPrimeB{index},...
                    2,[dimsA(index),dimsAPrime(index)]);
                rhoAB{index} = 1/2*(rhoOp + rhoOp'); %remove small non-hermitian parts
            end
            
            %Concatenate rhoAB and keyProj and krausOps as repeating arguments,
            %i.e. (a1,b1,a2,b2,...)
            repeatingArgs = [rhoAB;keyProj;krausOps];                        

            %rhoCtest*2^((alpha-1)*f)         
            sumTraceConTest = (rhoCTestResults.')*2.^((alpha-1)*gradQES(1:end-1));
            
            %Value of QES entropy
            val = -1/(alpha-1)*log2((1-testProb)*2^(-(alpha-1)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,probBlockConGen,repeatingArgs{:}))*2^((alpha-1)*gradQES(end)) ...
                + testProb*sumTraceConTest);
        end

        function grad = gradFuncFW(vecFWQes,genStatesAAprime,keyProj,krausOps,testCons, ...
                testProb,gradQES,alpha,dimsAPrime,probBlockConGen,perturbation,perturbationAff)
            arguments
                vecFWQes (:,1) double
                genStatesAAprime (1,:) cell {mustBeCellOf(genStatesAAprime,"double")}
                keyProj (1,:) cell {mustBeCellOf(keyProj,"cell")}
                krausOps (1,:) cell {mustBeCellOf(krausOps,"cell")}
                testCons (:,1) EqualityConstraintChoiBlock
                testProb (1,1) double
                gradQES (:,1) double
                alpha (1,1) double
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
                probBlockConGen (:,1) double {mustBeNonnegative}
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
                perturbationAff(1,1) double {mustBeInRange(perturbationAff,0,1)} 
            end

            %Dimensions of the Choi matrices
            dimsAB = testCons(1).rhoDim;
            dimsAAPrime = cellfun(@(x) size(x,1), genStatesAAprime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;

            %Number of observation
            numObs = numel(testCons);

            %Number of Blocks
            numBlocks = numel(dimsAPrimeB);
            
            %unpack choi states and outputs in test rounds
            [choiEveAPrimeB,rhoCTestResults] = unpackVecFWQes(vecFWQes,dimsAPrimeB,numObs);         
                        
            %Find corresponding states rhoAB
            rhoAB = cell(1,numel(probBlockConGen));
            for index=1:numel(probBlockConGen)                
                rhoOp = PartialMap(genStatesAAprime{index},choiEveAPrimeB{index},...
                    2,[dimsA(index),dimsAPrime(index)]);
                rhoAB{index} = 1/2*(rhoOp + rhoOp'); %remove small non-hermitian parts
            end
            
            %Concatenate rhoAB and keyProj and krausOps as repeating arguments,
            %i.e. (a1,b1,a2,b2,...)
            repeatingArgs = [rhoAB;keyProj;krausOps];                        

            %rhoCtest*2^((alpha-1)*f)         
            sumTraceConTest = (rhoCTestResults.')*2.^((alpha-1)*gradQES(1:end-1));

            %Calculate inner term of log
            innerTerm = (1-testProb)*2^(-(alpha-1)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,probBlockConGen,repeatingArgs{:}))*2^((alpha-1)*gradQES(end)) ...
                + testProb*sumTraceConTest;
            
            %Gradient by RhoCTest
            gradByRhoCTest = -testProb/(alpha-1)/log(2)*1/innerTerm*2.^((alpha-1)*gradQES(1:end-1));

            %Gradient by Choi
            gradByChoi = RenyiTildeDown.gradConEnt(alpha,perturbation,perturbationAff,probBlockConGen,repeatingArgs{:});
            for index=1:numel(probBlockConGen)
                currentdimA = size(genStatesAAprime{index},1)/dimsAPrime(index);
                
                % The choi matrix must be ordered input system tensor output
                % system.
                genChoiAPrimeA = Swap(genStatesAAprime{index},[1,2],[currentdimA,dimsAPrime(index)]);        
                tempDualMap = DualMap(genChoiAPrimeA,[dimsAPrime(index),currentdimA]);

                gradConEntByChoi = PartialMap(gradByChoi{index},tempDualMap,1,[currentdimA,dimsB(index)]);
                
                gradByChoi{index} = (1-testProb)/innerTerm*2^(-(alpha-1)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,probBlockConGen,repeatingArgs{:}))...
                *2^((alpha-1)*gradQES(end))*gradConEntByChoi;
            end
            %Reshape each gradient in the cell gradbyChoi to be linearized
            gradByChoi = cell2mat(cellfun(@(x) x(:),gradByChoi, "uniformoutput", false));

            %Final gradient
            grad = [gradByChoi(:);gradByRhoCTest];
        end

        function [deltaVecFW,exitFlag,statusMessage,tol] = subproblemQes(vecFWQes,vecFWQesGrad,...
                testStatesAAPrime,genStatesAAprime,dimsAPrime,testCons,squashingConsTest,squashingConsGen, ...
                probBlockConTest,probSignalConTestAndBlock,probRemaining, epsilonProb, epsilonBlock,options,feasibleMode)
            arguments
                vecFWQes (:,1) double
                vecFWQesGrad (:,1) double
                testStatesAAPrime (1,:) cell {mustBeNonempty,mustBeCellOf(testStatesAAPrime,"double")}% cell of (:,:) double
                genStatesAAprime (1,:) cell {mustBeNonempty,mustBeCellOf(genStatesAAprime,"double")}
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
                testCons (:,1) EqualityConstraintChoiBlock
                squashingConsTest (:,:) InequalityConstraint
                squashingConsGen (:,:) InequalityConstraint
                probBlockConTest (1,:) double {mustBeNonnegative}
                probSignalConTestAndBlock (:,:) double {mustBeNonnegative}
                probRemaining (:,1) double {mustBeNonnegative,mustBeInRange(probRemaining, 0, 1)}
                epsilonProb (1,1) double {mustBeNonnegative}
                epsilonBlock (1,:) double {mustBeNonnegative}
                options (1,1) struct
                feasibleMode (1,1) logical = false;
            end

            %Linear constraint tolerance
            linConTol = options.linearConstraintTolerance;

            %Dimensions of the Choi matrices
            dimsAB = testCons(1).rhoDim;
            dimsAAPrime = cellfun(@(x) size(x,1), genStatesAAprime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;

            %Number of observation
            numObs = numel(testCons);

            %Number of Blocks
            numBlocks = numel(dimsAPrimeB);
            
            %unpack choi states and outputs in test rounds and 
            % gradient of choi states and gradient outputs in test rounds
            [choiEveAPrimeB,rhoCTestResults]    = unpackVecFWQes(vecFWQes,      dimsAPrimeB,numObs);
            [gradChoiEveAPrimeB,gradRhoC]       = unpackVecFWQes(vecFWQesGrad,  dimsAPrimeB,numObs);

            %% CVX
            cvx_begin
            cvx_solver(convertStringsToChars(options.cvxSolver));
            cvx_precision(convertStringsToChars(options.cvxPrecision));
            cvx_quiet(options.verboseLevel<2);
            
            %% Define variables
            deltaRhoC = variable(sprintf('deltaRhoC(%d)',numObs)); 
                        
            %Sum of current rhoCTestResults and deltaRhoC is more
            %convenient still needs to be nonnegative
            rhoCTestMat = deltaRhoC + rhoCTestResults;
            rhoCTestMat == nonnegative([numObs,1]);            

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

            % block constraints
            blockConstraints(cvx_problem, rhoCTestMat, rhoABTests,...
                    probBlockConTest, probSignalConTestAndBlock, probRemaining, testCons, ...
                    epsilonProb, epsilonBlock, linConTol);
            
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

        function [deltaVecFW,exitFlag,statusMessage] = subproblemQesYalmip(vecFWQes,vecFWQesGrad,...
            testStatesAAPrime,genStatesAAprime,dimsAPrime,testCons,squashingConsTest,squashingConsGen, ...
            probDistPhotonConMuTest,decoyProbs,probSignalConTest,blockPhotonNum,options,feasibleMode)
            arguments
                vecFWQes (:,1) double
                vecFWQesGrad (:,1) double
                testStatesAAPrime (1,:) cell {mustBeNonempty,mustBeCellOf(testStatesAAPrime,"double")}% cell of (:,:) double 
                genStatesAAprime (1,:) cell {mustBeNonempty,mustBeCellOf(genStatesAAprime,"double")}
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
                testCons (:,1) EqualityConstraintChoiBlock
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
    
            %unpack choi states and outputs in test rounds and 
            % gradient of choi states and gradient outputs in test rounds
            [choiEveAPrimeB,rhoCTestResults]    = unpackVecFWQes(vecFWQes,      dimsAPrimeB,numObs,numTestInt);
            [gradChoiEveAPrimeB,gradRhoC]       = unpackVecFWQes(vecFWQesGrad,  dimsAPrimeB,numObs,numTestInt);
    
            %RhoCTest reshaped into correct form for optimization
            %size(numobs,numint)
            rhoCTestResults = reshape(rhoCTestResults,[],numTestInt);
    
            %% Yalmip
                    
            %% Define variables
            deltaRhoC = sdpvar(numObs,numTestInt); 
            YieldMat = sdpvar(numObs,photonCutoff+1);
            
            %Sum of current rhoCTestResults and deltaRhoC is more
            %convenient still needs to be nonnegative
            rhoCTestMat = deltaRhoC + rhoCTestResults;            
            consRhoCPos = [rhoCTestMat >= 0];            
    
            if options.blockDiagonal
                error("Not Ready for this!");
            else
                deltaChoiBlocks = cell(1,numBlocks);
                for index = 1:numBlocks
                    deltaChoiString = sprintf("deltaChoi%1d = sdpvar(%d,%d,'hermitian','complex');",index,dimsAPrimeB(index),dimsAPrimeB(index));
                    eval([deltaChoiString]);
                    deltaChoiBlockString = sprintf("deltaChoiBlocks{%1d} = deltaChoi%d;",index,index);
                    eval([deltaChoiBlockString]);
                end
            end
    
    
            %% objective
            if feasibleMode
                obj = 1;
            else
                tempChoi = cellfun(@(x,y) real(trace(x*y)),gradChoiEveAPrimeB,deltaChoiBlocks,"UniformOutput",false);           
                obj = gradRhoC.'*reshape(deltaRhoC,[],1) + sum([tempChoi{:}]);
            end
    
                   
            %% constraints
            
            % constraints from test states and channel
            rhoABTests = cell(size(testStatesAAPrime));

            %For squashing we also need states from generation rounds
            rhoABGens = cell(size(genStatesAAprime));
            
            consChoiPsd = [];
            consChoiPartialTrace = [];
    
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
                consChoiPsd = [consChoiPsd, choiTemp >= 0];
                consChoiPartialTrace = [consChoiPartialTrace, ...
                PartialTrace(choiTemp, 2,[currentdimAPrime,currentdimB]) == eye(currentdimAPrime)];
    
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
            consDecoy = decoyConstraintsYalmip(YieldMat, rhoCTestMat, rhoABTests,...
                probDistPhotonConMuTest, probSignalConTest, testCons, blockPhotonNum, linConTol);
            
            % Squashing constraints (apply only if they exist)
            consSquash = [];
            if ~isempty(squashingConsTest) || ~isempty(squashingConsGen)
                consSquash = squashingConstraintsYalmip(rhoABTests, rhoABGens, ...
                    squashingConsTest, squashingConsGen, numBlocks, linConTol);
            end          
           
            %% Run optimization
            Constraints = [consRhoCPos, consChoiPsd, consChoiPartialTrace, consDecoy, consSquash];       
            % solverOptions = sdpsettings('solver', 'Mosek', 'verbose', 0,'savesolveroutput',1);
            solverOptions = sdpsettings('solver', 'sdpa-gmp-neos', 'verbose', 0,'savesolveroutput',1);
            diagnostics = optimize(Constraints,obj,solverOptions);
    
            % reconstruct rho from blocks incase it didn't get filled in by CVX
            if options.blockDiagonal
                error("Not Ready for this!");
                % deltaChoi = zeros(dimAPrimeB);
                % deltaChoiBlocks = eval(deltaChoiBlocksString);
                % deltaChoi = directSumWorkArround(deltaChoi, deltaChoiBlocks);
                % deltaChoi = permMat'*deltaChoi*permMat;
            end
            
            %Reconstruct Blocks
            for indexBlock = 1:numBlocks
                deltaChoiBlockString = sprintf('value(deltaChoi%1d);',indexBlock);
                deltaChoiBlocks{indexBlock} = eval(deltaChoiBlockString);
            end
            
            %Reshape deltaChoiBlocks into vectors
            deltaChoiBlocks = cell2mat(cellfun(@(x) x(:),deltaChoiBlocks, "uniformoutput", false).');
    
            %Output overall delta vector
            deltaVecFW = [deltaChoiBlocks;reshape(value(deltaRhoC),[],1)];
    
            %Report status
            statusMessage = string(diagnostics.info);
            statusString = split(statusMessage,' (');
            
            switch statusString(1)
                case "Successfully solved"
                    exitFlag = SubProblemExitFlag.solved;
                % case "Inaccurate/Solved"
                %     exitFlag = SubProblemExitFlag.inaccurate;
                % case "Failed"
                %     exitFlag = SubProblemExitFlag.failed;
                otherwise
                    exitFlag = SubProblemExitFlag.failed;
            end
            
        end

        function [vecQesInit,exitFlag,statusMessage] = initPoint(vecFW,testStatesAAPrime,genStatesAAprime,dimsAPrime,...
                testCons,squashingConsTest,squashingConsGen,probBlockConTest,probSignalConTestAndBlock,probRemaining,...
                epsilonProb, epsilonBlock, options)
            arguments
                vecFW (:,1) double
                testStatesAAPrime (1,:) cell {mustBeNonempty,mustBeCellOf(testStatesAAPrime,"double")}% cell of (:,:) double
                genStatesAAprime (1,:) cell {mustBeNonempty,mustBeCellOf(genStatesAAprime,"double")}
                dimsAPrime (1,:) double {mustBeInteger,mustBePositive}
                testCons (:,1) EqualityConstraintChoiBlock
                squashingConsTest (:,:) InequalityConstraint
                squashingConsGen (:,:) InequalityConstraint
                probBlockConTest (1,:) double {mustBeNonnegative}
                probSignalConTestAndBlock (:,:) double {mustBeNonnegative}
                probRemaining (:,1) double {mustBeNonnegative,mustBeInRange(probRemaining, 0, 1)}
                epsilonProb (1,1) double {mustBeNonnegative}
                epsilonBlock (1,:) double {mustBeNonnegative}
                options (1,1) struct
            end

            %Linear constraint tolerance
            linConTol = options.linearConstraintTolerance;

            %Dimensions of the Choi matrices
            dimsAB = testCons(1).rhoDim;
            dimsAAPrime = cellfun(@(x) size(x,1), genStatesAAprime);
            dimsA = dimsAAPrime./dimsAPrime;
            dimsAPrimeB = dimsAB./dimsA .* dimsAPrime;
            dimsB = dimsAB./dimsA;

            %Number of observation
            numObs = numel(testCons);

            %Number of Blocks
            numBlocks = numel(dimsAPrimeB);

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
            rhoCTestMat = variable(sprintf('rhoCTestMat(%d)',numObs), 'nonnegative');


            %% objective
            % minimize(1)
            minimize( norm(reshape(rhoCTestMat,[],1) - reshape(cat(1,testCons.scalar),[],1), 1))
            
            %% constraints

            % block constraints
            blockConstraints(cvx_problem, rhoCTestMat, rhoABTests,...
                    probBlockConTest, probSignalConTestAndBlock, probRemaining, testCons,...
                    epsilonProb, epsilonBlock, linConTol);
            
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

function penaltyVal = penaltyTermVariableLengthChoiBlock(cvx_problem,rhoCTestMat,testProb,lambda)
    arguments
        cvx_problem
        rhoCTestMat (:,:) % double or CVX
        testProb (1,1) double
        lambda (:,:) % double or CVX
    end
    
    % rho_C = gamma* rho_C|Test
    rhoCTestMatJoint = reshape(testProb*rhoCTestMat,[],1);
    
    %D(lambda||rhoC)
    penaltyVal = sum(rel_entr(lambda(1:end-1),rhoCTestMatJoint))/log(2)...
                + rel_entr(lambda(end),(1-testProb))/log(2);
end

function blockConstraints(cvx_problem, rhoCTestMat, rhoABTests, probBlockConTest,...
    probSignalConTestAndBlock, probRemaining ,testCons, epsilonProb, epsilonBlock, linConTol)
    arguments
        cvx_problem (1,1)
        rhoCTestMat (:,:)
        rhoABTests
        probBlockConTest (:,1) double {mustBeNonnegative,mustBeNonempty}
        probSignalConTestAndBlock(:,:) double {mustBeNonnegative,mustBeNonempty}
        probRemaining (:,1) double {mustBeNonnegative,mustBeInRange(probRemaining, 0, 1)}                
        testCons (:,1) EqualityConstraintChoiBlock
        epsilonProb (1,1) double {mustBeNonnegative}
        epsilonBlock (1,:) double {mustBeNonnegative}
        linConTol (1,1) double {mustBeNonnegative}
    end
    
    %various sizes
    blockCutoff = size(probBlockConTest,1);
    numSignalsAlice = size(probSignalConTestAndBlock,1);
    numObs = numel(testCons);

    %Construct Tr[Gamma_{a,b} xi_{t,k}(J_k)] for each block k
    for indexBlock = 1:blockCutoff
        for indexCons = 1:numObs
            operatorCell = testCons(indexCons).operatorCell;
            obsMatBlock(indexCons,indexBlock) = ...
                real(trace(operatorCell{indexBlock}*rhoABTests{indexBlock}));
        end
    end

    %Lower bound on probBlockConTest
    probBlockConTestLower = max(probBlockConTest-epsilonProb*ones(blockCutoff,1),0);

    %Upper bound on probBlockConTest
    probBlockConTestUpper = min(probBlockConTest+epsilonProb*ones(blockCutoff,1),1);

    %recast epsilonBlock into same format as obsMatBlock
    epsilonBlockMat = repmat(epsilonBlock,numObs,1);
    
    % equivalent to 0 <= rhoCTestMat_{a,b} - sum_k p(k) Tr[Gamma_{a,b} xi_{t,k}(J_k)] <= ones(numObs,1)*(1-pTot);
    (obsMatBlock - epsilonBlockMat)*probBlockConTestLower <= rhoCTestMat                 + linConTol;
    (obsMatBlock + epsilonBlockMat)*probBlockConTestUpper >= rhoCTestMat - probRemaining - linConTol;
    
    %sum_b Y_n(a,b) = p(a|test,n)
    obsBlockTensor = reshape(obsMatBlock,numSignalsAlice,[],blockCutoff);
    abs(squeeze(sum(obsBlockTensor,2)) - probSignalConTestAndBlock) <= linConTol;
    
    %Y_n(a,b) <= p(a|test,n) Not needed, but sometimes it works better. No clue
    %why.
    obsMatBlock <= repmat(probSignalConTestAndBlock,size(obsBlockTensor,2),1);
    obsMatBlock == nonnegative(size(obsMatBlock));
end

function consDecoy = decoyConstraintsYalmip(YieldMat, rhoCTestMat, rhoABTests,...
        probDistPhotonConMuTest, probSignalConTest, testCons, blockPhotonNum, linConTol)
    arguments
        YieldMat (:,:)
        rhoCTestMat (:,:)
        rhoABTests
        probDistPhotonConMuTest
        probSignalConTest (:,1) double
        testCons
        blockPhotonNum (:,1) uint64
        linConTol (1,1) double {mustBeNonnegative}
    end
    %Store decoy constraints
    consDecoy = [];
    
    photonCutoff = size(probDistPhotonConMuTest,1)-1;
    numSignalsAlice = numel(probSignalConTest);
    numObs = numel(testCons);
    
    %sum of yields is close (<= 1- ptot(mu)) to statistics generated from state
    pTotMu = sum(probDistPhotonConMuTest,1);
    probRemaining = ones(numObs,1)*(1-pTotMu);
    
    % equivalent to 0 <= rhoCTestMat - YieldMat*probDistPhotonConMuTest <= ones(numObs,1)*(1-pTotMu);
    consDecoy = [consDecoy, YieldMat*probDistPhotonConMuTest <= rhoCTestMat + linConTol ];
    consDecoy= [consDecoy, YieldMat*probDistPhotonConMuTest >= rhoCTestMat - probRemaining - linConTol ];
    
    %sum_b Y_n(a,b) = p(a|test,n)
    YieldTensor = reshape(YieldMat,numSignalsAlice,[],photonCutoff+1);
    consDecoy = [consDecoy, abs(squeeze(sum(YieldTensor,2)) - probSignalConTest*ones(1,photonCutoff+1)) <= linConTol];
    
    %Y_n(a,b) <= p(a|test,n) Not needed, but sometimes it works better. No clue
    %why.
    consDecoy = [consDecoy, YieldTensor <= repmat(probSignalConTest,1,size(YieldTensor,2),size(YieldTensor,3))];
    
    
    % Trace constraints on e.g. single photons
    for index=1:numel(blockPhotonNum)
        rhoABTestTemp = rhoABTests{index};
        for indexCons = 1:numObs
            consDecoy = [consDecoy, YieldMat(indexCons,double(blockPhotonNum(index)+1)) == ... 
                real(trace(testCons(indexCons).operator*rhoABTestTemp))];
        end
    end
end

function squashingConstraints(cvx_problem, rhoABTests, rhoABGens, ...
        squashingConsTest, squashingConsGen, numBlocks, linConTol)
    arguments
        cvx_problem (1,1)            
        rhoABTests
        rhoABGens
        squashingConsTest (:,:) InequalityConstraint
        squashingConsGen (:,:) InequalityConstraint
        numBlocks (1,1) uint64
        linConTol (1,1) double {mustBeNonnegative}
    end
    % Squashing constraints

    % Create squashing constraints for test rounds
    for indexBlock = 1:numBlocks
        %iterate over constraints of Gen/Test block
        for indexCons = 1:size(squashingConsTest,1)            
            if ~isinf(squashingConsTest(indexCons,indexBlock).upperBound) %upper bound
                real(trace(squashingConsTest(indexCons,indexBlock).operator'*(rhoABTests{indexBlock}))) <= ...
                    squashingConsTest(indexCons,indexBlock).upperBound + linConTol;
            end
            if ~isinf(squashingConsTest(indexCons,indexBlock).lowerBound) %lower bound
                real(trace(squashingConsTest(indexCons,indexBlock).operator'*(rhoABTests{indexBlock}))) >= ... 
                    squashingConsTest(indexCons,indexBlock).lowerBound - linConTol;
            end
        end
    end

    % Create squashing constraints for gen rounds
    for indexBlock = 1:numBlocks
        %iterate over constraints of Gen/Test block
        for indexCons = 1:size(squashingConsGen,1)
            if ~isinf(squashingConsGen(indexCons,indexBlock).upperBound) %upper bound
                real(trace(squashingConsGen(indexCons,indexBlock).operator'*(rhoABGens{indexBlock}))) <= ...
                    squashingConsGen(indexCons,indexBlock).upperBound + linConTol;
            end
            if ~isinf(squashingConsGen(indexCons,indexBlock).lowerBound) %lower bound
                real(trace(squashingConsGen(indexCons,indexBlock).operator'*(rhoABGens{indexBlock}))) >= ... 
                    squashingConsGen(indexCons,indexBlock).lowerBound - linConTol;
            end
        end
    end
end

function consSquash = squashingConstraintsYalmip(rhoABTests, rhoABGens, ...
        squashingConsTest, squashingConsGen, numBlocks, linConTol)
    arguments      
        rhoABTests
        rhoABGens
        squashingConsTest
        squashingConsGen
        numBlocks (1,1) uint64
        linConTol (1,1) double {mustBeNonnegative}
    end
    % Squashing constraints
    consSquash = [];

    % Create squashing constraints for test rounds
    for indexBlock = 1:numBlocks
        %iterate over constraints of Gen/Test block
        for indexCons = 1:numel(squashingConsTest)
            if ~isinf(squashingConsTest(indexCons).upperBound) %upper bound
                consSquash = [consSquash,real(trace(squashingConsTest(indexCons).operator'*(rhoABTests{indexBlock}))) <= ...
                    squashingConsTest(indexCons).upperBound + linConTol];
            end
            if ~isinf(squashingConsTest(indexCons).lowerBound) %lower bound
                consSquash = [consSquash,real(trace(squashingConsTest(indexCons).operator'*(rhoABTests{indexBlock}))) >= ... 
                    squashingConsTest(indexCons).lowerBound - linConTol];
            end
        end
    end

    % Create squashing constraints for gen rounds
    for indexBlock = 1:numBlocks
        %iterate over constraints of Gen/Test block
        for indexCons = 1:numel(squashingConsGen)
            if ~isinf(squashingConsGen(indexCons).upperBound) %upper bound
                consSquash = [consSquash,real(trace(squashingConsGen(indexCons).operator'*(rhoABGens{indexBlock}))) <= ...
                    squashingConsGen(indexCons).upperBound + linConTol];
            end
            if ~isinf(squashingConsGen(indexCons).lowerBound) %lower bound
                consSquash = [consSquash,real(trace(squashingConsGen(indexCons).operator'*(rhoABGens{indexBlock}))) >= ... 
                    squashingConsGen(indexCons).lowerBound - linConTol];
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

function [choiMatsEveAPrimeB,rhoCTestResults] = unpackVecFWQes(vecFWQes,dimsAPrimeB,numObs)
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
    rhoCTestResults = vecFWQes(end-numObs+1:end);
end