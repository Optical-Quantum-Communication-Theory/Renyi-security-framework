classdef QESFinder
    methods (Static)
        function [gradQES,tol] = getGradQES(vecFW,vecFWGrad,...
                testStateAAPrime,dimAPrime,testProb,testCons,options)
            arguments
                vecFW (:,1) double
                vecFWGrad (:,1) double
                testStateAAPrime (:,:) double
                dimAPrime (1,1) uint64
                testProb (1,1) double {mustBeInRange(testProb, 0, 1)}
                testCons (:,1) EqualityConstraint
                options (1,1) struct
            end

            %Number of observations
            numObs = numel(testCons);
            
            %phi and gradient
            phi = vecFW(end);
            gradPhi = vecFWGrad(end);
            
            %Choi state and gradient
            dimAPrimeB = sqrt(numel(vecFW)-1);
            choiEveAPrimeB = reshape(vecFW(1:end-1),[dimAPrimeB,dimAPrimeB]);            
            gradChoiEveAPrimeB = reshape(vecFWGrad(1:end-1),[dimAPrimeB,dimAPrimeB]);
            
            %Additional dimensions required
            dimB = dimAPrimeB/dimAPrime;
            dimAAPrime = size(testStateAAPrime,1);
            dimA = dimAAPrime/dimAPrime;


            %% CVX
            cvx_begin sdp
            cvx_solver(convertStringsToChars(options.cvxSolver));
            cvx_precision(convertStringsToChars(options.cvxPrecision));
            cvx_quiet(options.verboseLevel<2);


            deltaPhi = variable('deltaPhi(1,1)');

            lambda = variable('lambda(1,numObs+1)');

            dualLambda = dual('variable','dualLambda'); %declares dual variable

            if options.blockDiagonal

                %% block diagonal set up which cant be placed in a separate function
                permMat = options.blockP;
                newDims = options.newDims.';

                deltaChoiStrings = compose("deltaChoi%1d(%d,%d)",(1:numel(newDims)).',newDims,newDims);
                % CVX MUST generates variables in this work space. We have to do each
                % separately.
                for index = 1:numel(deltaChoiStrings)
                    variable(convertStringsToChars(deltaChoiStrings(index)),'hermitian')
                end
                % Get a single cell array to capture those variables.
                deltaChoiBlocksString = "{"+strjoin(compose("deltaChoi%d",1:numel(newDims)),",")+"}";
                deltaChoiBlocks = eval(deltaChoiBlocksString);


                % there is a bug in cvx's blkdiag override that messes up some times.
                % We can get around it by doing block diagonal manually.
                deltaChoi = expression(sprintf('deltaChoi(%1$d,%1$d)',dimAPrimeB));
                deltaChoi = directSumWorkArround(deltaChoi, deltaChoiBlocks);

                % unscramble from block diagonal form.
                deltaChoi = permMat'*deltaChoi*permMat;
                deltaChoi = (deltaChoi+deltaChoi')/2;
            else
                deltaChoi = variable(sprintf('deltaChoi(%1$d,%1$d)',dimAPrimeB),'hermitian');
            end


            %% objective
            minimize(gradPhi*deltaPhi + real(trace(gradChoiEveAPrimeB*deltaChoi)));

            %% constraints

            % Eve's channel is CPTP
            deltaChoi+choiEveAPrimeB == hermitian_semidefinite(dimAPrimeB);
            PartialTrace(deltaChoi+choiEveAPrimeB,2,[dimAPrime,dimB]) == eye(dimAPrime);

            %Accept point
            acceptPoint = testProb*[testCons.scalar];
            acceptPoint = [acceptPoint,1-testProb];

            % unique acceptance penalty
            phi+deltaPhi >= penaltyTermVariableLength(cvx_problem,testStateAAPrime,deltaChoi+choiEveAPrimeB,dimA,testCons,testProb,lambda);
            phi+deltaPhi >= 0;
            
            %% Get gradQES from dual variables
            acceptPoint - lambda == 0 : dualLambda;

            cvx_end

            %%
            %(-1 input as cvx_problem since we are outside of optimization routine)
            assert(penaltyTermVariableLength(-1,testStateAAPrime,deltaChoi+choiEveAPrimeB,dimA,testCons,testProb,lambda)>=0);
            % assert(gradPhi*deltaPhi + real(trace(gradChoiEveAPrimeB*deltaChoi)) <=0);

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
            gradQES = -full(dualLambda).';

            %store tolerance of solution
            tol = cvx_slvtol;
        end       
    
        function val = funcFW(vecFW,genStateAAPrime,keyProj,krausOps,testCons,testProb,gradQES,alpha,dimAPrime,perturbation,perturbationAff)
            arguments
                vecFW (:,1) double
                genStateAAPrime (:,:) double
                keyProj (:,1) cell
                krausOps (:,1) cell
                testCons (:,1) EqualityConstraint
                testProb (1,1) double
                gradQES (:,1) double
                alpha (1,1) double
                dimAPrime (1,1) uint64
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
                perturbationAff (1,1) double {mustBeInRange(perturbationAff,0,1)}
            end

            %Number of observations
            numObs = numel(testCons);

            %Reshape FW vector to Choi matrix
            dimAPrimeB = sqrt(numel(vecFW) - numObs);
            choiEveAPrimeB = reshape(vecFW(1:end-numObs),[dimAPrimeB,dimAPrimeB]);
            
            %Calculate rho conditioned on gen and test rounds
            dimA = size(genStateAAPrime,1)/dimAPrime;
            rhoABGen = PartialMap(genStateAAPrime,choiEveAPrimeB,2,[dimA,dimAPrime]);
            
            %Extract outputs in test rounds
            rhoCTestResults = vecFW(end-numObs+1:end).';

            sumTraceConTest = rhoCTestResults*2.^((alpha-1)*gradQES(1:end-1));
            
            %Value of QES entropy
            val = -1/(alpha-1)*log2((1-testProb)*2^(-(alpha-1)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,1,rhoABGen,keyProj,krausOps))*2^((alpha-1)*gradQES(end)) ...
                + testProb*sumTraceConTest);
        end

        function val = funcFW2(vecFW,genStateAAPrime,keyProj,krausOps,testCons,testProb,gradQES,alpha,dimAPrime,perturbation,perturbationAff)
            arguments
                vecFW (:,1) double
                genStateAAPrime (:,:) double
                keyProj (:,1) cell
                krausOps (:,1) cell
                testCons (:,1) EqualityConstraint
                testProb (1,1) double
                gradQES (:,1) double
                alpha (1,1) double
                dimAPrime (1,1) uint64
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
                perturbationAff (1,1) double {mustBeInRange(perturbationAff,0,1)}
            end

            %Number of observations
            numObs = numel(testCons);

            %Reshape FW vector to Choi matrix
            dimAPrimeB = sqrt(numel(vecFW) - numObs);
            choiEveAPrimeB = reshape(vecFW(1:end-numObs),[dimAPrimeB,dimAPrimeB]);
            
            %Calculate rho conditioned on gen and test rounds
            dimA = size(genStateAAPrime,1)/dimAPrime;
            rhoABGen = PartialMap(genStateAAPrime,choiEveAPrimeB,2,[dimA,dimAPrime]);
            
            %Extract outputs in test rounds
            rhoCTestResults = vecFW(end-numObs+1:end).';

            sumTraceConTest = rhoCTestResults*2.^((alpha-1)*gradQES(1:end-1));
            
            %Inner term of QES entropy
            val = -((1-testProb)*2^(-(alpha-1)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,1,rhoABGen,keyProj,krausOps))*2^((alpha-1)*gradQES(end)) ...
                + testProb*sumTraceConTest);
        end

        function grad = gradFuncFW(vecFW,genStateAAPrime,keyProj,krausOps,testCons,testProb,gradQES,alpha,dimAPrime,perturbation,perturbationAff)
            arguments
                vecFW (:,1) double
                genStateAAPrime (:,:) double
                keyProj (:,1) cell
                krausOps (:,1) cell
                testCons (:,1) EqualityConstraint
                testProb (1,1) double
                gradQES (:,1) double
                alpha (1,1) double
                dimAPrime (1,1) uint64
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
                perturbationAff (1,1) double {mustBeInRange(perturbationAff,0,1)}
            end

            %Number of observations
            numObs = numel(testCons);

            %Reshape FW vector to Choi matrix
            dimAPrimeB = sqrt(numel(vecFW) - numObs);
            choiEveAPrimeB = reshape(vecFW(1:end-numObs),[dimAPrimeB,dimAPrimeB]);
            
            %Calculate rho conditioned on gen and test rounds
            dimA = size(genStateAAPrime,1)/dimAPrime;
            rhoABGen = PartialMap(genStateAAPrime,choiEveAPrimeB,2,[dimA,dimAPrime]);
            dimB = size(rhoABGen,1)/dimA;
            
            %Extract outputs in test rounds
            rhoCTestResults = vecFW(end-numObs+1:end).';
            sumTraceConTest = rhoCTestResults*2.^((alpha-1)*gradQES(1:end-1));

            %Calculate inner term of log
            innerTerm = (1-testProb)*2^(-(alpha-1)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,1,rhoABGen,keyProj,krausOps))*2^((alpha-1)*gradQES(end)) ...
                + testProb*sumTraceConTest;
            
            %Gradient by RhoCTest
            gradByRhoCTest = -testProb/(alpha-1)/log(2)*1/innerTerm*2.^((alpha-1)*gradQES(1:end-1));

            %Gradient by Choi
            gradConEntByRho = RenyiTildeDown.gradConEnt(alpha,perturbation,perturbationAff,1,rhoABGen,keyProj,krausOps);
            
            % The choi matrix must be ordered input system tensor output system.
            genChoiAPrimeA = Swap(genStateAAPrime,[1,2],[dimA,dimAPrime]);
            tempDualMapGen = DualMap(genChoiAPrimeA,[dimAPrime,dimA]); %ordered [output,input] from taking dual
            
            gradConEntByChoi = PartialMap(gradConEntByRho{1},tempDualMapGen,1,[dimAPrime,dimB]);
            gradByChoi = (1-testProb)/innerTerm*2^(-(alpha-1)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,1,rhoABGen,keyProj,krausOps))...
                *2^((alpha-1)*gradQES(end))*gradConEntByChoi;           

            %Final gradient
            grad = [gradByChoi(:);gradByRhoCTest];
        end

        function grad = gradFuncFW2(vecFW,genStateAAPrime,keyProj,krausOps,testCons,testProb,gradQES,alpha,dimAPrime,perturbation,perturbationAff)
            arguments
                vecFW (:,1) double
                genStateAAPrime (:,:) double
                keyProj (:,1) cell
                krausOps (:,1) cell
                testCons (:,1) EqualityConstraint
                testProb (1,1) double
                gradQES (:,1) double
                alpha (1,1) double
                dimAPrime (1,1) uint64
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
                perturbationAff (1,1) double {mustBeInRange(perturbationAff,0,1)}
            end

            %Number of observations
            numObs = numel(testCons);

            %Reshape FW vector to Choi matrix
            dimAPrimeB = sqrt(numel(vecFW) - numObs);
            choiEveAPrimeB = reshape(vecFW(1:end-numObs),[dimAPrimeB,dimAPrimeB]);
            
            %Calculate rho conditioned on gen and test rounds
            dimA = size(genStateAAPrime,1)/dimAPrime;
            rhoABGen = PartialMap(genStateAAPrime,choiEveAPrimeB,2,[dimA,dimAPrime]);
            dimB = size(rhoABGen,1)/dimA;

            %Gradient by RhoCTest
            gradByRhoCTest = -testProb*2.^((alpha-1)*gradQES(1:end-1));

            %Gradient by Choi
            gradConEntByRho = RenyiTildeDown.gradConEnt(alpha,perturbation,perturbationAff,1,rhoABGen,keyProj,krausOps);
            
            % The choi matrix must be ordered input system tensor output system.
            genChoiAPrimeA = Swap(genStateAAPrime,[1,2],[dimA,dimAPrime]);
            tempDualMapGen = DualMap(genChoiAPrimeA,[dimAPrime,dimA]); %ordered [output,input] from taking dual

            gradConEntByChoi = PartialMap(gradConEntByRho{1},tempDualMapGen,1,[dimAPrime,dimB]);
            gradByChoi = (alpha-1)*log(2)*(1-testProb)*2^((alpha-1)*gradQES(end))*...
                        2^(-(alpha-1)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,1,rhoABGen,keyProj,krausOps))...
                        *gradConEntByChoi;           

            %Final gradient
            grad = [gradByChoi(:);gradByRhoCTest];
        end


        function vecQesInit = initPoint(vecFW,testStateAAPrime,testCons,dimAPrime)
            arguments
                vecFW (:,1) double
                testStateAAPrime (:,:) double                
                testCons (:,1) EqualityConstraint
                dimAPrime (1,1) uint64
            end

            %Reshape FW vector to Choi matrix
            dimAPrimeB = sqrt(numel(vecFW));
            choiEveAPrimeB = reshape(vecFW(1:end),[dimAPrimeB,dimAPrimeB]);
            
            %Calculate rho conditioned on test rounds
            dimA = size(testStateAAPrime,1)/dimAPrime;
            rhoABTest = PartialMap(testStateAAPrime,choiEveAPrimeB,2,[dimA,dimAPrime]);

            for index = numel(testCons):-1:1
                rhoCTestResults(index) = real(trace(testCons(index).operator*rhoABTest));
            end
            
            %Stack initial point as vector together
            vecQesInit = [vecFW;rhoCTestResults.'];
        end

        function [deltaVecFW,exitFlag,statusMessage,tol] = subproblemQes(vecFW,vecFWGrad,...
                testStateAAPrime,dimAPrime,testCons,options)
            arguments
                vecFW (:,1) double
                vecFWGrad (:,1) double
                testStateAAPrime (:,:) double
                dimAPrime (1,1) uint64
                testCons (:,1) EqualityConstraint
                options (1,1) struct
            end

            %Linear constraint tolerance
            linConTol = 1/2*options.linearConstraintTolerance;

            %Number of observations
            numObs = numel(testCons);

            %Reshape FW vector and gradient to Choi matrix and gradient
            dimAPrimeB = sqrt(numel(vecFW) - numObs);
            choiEveAPrimeB = reshape(vecFW(1:end-numObs),[dimAPrimeB,dimAPrimeB]);
            gradChoiEveAPrimeB = reshape(vecFWGrad(1:end-numObs),[dimAPrimeB,dimAPrimeB]);          
            

            %Extract outputs in test rounds and gradient from FW vector
            rhoCTestResults = vecFW(end-numObs+1:end);
            gradRhoC = vecFWGrad(end-numObs+1:end);

            %Additional dimensions required
            dimB = dimAPrimeB/dimAPrime;
            dimAAPrime = size(testStateAAPrime,1);
            dimA = dimAAPrime/dimAPrime;

            %% CVX
            cvx_begin sdp
            cvx_solver(convertStringsToChars(options.cvxSolver));
            cvx_precision(convertStringsToChars(options.cvxPrecision));
            cvx_quiet(options.verboseLevel<2);

            deltaRhoC = variable('deltaRhoC(numObs,1)');

            if options.blockDiagonal

                %% block diagonal set up which cant be placed in a separate function
                permMat = options.blockP;
                newDims = options.newDims.';

                deltaChoiStrings = compose("deltaChoi%1d(%d,%d)",(1:numel(newDims)).',newDims,newDims);
                % CVX MUST generates variables in this work space. We have to do each
                % separately.
                for index = 1:numel(deltaChoiStrings)
                    variable(convertStringsToChars(deltaChoiStrings(index)),'hermitian')
                end
                % Get a single cell array to capture those variables.
                deltaChoiBlocksString = "{"+strjoin(compose("deltaChoi%d",1:numel(newDims)),",")+"}";
                deltaChoiBlocks = eval(deltaChoiBlocksString);


                % there is a bug in cvx's blkdiag override that messes up some times.
                % We can get around it by doing block diagonal manually.
                deltaChoi = expression(sprintf('deltaChoi(%1$d,%1$d)',dimAPrimeB));
                deltaChoi = directSumWorkArround(deltaChoi, deltaChoiBlocks);

                % unscramble from block diagonal form.
                deltaChoi = permMat'*deltaChoi*permMat;
                deltaChoi = (deltaChoi+deltaChoi')/2;
            else
                deltaChoi = variable(sprintf('deltaChoi(%1$d,%1$d)',dimAPrimeB),'hermitian');
            end

            %% objective
            minimize(gradRhoC.'*deltaRhoC + real(trace(gradChoiEveAPrimeB*deltaChoi)));

            %% constraints

            % Eve's channel is CPTP
            deltaChoi+choiEveAPrimeB == hermitian_semidefinite(dimAPrimeB);
            PartialTrace(deltaChoi+choiEveAPrimeB,2,[dimAPrime,dimB]) == eye(dimAPrime);

            % RhoCTest(k) = Tr(Gamma_k rho_AB|Test)
            %Find rhoAB conditioned on test rounds 
            rhoABTest = PartialMap(testStateAAPrime,deltaChoi+choiEveAPrimeB,2,[dimA,dimAPrime]);

            for index = 1:numel(testCons)
                abs(deltaRhoC(index) + rhoCTestResults(index) - real(trace(testCons(index).operator*rhoABTest))) <= linConTol;
            end
            
            % rho_C + deltaRho_C in [0,1]
            deltaRhoC + rhoCTestResults == nonnegative(numObs);
            deltaRhoC + rhoCTestResults <= ones(numObs,1);

            cvx_end

            % reconstruct rho from blocks incase it didn't get filled in by CVX
            if options.blockDiagonal
                deltaChoi = zeros(dimAPrimeB);
                deltaChoiBlocks = eval(deltaChoiBlocksString);
                deltaChoi = directSumWorkArround(deltaChoi, deltaChoiBlocks);
                deltaChoi = permMat'*deltaChoi*permMat;
            end
            deltaChoi = (deltaChoi + deltaChoi')/2;
            
            %Return final delta vector
            deltaVecFW = [deltaChoi(:);deltaRhoC];
            
            %Store status and exit flag
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
    end
end

function penaltyVal = penaltyTermVariableLength(cvx_problem,testStateAAPrime,ChoiEveAPrimeB,dimA,testCons,testProb,lambda)
    arguments
        cvx_problem
        testStateAAPrime (:,:) double
        ChoiEveAPrimeB (:,:) % double or CVX
        dimA (1,1) double
        testCons (:,1) EqualityConstraint                
        testProb (1,1) double
        lambda (:,:) % double or CVX
    end
    
    dimAPrime = size(testStateAAPrime,1)/dimA;
    
    % rhoABTest is conditioned on test rounds!
    rhoABTest = PartialMap(testStateAAPrime,ChoiEveAPrimeB,2,[dimA,dimAPrime]);
    
    for index = numel(testCons):-1:1
        rhoCTestResults(index) = real(trace(testCons(index).operator*rhoABTest));
    end            
    
    %D(lambda||rhoC)
    penaltyVal = sum(rel_entr(lambda(1:end-1),testProb*rhoCTestResults))/log(2)...
                + rel_entr(lambda(end),(1-testProb))/log(2);
end