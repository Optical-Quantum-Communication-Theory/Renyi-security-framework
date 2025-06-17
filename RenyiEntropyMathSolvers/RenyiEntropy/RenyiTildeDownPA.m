classdef RenyiTildeDownPA

    methods (Static)
        
        %% Frank Wolfe functions
        function val = funcFW(vecFW,genStateAAPrime,keyProj,krausOps,testProb,alpha,dimAPrime,perturbation,perturbationAff)
            arguments
                vecFW (:,1) double
                genStateAAPrime (:,:) double
                keyProj (:,1) cell
                krausOps (:,1) cell
                testProb (1,1) double
                alpha (1,1) double
                dimAPrime (1,1) uint64
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
                perturbationAff (1,1) double {mustBeInRange(perturbationAff,0,1)}
            end
            phi = vecFW(end);

            dimAPrimeB = sqrt(numel(vecFW)-1);
            choiEveAPrimeB = reshape(vecFW(1:end-1),[dimAPrimeB,dimAPrimeB]);

            val = privacyAmpSplit(phi,choiEveAPrimeB,genStateAAPrime, ...
                keyProj,krausOps,testProb,alpha,dimAPrime,perturbation,perturbationAff);
        end


        function gradValVec = gradFuncFW(vecFW,genStateAAPrime,keyProj,krausOps,testProb,alpha,dimAPrime,perturbation,perturbationAff)
            arguments
                vecFW (:,1) double
                genStateAAPrime (:,:) double
                keyProj (:,1) cell
                krausOps (:,1) cell
                testProb (1,1) double
                alpha (1,1) double
                dimAPrime (1,1) uint64
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
                perturbationAff (1,1) double {mustBeInRange(perturbationAff,0,1)}
            end
            dimAPrimeB = sqrt(numel(vecFW)-1);
            choiEveAPrimeB = reshape(vecFW(1:end-1),[dimAPrimeB,dimAPrimeB]);
            phi = vecFW(end);

            [gradByChoi,gradByPhi] = gradPrivacyAmpSplit(phi, ...
                choiEveAPrimeB,genStateAAPrime,keyProj,krausOps,testProb,...
                alpha,dimAPrime,perturbation,perturbationAff);

            gradValVec = [gradByChoi(:); gradByPhi];
        end


        function [deltaVecFW,exitFlag,statusMessage] = subproblemUnique(vecFW,vecFWGrad,...
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
            
            %Phi and gradient
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

            if options.blockDiagonal

                %% block diagonal set up which cant be placed in a separate function
                permMat = options.blockP;
                newDims = options.newDims.';

                deltaChoiStrings = compose("deltaChoi%1d(%d,%d)",(1:numel(newDims)).',newDims,newDims);
                % CVX MUST generate variables in this work space. We have to do each
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

            % unique acceptance penalty
            phi+deltaPhi >= penaltyTermUnique(testStateAAPrime,deltaChoi+choiEveAPrimeB,dimA,testCons,testProb);
            phi+deltaPhi >= 0;
            cvx_end


            assert(penaltyTermUnique(testStateAAPrime,deltaChoi+choiEveAPrimeB,dimA,testCons,testProb)>=0);
            % assert(gradPhi*deltaPhi + real(trace(gradChoiEveAPrimeB*deltaChoi)) <=0);

            % reconstruct rho from blocks incase it didn't get filled in by CVX
            if options.blockDiagonal
                deltaChoi = zeros(dimAPrimeB);
                deltaChoiBlocks = eval(deltaChoiBlocksString);
                deltaChoi = directSumWorkArround(deltaChoi, deltaChoiBlocks);
                deltaChoi = permMat'*deltaChoi*permMat;
            end

            deltaChoi = (deltaChoi + deltaChoi')/2;

            deltaVecFW = [deltaChoi(:);deltaPhi];

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

        %% Initial point
        function vecFW = closestVecFW(choi0,testStateAAPrime,testProb,testCons,dimA,dimAPrime,dimB,options)
            %solving for closest Choi matrix compatible with our observations


            % the penalty at the starting choi matrix.
            % phi0 = penaltyTermUnique(testStateAAPrime,choi0,dimA,testCons,testProb);

            cvx_begin sdp
            cvx_solver(convertStringsToChars(options.cvxSolver));
            cvx_precision(convertStringsToChars(options.cvxPrecision));
            cvx_quiet(options.verboseLevel<2);

            %Dimension of Choi matrix
            dimAPrimeB = dimAPrime*dimB;

            phi = variable('phi(1,1)', 'nonnegative');

            if options.blockDiagonal

                %% block diagonal set up which cant be placed in a separate function
                permMat = options.blockP;
                newDims = options.newDims.';

                choiEveAPrimeBStrings = compose("choiEveAPrimeB%1d(%d,%d)",(1:numel(newDims)).',newDims,newDims);
                % CVX MUST generates variables in this work space. We have to do each
                % separately.
                for index = 1:numel(choiEveAPrimeBStrings)
                    variable(convertStringsToChars(choiEveAPrimeBStrings(index)),'hermitian','semidefinite')
                end
                % Get a single cell array to capture those variables.
                choiEveAPrimeBBlocksString = "{"+strjoin(compose("choiEveAPrimeB%d",1:numel(newDims)),",")+"}";
                choiEveAPrimeBBlocks = eval(choiEveAPrimeBBlocksString);


                % there is a bug in cvx's blkdiag override that messes up some times.
                % We can get around it by doing block diagonal manually.
                choiEveAPrimeB = expression(sprintf('choiEveAPrimeB(%1$d,%1$d)',dimAPrimeB));
                choiEveAPrimeB = directSumWorkArround(choiEveAPrimeB, choiEveAPrimeBBlocks);

                % unscramble from block diagonal form.
                choiEveAPrimeB = permMat'*choiEveAPrimeB*permMat;
                choiEveAPrimeB = (choiEveAPrimeB+choiEveAPrimeB')/2;
            else
                choiEveAPrimeB = variable(sprintf('choiEveAPrimeB(%1$d,%1$d)',dimAPrimeB),'hermitian','semidefinite');
            end

            %% objective
            % minimize( norm(choi0 - choiEveAPrimeB,"fro") + norm(phi0 - phi))
            minimize(phi)

            %% Constraints

            % Eve's channel is TP
            PartialTrace(choiEveAPrimeB,2,[dimAPrime,dimB]) == eye(dimAPrime);

            % unique acceptance test penalty term
            phi >= penaltyTermUnique(testStateAAPrime,choiEveAPrimeB,dimA,testCons,testProb);

            % % rhoABTest is conditioned on test rounds!
            % rhoABTest = PartialMap(testStateAAPrime,choiEveAPrimeB,2,[dimA,dimAPrime]);
            % 
            % for index = 1:numel(testCons)
            %     testCons(index).scalar == real(trace(testCons(index).operator*rhoABTest));
            % end

            cvx_end

            if options.blockDiagonal
                choiEveAPrimeB = zeros(dimAPrimeB);
                choiEveAPrimeBBlocks = eval(choiEveAPrimeBBlocksString);
                choiEveAPrimeB = directSumWorkArround(choiEveAPrimeB, choiEveAPrimeBBlocks);
                choiEveAPrimeB = permMat'*choiEveAPrimeB*permMat;
            end

            %Make sure resulting Choi state is hermitian
            choiEveAPrimeB = full(choiEveAPrimeB + choiEveAPrimeB')/2;

            % phi = penaltyTermUnique(testStateAAPrime,choiEveAPrimeB,dimA,testCons,testProb);

            % return vecFW
            vecFW = [choiEveAPrimeB(:); phi];
        end


        function [lowerBound,exitFlag] = simpleStep2(constOffsets,vecsFWGrad,...
                testStateAAPrime,dimAPrime,testProb,testCons,options)
            arguments
                constOffsets (:,1) double
                vecsFWGrad (:,1) cell % of (:,1) doubles
                testStateAAPrime (:,:) double
                dimAPrime (1,1) uint64
                testProb (1,1) double {mustBeInRange(testProb, 0, 1)}
                testCons (:,1) EqualityConstraint
                options (1,1) struct
            end

            dimAPrimeB = sqrt(numel(vecsFWGrad{1})-1);

            gradPhi = cellfun(@(x) x(end),vecsFWGrad);
            gradChoisEveAPrimeB = cellfun(@(x) reshape(x(1:end-1),[dimAPrimeB,dimAPrimeB]),...
                vecsFWGrad,"UniformOutput",false);

            dimB = dimAPrimeB/dimAPrime;
            dimAAPrime = size(testStateAAPrime,1);
            dimA = dimAAPrime/dimAPrime;


            %% CVX
            cvx_begin sdp
            cvx_solver(convertStringsToChars(options.cvxSolver));
            cvx_precision(convertStringsToChars(options.cvxPrecision));
            cvx_quiet(options.verboseLevel<2);


            phi = variable('phi(1,1)','nonnegative');

            if options.blockDiagonal

                %% block diagonal set up which cant be placed in a separate function
                permMat = options.blockP;
                newDims = options.newDims.';

                choiStrings = compose("choi%1d(%d,%d)",(1:numel(newDims)).',newDims,newDims);
                % CVX MUST generates variables in this work space. We have to do each
                % separately.
                for index = 1:numel(choiStrings)
                    variable(convertStringsToChars(choiStrings(index)),'hermitian','semidefinite')
                end
                % Get a single cell array to capture those variables.
                choiBlocksString = "{"+strjoin(compose("choi%d",1:numel(newDims)),",")+"}";
                choiBlocks = eval(choiBlocksString);


                % there is a bug in cvx's blkdiag override that messes up some times.
                % We can get around it by doing block diagonal manually.
                choi = expression(sprintf('choi(%1$d,%1$d)',dimAPrimeB));
                choi = directSumWorkArround(choi, choiBlocks);

                % unscramble from block diagonal form.
                choi = permMat'*choi*permMat;
                choi = (choi+choi')/2;
            else
                choi = variable(sprintf('choi(%1$d,%1$d)',dimAPrimeB),'hermitian','semidefinite');
            end


            %% objective
            temp = cellfun(@(x) real(trace(x*choi)),gradChoisEveAPrimeB,"UniformOutput",false);
            objective = constOffsets + gradPhi.*phi + [temp{:}].';
            minimize(max(objective));

            %% constraints

            % Eve's channel is CPTP
            PartialTrace(choi,2,[dimAPrime,dimB]) == eye(dimAPrime);

            % unique acceptance penalty
            phi >= penaltyTermUnique(testStateAAPrime,choi,dimA,testCons,testProb);
            cvx_end

            lowerBound = cvx_optval;

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


function cost = privacyAmpSplit(phi,choiEveAPrimeB,genStateAAPrime,keyProj,krausOps,testProb,alpha,dimAPrime,perturbation,perturbationAff)
arguments
    % minimal checks for shape and type
    phi (1,1) double
    choiEveAPrimeB (:,:) double
    genStateAAPrime (:,:) double
    keyProj (:,1) cell
    krausOps (:,1) cell
    testProb (1,1) double
    alpha (1,1) double
    dimAPrime (1,1) uint64
    perturbation (1,1) double
    perturbationAff (1,1) double
end

dimA = size(genStateAAPrime,1)/dimAPrime;

rhoAB = PartialMap(genStateAAPrime,choiEveAPrimeB,2,[dimA,dimAPrime]);

%version with alphaHat
% cost = phi/(alpha-1) + (1-testProb)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,1,rhoAB,keyProj,krausOps);

%version with alpha
cost = alpha/(alpha-1)*phi + (1-testProb)*RenyiTildeDown.conEnt(alpha,perturbation,perturbationAff,1,rhoAB,keyProj,krausOps);
end

function [gradByChoi,gradByPhi] = gradPrivacyAmpSplit(~,choiEveAPrimeB,genStateAAPrime,keyProj,krausOps,testProb,alpha,dimAPrime,perturbation,perturbationAff)
    arguments
        % minimal checks for shape and type
        ~ %phi (1,1) double Not needed
        choiEveAPrimeB (:,:) double
        genStateAAPrime (:,:) double
        keyProj (:,1) cell
        krausOps (:,1) cell
        testProb (1,1) double
        alpha (1,1) double
        dimAPrime (1,1) uint64
        perturbation (1,1) double
        perturbationAff (1,1) double
    end
    %gradient of phi term
    %version with alphaHat
    % gradByPhi = 1/(alpha-1);

    %version with alpha
    gradByPhi = alpha/(alpha-1);
    
    %gradient by Choi
    dimA = size(genStateAAPrime,1)/dimAPrime;
    
    rhoAB = PartialMap(genStateAAPrime,choiEveAPrimeB,2,[dimA,dimAPrime]);
    dimB = size(rhoAB,1)/dimA;
    
    gradByChoi = RenyiTildeDown.gradConEnt(alpha,perturbation,perturbationAff,1,rhoAB,keyProj,krausOps);
    
    % The choi matrix must be ordered input system tensor output system,
    % hence we swap
    genChoiAPrimeA = Swap(genStateAAPrime,[1,2],[dimA,dimAPrime]);
    
    dualGenChoi = DualMap(genChoiAPrimeA,[dimAPrime,dimA]); %ordered [output,input] from taking dual
    gradByChoi = (1-testProb)*PartialMap(gradByChoi{1},dualGenChoi,1,[dimA,dimB]);

end

function penaltyVal = penaltyTermUnique(testStateAAPrime,ChoiEveAPrimeB,dimA,testCons,testProb)
arguments
    testStateAAPrime (:,:) double
    ChoiEveAPrimeB (:,:) % double or CVX
    dimA (1,1) double
    testCons (:,1) EqualityConstraint
    testProb (1,1) double
end

dimAPrime = size(testStateAAPrime,1)/dimA;

% rhoABTest is conditioned on test rounds!
rhoABTest = PartialMap(testStateAAPrime,ChoiEveAPrimeB,2,[dimA,dimAPrime]);

for index = numel(testCons):-1:1
    rhoCTestResults(index) = real(trace(testCons(index).operator*rhoABTest));
end

acceptPoint = testProb*[testCons.scalar];

penaltyVal = sum(rel_entr(acceptPoint,testProb*rhoCTestResults))/log(2)...
    ...% this part should always be zero for unique acceptance
    + (1-testProb)*log2((1-testProb)/(1-testProb));
end