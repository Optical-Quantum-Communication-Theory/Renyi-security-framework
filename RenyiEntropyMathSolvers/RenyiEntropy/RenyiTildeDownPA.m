classdef RenyiTildeDownPA
    % Collection of the functions needed to set up the basic Renyi sandwich
    % (tilde) down objective function, including the classical relative
    % entropy soft penalty. See RenyiTildeDown for more details on the
    % Renyi term. The objective function is constructed as an optimization
    % over Eve's channel from A' to B. For more details see
    % https://arxiv.org/abs/1710.05511v2.
    %
    % This simple version does not support any additional block structure.
    %
    % The Frank Wolfe vectorization of the inputs is
    %
    % vecFW = [choiEveAPrimeB(:); phi]
    %
    % where choiEveAPrimeB is the Eve's choi matrix (systems ordered A'B),
    % and phi is the slack variable for the classical relative entropy cone
    % soft penalty. The gradient similarly follows
    %
    % gradVecFW = [gradChoiEveAPrimeB(:); gradByPhi]
    %
    %
    % see also RenyiTildeDown, FrankWolfe

    methods (Static)

        %% Frank Wolfe functions
        function val = funcFW(vecFW,genStateAAPrime,keyProj,krausOps,...
                testProb,alpha,dimAPrime,perturbation,perturbationAff)
            % Frank Wolfe objective function
            %
            % NOTE that inputs have very few checks on them for
            % performance. Take care when using this function.
            %
            % Inputs:
            % * vecFW: Vectorized channel Choi matrix and the classical
            %   relative entropy cone soft penalty.
            % * genStateAAPrime: Alice's source replaced state conditioned
            %   on generation rounds. A density matrix (and likely pure).
            % * keyProj: A cell array of the Projection operators of the
            %   pinching channel, Z, acting on the key register used in the
            %   solver.
            % * krausOps: A cell array of the Kraus operators for the
            %   post-selection map, G, of Alice and Bob. Kraus operators
            %   that form a CPTNI map are effectively completed to a CPTP
            %   map where all remaining weight is mapped to a discard
            %   state.
            % * testProb: Probability of selecting a test round, (p(test) =
            %   1 - p(gen)).
            % * alpha: The Renyi entropy parameter. Alpha must be in the
            %   half open interval (1,2].
            % * dimAPrime: Dimension of the signal sent from Alice into the
            %   channel.
            % * perturbation: Linear perturbation for G(rho). In the
            %   range [0,1].
            % * perturbationAff: Affine perturbation for G(rho). In the
            %   range [0,1].
            %
            % Output:
            % * val: Value of the objective function.
            %
            % see also RenyiTildeDownPA, RenyiTildeDown, FrankWolfe
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


        function gradVecFW = gradFuncFW(vecFW,genStateAAPrime,keyProj,...
                krausOps,testProb,alpha,dimAPrime,perturbation,perturbationAff)
            % Gradient of the Frank Wolfe objective function
            %
            % NOTE that inputs have very few checks on them for
            % performance. Take care when using this function.
            %
            % Inputs:
            % * vecFW: Vectorized channel Choi matrix and the classical
            %   relative entropy cone soft penalty. Gradient is evaluated
            %   at this point.
            % * genStateAAPrime: Alice's source replaced state conditioned
            %   on generation rounds. A density matrix (and likely pure).
            % * keyProj: A cell array of the Projection operators of the
            %   pinching channel, Z, acting on the key register used in the
            %   solver.
            % * krausOps: A cell array of the Kraus operators for the
            %   post-selection map, G, of Alice and Bob. Kraus operators
            %   that form a CPTNI map are effectively completed to a CPTP
            %   map where all remaining weight is mapped to a discard
            %   state.
            % * testProb: Probability of selecting a test round, (p(test) =
            %   1 - p(gen)).
            % * alpha: The Renyi entropy parameter. Alpha must be in the
            %   half open interval (1,2].
            % * dimAPrime: Dimension of the signal sent from Alice into the
            %   channel.
            % * perturbation: Linear perturbation for G(rho). In the
            %   range [0,1].
            % * perturbationAff: Affine perturbation for G(rho). In the
            %   range [0,1].
            %
            % Output:
            % * gradVecFW: Gradient of the objective function evaluated at
            %   vecFW. Has the same formatting as vecFW.
            %
            % see also RenyiTildeDownPA, RenyiTildeDown, FrankWolfe
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

            [gradChoiEveAPrimeB,gradByPhi] = gradPrivacyAmpSplit(phi, ...
                choiEveAPrimeB,genStateAAPrime,keyProj,krausOps,testProb,...
                alpha,dimAPrime,perturbation,perturbationAff);

            gradVecFW = [gradChoiEveAPrimeB(:); gradByPhi];
        end


        function [deltaVecFW,exitFlag,statusMessage] = subproblemUnique(vecFW,vecFWGrad,...
                testStateAAPrime,dimAPrime,testProb,testCons,options)
            % Direction finding convex subproblem for the Frank Wolfe
            % optimization. This one is for a unique acceptance frequency
            % distribution.
            %
            % See FrankWolfe for more details.
            %
            % NOTE that inputs have very few checks on them for
            % performance. Take care when using this function.
            %
            % Inputs:
            % * vecFW: Vectorized channel Choi matrix and the classical
            %   relative entropy cone soft penalty.
            % * gradVecFW: Gradient of the objective function evaluated at
            %   vecFW.
            % * testStateAAPrime: Alice's source replaced state conditioned
            %   on test rounds. A density matrix (and likely pure).
            % * dimAPrime: Dimension of the signal sent from Alice into the
            %   channel.
            % * testProb: Probability of selecting a test round, (p(test) =
            %   1 - p(gen)).
            % * testCons: An array of class EqualityConstraint.
            %   Represents the constraints from test rounds. For each con
            %   in testCons, con.scalar represents the accepted
            %   frequency conditioned on test rounds, and con.operator is
            %   the corresponding operator of Alice and Bob's joint POVM
            %   element conditioned on test rounds. The sum of all
            %   operators should be identity and the sum of all scalars
            %   should be 1.
            % * options: Options struct taken from the math solver. Same
            %   one we always use.
            %
            % Outputs:
            % * deltaVecFW: Optimal step direction (and size) from the
            %   convex subproblem.
            % * exitFlag: SubProblemExitFlag for the convex problem.
            % * statusMessage: CVX status message string from the
            %   optimization.
            %
            % see also RenyiTildeDownPA, FrankWolfe, SubProblemExitFlag
            arguments
                vecFW (:,1) double
                vecFWGrad (:,1) double
                testStateAAPrime (:,:) double
                dimAPrime (1,1) uint64
                testProb (1,1) double {mustBeInRange(testProb, 0, 1)}
                testCons (:,1) EqualityConstraint
                options (1,1) struct

                % options:
                % options.cvxSolver
                % options.cvxPrecision
                % options.verboseLevel
                % options.blockDiagonal
                % options.blockP
                % options.newDims
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

            % sanity check
            assert(penaltyTermUnique(testStateAAPrime,deltaChoi+choiEveAPrimeB,dimA,testCons,testProb)>=0);
            % assert(gradPhi*deltaPhi + real(trace(gradChoiEveAPrimeB*deltaChoi)) <=0);

            % reconstruct rho from blocks encase it didn't get filled in by CVX
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
        function [vecFW,exitFlag] = closestVecFW(choi0,testStateAAPrime,testProb,testCons,dimA,dimAPrime,dimB,options)
            % Solving for closest Choi matrix compatible with our
            % observations. Used to Find an initial point for the Frank
            % Wolfe routine. The initial point is selected by finding a
            % Choi matrix that minimizes the classical relative entropy
            % soft penalty term (which gives phi).
            %
            % The argument choi0 is from an earlier initial point routine.
            % It is now ignored, but you can still find the commented out
            % code in here if you want to check it out.
            %
            % Inputs:
            % * choi0: This term is now ignored. Initial Choi matrix to
            %   minimize distance to.
            % * testStateAAPrime: Alice's source replaced state conditioned
            %   on test rounds. A density matrix (and likely pure).
            % * testProb: Probability of selecting a test round, (p(test) =
            %   1 - p(gen)).
            % * testCons: An array of class EqualityConstraint.
            %   Represents the constraints from test rounds. For each con
            %   in testCons, con.scalar represents the accepted
            %   frequency conditioned on test rounds, and con.operator is
            %   the corresponding operator of Alice and Bob's joint POVM
            %   element conditioned on test rounds. The sum of all
            %   operators should be identity and the sum of all scalars
            %   should be 1.
            % * dimA: Dimension of Alice's system she measures.
            % * dimAPrime: Dimension of the signal sent from Alice into the
            %   channel.
            % * dimB: Dimension of Bob's system he measures.
            % * options: Options struct taken from the math solver. Same
            %   one we always use.
            %
            % Outputs:
            % * vecFW: An initial point for the Frank Wolfe routine.
            %   Vectorized channel Choi matrix and the classical relative
            %   entropy cone soft penalty.
            % * exitFlag: SubProblemExitFlag for the convex problem.
            %
            % see also RenyiTidleDownPA, FrankWolfe, SubProblemExitFlag

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


        function [lowerBound,exitFlag] = simpleStep2(constOffsets,gradVecsFW,...
                testStateAAPrime,dimAPrime,testProb,testCons,options)
            % Experimental step 2 solver for achieving a lower bound from
            % Frank Wolfe. Instead of running the subproblem one more time
            % to get a lower bound, it instead minimizes the maximum value
            % of multiple linearizations. Sometimes this worked way better.
            % Sometimes it just wasted time.
            %
            % Inputs:
            % * constOffsets: vector of the constant offsets "up" from the
            %   origin for each affine hyperplane.
            % * gradVecFWs: Cell array of gradients of the objective
            %   function used for the affine hyperplanes.
            % * dimAPrime: Dimension of the signal sent from Alice into the
            %   channel.
            % * testProb: Probability of selecting a test round, (p(test) =
            %   1 - p(gen)).
             % * testCons: An array of class EqualityConstraint.
            %   Represents the constraints from test rounds. For each con
            %   in testCons, con.scalar represents the accepted
            %   frequency conditioned on test rounds, and con.operator is
            %   the corresponding operator of Alice and Bob's joint POVM
            %   element conditioned on test rounds. The sum of all
            %   operators should be identity and the sum of all scalars
            %   should be 1.
            % * options: Options struct taken from the math solver. Same
            %   one we always use.
            %
            % Outputs:
            % * lowerBound: Lower bound on the Frank Wolfe objective
            %   function found by CVX.
            % * exitFlag: SubProblemExitFlag for the convex problem.
            %
            % see also RenyiTidleDownPA, FrankWolfe, SubProblemExitFlag
            arguments
                constOffsets (:,1) double
                gradVecsFW (:,1) cell % of (:,1) doubles
                testStateAAPrime (:,:) double
                dimAPrime (1,1) uint64
                testProb (1,1) double {mustBeInRange(testProb, 0, 1)}
                testCons (:,1) EqualityConstraint
                options (1,1) struct
            end

            dimAPrimeB = sqrt(numel(gradVecsFW{1})-1);

            gradPhi = cellfun(@(x) x(end),gradVecsFW);
            gradChoisEveAPrimeB = cellfun(@(x) reshape(x(1:end-1),[dimAPrimeB,dimAPrimeB]),...
                gradVecsFW,"UniformOutput",false);

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
% privacy amplification function with vecFW decomposed back into phi and
% choiEveAPrimeB.
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
% gradient of the privacy amplification function with vecFW decomposed back
% into phi and choiEveAPrimeB.
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
% classical relative entropy soft penalty term specialized for unique
% acceptance.
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