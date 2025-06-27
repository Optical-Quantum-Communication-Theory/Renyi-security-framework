classdef RenyiTildeDown
    % Class containing functions to evaluate the sandwiched conditional
    % Renyi entropy (tilde down arrow in https://arxiv.org/abs/1504.00233)
    % of QKD protocols by applying CPTNI Kraus operators and a projection
    % map. Specifically, the CPTNI map and projection map are the G and Z
    % maps from https://arxiv.org/abs/1710.05511v2.
    %
    % The conditional entropy restricted to QKD protocols realized via the
    % G and Z maps is given by:
    %
    % \bar{g}_\alpha(\rho_{AB}) :=
    % \tilde{H}^\downarrow_\alpha(S|EC)_{\mathcal{M}_{\text{QKD}_{AB
    % \rightarrow SQC}} (\ketbra{\rho}_{ABE}) } = \frac{-1}{\alpha-1} \log
    % \circ \theta \circ \mathcal{G}(\rho_{AB}),
    %
    % with
    %
    % \theta(\sigma) = \text{Tr} \circ \text{Pow}_\alpha \circ \mathcal{Z}
    % \circ \text{Pow}_{1/\alpha}(\sigma) + 1 - \text{Tr}[\sigma],
    %
    % where S is the secret key register, E is Eve's quantum register, C is
    % the classical announcements and AB is Alice and Bob's registers. Q is
    % what ever remaining registers are needed by the G map.
    %
    % All functions support perturbations of the G map to improve stability
    % while keeping a lower bound. The G map is replaced with
    %
    % \Psi_{\epsilon_2} \circ \Delta_{\epsilon_1} \circ \mathcal{G},
    %
    % where \Delta_{\epsilon_1} is the perturbation channel
    %
    % \Delta_{\epsilon_1}(\rho_F) = (1-\epsilon_1)\rho_F + \epsilon_1
    % \frac{\Tr[\rho_F] I_F}{\dim(F)},
    %
    % and \Psi_{\epsilon_2} is the affine perturbation map
    %
    % \Psi_{\epsilon_2}(\rho_F) = (1-\epsilon_2)\rho_F + \epsilon_2
    % \frac{I_F}{\dim(F)}.
    %
    % These perturbations are then used to compute a lower bound on the
    % eigenvalues for the matrix powers used in the conditional entropy and
    % gradients.
    %
    % Finally, `RenyiTildeDown.conEnt` and `RenyiTildeDown.gradConEnt` support
    % weighted block diagonal structure. Let \{(p_i, \rho_{|i},
    % \mathcal{G}_{|i}, \mathcal{Z}_{|i})\}_{i=1}^n be a set of n tuples
    % where p_i \in [0,1] and \sum_i p_i \leq 1, and \rho_{|i},
    % \mathcal{G}_{|i}, and \mathcal{Z}_{|i} are defined as above. For the
    % subnormalized block diagonal state \bigoplus_{i=1}^n p_i \rho_i, and
    % block diagonal G and Z maps, the conditional Renyi sandwich entropy
    % is given by
    % 
    % \bar{g}_\alpha(\rho) = \frac{-1}{\alpha-1} \log \left(\sum_{i=1}^n
    % \theta_{|i} \circ \mathcal{G_{|i}}(\rho_{|i}) + 1 - \sum_{i=1}^n p_i
    % \right),
    %
    % where \theta{|i} uses \mathcal{Z}_{|i}. Here if \sum_i p_i < 1 then
    % the remaining weight is assumed to be discarded (hence the 1-\sum_i
    % p_i term).
    %
    % An important observation is that the perturbation and the block
    % structure can be combined together. First, we rewrite the weighted
    % block diagonal structure using an affine function \mathcal{B}
    %
    % \bar{g}_\alpha(\rho) = \bar{g}_\alpha \circ \mathcal{B}(\rho_{|1},
    % \dots, \rho_{|n}),
    %
    % where
    %
    % \mathcal{B}(\rho_{F_1|1}, \dots, \rho_{F_n|n}) = \bigoplus_{i=1}^n
    % \rho_{F_i|i} \oplus (1-\sum_{i=1}^n p_i) \ketbra{\bot}.
    %
    % Then (using the same values of \epsilon_1 and \epsilon_2) we can
    % apply the perturbation maps to each sub space F_i individually. In
    % other words,
    %
    % \Delta_{\epsilon_1|i}(\rho_{F_i|i}) = (1-\epsilon_1)\rho_{F_i|i} +
    % \epsilon_1 \frac{\Tr[\rho_{F_i|i}] I_{F_i}}{\dim(F_i)}.
    %
    % \Psi_{\epsilon_2|i} is defined analogously. The final result
    % still recovers the same perturbation lower bound
    %
    % \bar{g}_\alpha \circ \mathcal{B} \left( \Delta_{\epsilon_1|1} \circ
    % \Psi_{\epsilon_2|1} \circ G_{|1}(\rho_{F_1|1}), \dots
    % \Delta_{\epsilon_1|n} \circ \Psi_{\epsilon_2|n} \circ
    % G_{|n}(\rho_{F_n|n}) \right) \leq (1-\epsilon_1)(1-\epsilon_2)
    % \bar{g}_\alpha \circ \mathcal{B}(\rho_{F_1|1}, \dots, \rho_{F_n|n}).

    methods (Static)
        function ent = conEnt(alpha,perturbation,perturbationAff,probBlocks,rho,keyProj,krausOps)
            % The conditional Renyi sandwich (tilde) down entropy
            % restricted to QKD as seen in the class definition.
            %
            % When calling the repeating arguments more than once, the
            % repeating arguments are treated as conditional with the
            % probability that the block occurs given by probBlocks.
            %
            %
            % NOTE that inputs have very few checks on them for
            % performance. Take care when using this function.
            %
            % Input:
            % * alpha: The Renyi entropy parameter. Alpha must be in the
            %   half open interval (1,2].
            % * perturbation (0): Linear perturbation for G(rho). In the
            %   range [0,1].
            % * perturbationAff (0): Affine perturbation for G(rho). In the
            %   range [0,1].
            % * probBlocks (1): Probability weights for each block
            %   occurring. Stored as a vector. The number of elements in
            %   probBlocks must be equal to the number of times the
            %   repeated arguments are provided. If probBlocks is
            %   subnormalized, the remaining weight is assumed to be mapped
            %   to a discard symbol.
            %
            % Repeating Input:
            % * rho: The density matrix shared between Alice and Bob.
            % * keyProj: A cell array of the Projection operators of the
            %   pinching channel, Z, acting on the key register used in the
            %   solver.
            % * krausOps: A cell array of the Kraus operators for the
            %   post-selection map, G, of Alice and Bob. Kraus operators
            %   that form a CPTNI map are effectively completed to a CPTP
            %   map where all remaining weight is mapped to a discard
            %   state.
            %
            % Output:
            % * ent: The entropy returned by the above expression.
            %
            % see also RenyiTidleDown, RenyiTildeDown.gradConEnt
            arguments                              
                alpha (1,1) double {mustBeInRange(alpha,1,2,"exclude-lower")}
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)} = 0
                perturbationAff (1,1) double {mustBeInRange(perturbationAff,0,1)} = 0;
                probBlocks (1,:) double {mustBeNonnegative,mustBeNonempty} = 1 %note: check sum to <= 1, check same size as number of repeating blocks
            end
            arguments(Repeating)
                rho (:,:) double {mustBeHermitian} % must be a density matrix  
                keyProj (:,1) cell % cell of projection operators
                krausOps (:,1) cell % cell of kraus ops for a CPTNI map
            end
            mustBeEqualSize(probBlocks,keyProj); %check if probBlocks and keyProj have same number of elements

            innerFunc = @(x,y,z) innerTerm(x,y,z,alpha,perturbation,perturbationAff);
            sumInnerTerms = cellfun(innerFunc,rho,keyProj,krausOps);
            sumInnerTerms = sumInnerTerms*probBlocks.' + (1 - sum(probBlocks));
            
            ent = -log2(sumInnerTerms)/(alpha-1)/(1-perturbation)/(1-perturbationAff);
        end

        function gradients = gradConEnt(alpha,perturbation,perturbationAff, probBlocks,rho,keyProj,krausOps)
            % Gradient of the conditional Renyi sandwich (tilde) down
            % conditional entropy restricted to QKD as seen in the class
            % definition.
            %
            %
            % When calling the repeating arguments more than once, the
            % repeating arguments are treated as conditional with the
            % probability that the block occurs given by probBlocks.
            % Furthermore, the gradient is returned as a cell array of each
            % block.
            %
            % Note that inputs have very few checks on them for
            % performance. Take care when using this function.
            %
            % Input:
            % * alpha: The Renyi entropy parameter. Alpha must be in the
            %   half open interval (1,2].
            % * perturbation (0): Linear perturbation for G(rho). In the
            %   range [0,1].
            % * perturbationAff (0): Affine perturbation for G(rho). In the
            %   range [0,1].
            % * probBlocks (1): Probability weights for each block
            %   occurring. Stored as a vector. The number of elements in
            %   probBlocks must be equal to the number of times the
            %   repeated arguments are provided. If probBlocks is
            %   subnormalized, the remaining weight is assumed to be mapped
            %   to a discard symbol.
            %
            % Repeating Input:
            % * rho: The density matrix shared between Alice and Bob.
            % * keyProj: A cell array of the Projection operators of the
            %   pinching channel, Z, acting on the key register used in the
            %   solver.
            % * krausOps: A cell array of the Kraus operators for the
            %   post-selection map, G, of Alice and Bob. Kraus operators
            %   that form a CPTNI map are effectively completed to a CPTP
            %   map where all remaining weight is mapped to a discard
            %   state.
            %
            % Output:
            % * gradients: Cell array of the gradients evaluated at rho.
            %
            % see also RenyiTildeDown, RenyiTildeDown.conEnt
            arguments
                alpha (1,1) double {mustBeInRange(alpha,1,2,"exclude-lower")}
                perturbation (1,1) double {mustBeInRange(perturbation,0,1)} = 0
                perturbationAff (1,1) double {mustBeInRange(perturbationAff,0,1)} = 0
                probBlocks (1,:) double {mustBeNonnegative,mustBeNonempty} = 1 %note: check sum to <= 1, check same size as number of repeating blocks
            end
            arguments(Repeating)
                rho (:,:) double {mustBeHermitian} % must be a density matrix  
                keyProj (:,1) cell % cell of projection operators
                krausOps (:,1) cell % cell of kraus ops for a CPTNI map
            end

            mustBeEqualSize(probBlocks,keyProj); %check if probBlocks and keyProj have same number of elements

            innerFunc = @(x,y,z) innerTerm(x,y,z,alpha,perturbation,perturbationAff);
            sumInnerTerms = cellfun(innerFunc,rho,keyProj,krausOps);
            sumInnerTerms = sumInnerTerms*probBlocks.' + (1 - sum(probBlocks));
            

            constOutFront = -1/(alpha-1)/sumInnerTerms/log(2)...
                /(1-perturbation)/(1-perturbationAff);

            gradients = cell(numel(probBlocks),1);
            for index = 1:numel(probBlocks)
                gradients{index} = constOutFront*probBlocks(index)*...
                    gradInnerTerm(rho{index},keyProj{index},krausOps{index},alpha,perturbation,perturbationAff);
            end
        end
    end

        
end

function val = innerTerm(rho,projOps,krausOps,alpha,perturbation,perturbationAff)
% The \theta_\alpha \circ \Psi_{\epsilon_1} \circ \Delta_{\epsilon_2} \circ
% \mathcal{G} (\rho) term.
%
% NOTE perturbation lower bound correction of 1/(1-perturbation) is handled
% in gradConEnt
arguments
    rho (:,:) double % subnormalized density matrix
    projOps (:,1) cell
    krausOps (:,1) cell % CPTNI
    alpha (1,1) double
    perturbation (1,1) double
    perturbationAff (1,1) double
end
% Computation of the value inside the log for the Renyi tilde down arrow
% conditional entropy.

gRho = ApplyMap(rho,krausOps); %G(rho)
gRho = (gRho+gRho')/2;

gRhoSafeCutOff = powmSafeCutoff(gRho,perturbation,perturbationAff);

gRho = perturbationChannel(gRho,perturbation);
gRho = affinePerturbationMap(gRho,perturbationAff);

gRhoPow = powmSafe(gRho,1/alpha,gRhoSafeCutOff);

zRho = ApplyMap(gRhoPow,projOps); %Z((G(rho))^(1/alpha))
zRho = (zRho+zRho')/2;

% don't know how to find a safe cut off for zRho^alpha beyond 0. One
% probably isn't needed though.
zRhoPow = powmSafe(zRho,alpha,0);

% fix CPTNI maps by adding the value for the discard state and affine map
termCPTNI = 1 - trace(gRho); %note gRho is already perturbed

val = real(trace(zRhoPow))+termCPTNI;
end

%% Gradient of inner term related functions

function grad = gradInnerTerm(rho,projOps,krausOps,alpha,perturbation,perturbationAff)
arguments
    rho (:,:) double % density matrix
    projOps (:,1) cell
    krausOps (:,1) cell
    alpha (1,1) double
    perturbation (1,1) double
    perturbationAff (1,1) double
end
% we need the eigenvalues and vectors of gRho
gRho = ApplyMap(rho,krausOps);
gRho = (gRho+gRho')/2;

gRhoSafeCutOff = powmSafeCutoff(gRho,perturbation,perturbationAff);

gRho = perturbationChannel(gRho,perturbation);
gRho = affinePerturbationMap(gRho,perturbationAff);

gRhoPow = powmSafe(gRho,1/alpha,gRhoSafeCutOff);

zRho = ApplyMap(gRhoPow,projOps);
zRho = (zRho+zRho')/2;

% Jack 2025/03/06: note that eig_min(Z(sigma)) >= eig_min(sigma), so we
% could tighten this safety bound 0 -> eig_min(gRhoPow).^(1/alpha), which
% we know. However, I don't think it would do much and I don't want to
% rerun the plots.
zRhoPow = powmSafe(zRho,alpha-1,0); % should have alpha-1 >=0, so 0 is be safe.

powGradTerm = powGrad(zRhoPow,gRho,1/alpha,gRhoSafeCutOff);

% affine perturbation derivative and dual
powGradTerm = (1-perturbationAff)*powGradTerm;
% linear perturbation derivative and dual
powGradTerm = perturbationChannel(powGradTerm,perturbation);

% G map derivative and dual
dualKraus = DualMap(krausOps);
powGradTerm = ApplyMap(powGradTerm,dualKraus);
powGradTerm = (powGradTerm+powGradTerm')/2;

% fix CPTNI maps by adding the value for the discard state and affine map
termCPTNI = - (1-perturbationAff)*ApplyMap(eye(size(krausOps{1},1)),dualKraus);
termCPTNI = (termCPTNI+termCPTNI')/2;
grad = alpha*powGradTerm + termCPTNI;
end


function grad = powGrad(rho,diffPoint,beta,safeCutOff)
% gradient of the matrix power function rho^beta. It sucks to compute.
arguments
    rho (:,:) double % hermitian PSD
    diffPoint (:,:) double %same size PSD
    beta (1,1) double {mustBePositive}
    safeCutOff (1,1) double {mustBePositive} % 0 doesn't work because we need beta-1 as an exponent
end

[eigBasis,eigVals] = eig(diffPoint,"vector");

% calculate the power gradient coefficient matrix. I would love a not for
% loop way of doing it.
eigVals = real(eigVals);
eigVals(eigVals <safeCutOff) = safeCutOff;

powMatCoef = zeros(size(diffPoint));
for indexI = 1:numel(eigVals)
    eigValI = eigVals(indexI);
    % diag term
    powMatCoef(indexI,indexI) = beta*(eigValI^(beta-1));

    % off diag terms
    for indexJ = indexI+1:numel(eigVals)
        eigValJ = eigVals(indexJ);

        if ismembertol(eigValI,eigValJ)
            term = beta*(eigValI^(beta-1)+eigValJ^(beta-1))/2; %average of the derivatives.
        else
            term = (eigValI^beta-eigValJ^beta)/(eigValI-eigValJ);
        end
        powMatCoef(indexI,indexJ) = term;
        powMatCoef(indexJ,indexI) = term;
    end
end

grad = eigBasis*(powMatCoef.*(eigBasis'*rho*eigBasis))*eigBasis';
grad = (grad+grad')/2; % clean up any skew hermitian parts
end

%% matrix power with safe cut offs.
function powMat = powmSafe(mat, power, safeCutOff)
arguments
    mat (:,:) double {mustBeHermitian}
    power (1,1) double
    safeCutOff (1,1) double {mustBeNonnegative}
end

[eigVecs,eigVals] = eig(mat,'vector');
eigVals = real(eigVals); % this shouldn't be required
eigVals(eigVals<safeCutOff)=safeCutOff;

powMat = eigVecs * diag(eigVals.^power) * eigVecs';
powMat = (powMat+powMat')/2;

end

function safeCutOff = powmSafeCutoff(mat,perturbation,perturbationAff)
arguments
    mat (:,:) double
    perturbation (1,1) double
    perturbationAff (1,1) double
end
% Assumes the linear perturbation is applied first. Assumes that the matrix
% should be positive semidefinite, and any deviation is caused by numerical
% stability and round off errors.

dim = size(mat,1);
% linear perturbation
minEig = max(min(eig(mat)),0);
traceMat = max(trace(mat),0);
safeCutOff = (1-perturbation)*minEig + perturbation*traceMat/dim;
% and add the affine perturbation
safeCutOff = (1-perturbationAff)*safeCutOff + perturbationAff/dim;
end

function rhoPrime = affinePerturbationMap(rho,perturbationAff)
arguments
    rho (:,:) double % manually check that it's positive semi-definite
    perturbationAff (1,1) double {mustBeInRange(perturbationAff,0,1)}
end
dim = size(rho, 1);
rhoPrime = (1-perturbationAff) * rho + perturbationAff * eye(dim)/dim;
rhoPrime = (rhoPrime + rhoPrime')/2;
end