classdef RenyiTildeDown
    % Class containing functions to implement the conditional renyi tidle
    % down entropy For a set of CPTNI kraus operators and projection map.
    % Specifically, the CPTNI map and projection map are the G and Z maps
    % from https://arxiv.org/abs/1710.05511v2.

    methods (Static)
        function ent = conEnt(alpha,perturbation,perturbationAff,probBlocks,rho,keyProj,krausOps)
            % conEnt Renyi tilde down arrow conditional entropy for the key
            % register R given Eve's quantum system E and the classical
            % announcements C. Here we rexpressed it in terms of the G and
            % Z maps acting on Alice and Bob's marginal state $\rho_{AB}$.
            %
            % The expression is given by:
            % $\tilde{H}^\uparrow_\alpha(R|EC)_{\ket{\rho}}
            % = -\bar{H}^\downarrow_{1/\alpha}(R|R'ABC')_{\ket{\rho}}
            % = -\frac{1}{\alpha-1} \log(\Tr[ \text{Pow}_\alpha \circ
            %   \mathcal{Z} \circ \text{Pow}_{1/\alpha} \circ \mathcal{G}
            %   (\rho) + I - \mathcal{G}(\rho) ])$
            %
            %NOTE: FIX DESCRIPTION!!!!
            %
            % Note that inputs have very few checks on them for
            % performance. Take care when using this function.
            %
            % Input:
            % * rho: The density matrix shared between Alice and Bob.
            % * keyProj: A cell array of the Projection operators of the
            %   pinching channel, Z, acting on the key register used in the
            %   solver.
            % * krausOps: A cell array of the Kraus operators for the
            %   post-selection map, G, of Alice and Bob. Kraus operators
            %   that form a CPTNI map are effectively completed to a CPTP
            %   map where all remaining weight is mapped to a discard state
            %   (The $\rho - \mathcal{G}(\rho)$ term in the formula above).
            % * alpha: The Renyi entropy parameter. Alpha must be in the
            %   half open interval (1,2].
            %
            % Output:
            % * ent: The entropy returned by the above expression.
            %
            % see also RenyiTildeDown.gradConEnt
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
            % gradConEnt gradient of the Renyi tilde down arrow conditional
            % entropy for the key register R given Eve's quantum system E
            % and the classical announcements C. Here we rexpressed it in
            % terms of the G and Z maps acting on Alice and Bob's marginal
            % state $\rho_{AB}$.
            %
            % Gradient of this:
            % $\tilde{H}^\uparrow_\alpha(R|EC)_{\ket{\rho}}
            % = -\bar{H}^\downarrow_{1/\alpha}(R|R'ABC')_{\ket{\rho}}
            % = -\frac{1}{\alpha-1} \log(\Tr[ \text{Pow}_\alpha \circ
            %   \mathcal{Z} \circ \text{Pow}_{1/\alpha} \circ \mathcal{G}
            %   (\rho) + I - \mathcal{G}(\rho) ])$
            %
            % Note that inputs have very few checks on them for
            % performance. Take care when using this function.
            %
            % Input:
            % * rho: The density matrix shared between Alice and Bob.
            % * keyProj: A cell array of the Projection operators of the
            %   pinching channel, Z, acting on the key register used in the
            %   solver.
            % * krausOps: A cell array of the Kraus operators for the
            %   post-selection map, G, of Alice and Bob. Kraus operators
            %   that form a CPTNI map are effectively completed to a CPTP
            %   map where all remaining weight is mapped to a discard state
            %   (The $\rho - \mathcal{G}(\rho)$ term in the formula above).
            % * alpha: The Renyi entropy parameter. Alpha must be in the
            %   half open interval (1,2].
            %
            % Output:
            % * grad: Gradient of the expression above with respect to
            %   $\rho$.
            %
            % see also RenyiTildeDown.gradConEnt
            % Derivative of this:
            % \tilde{H}^\uparrow_\alpha(R|EC)_{\ket{\rho}}
            % = -\bar{H}^\downarrow_{1/\alpha}(R|R'ABC')_{\ket{\rho}}
            % = -\frac{1}{\alpha-1} \log(\Tr[ \text{Pow}_\alpha \circ
            %   \mathcal{Z} \circ \text{Pow}_{1/\alpha} \circ \mathcal{G}
            %   (\rho) + \rho - \mathcal{G}(\rho) ])
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
            %****** Now outputs cell of gradients. Fix output call everywhere!!!!*****
            % gradients = constOutFront*( gradSumInnerTerms ); % moved alpha from constOutFront to only the one term.
        end
    end

        
end

function val = innerTerm(rho,projOps,krausOps,alpha,perturbation,perturbationAff)

% perturbation lower bound correction of 1/(1-perturbation) is handled in
% gradConEnt
arguments
    rho (:,:) double % density matrix
    projOps (:,1) cell
    krausOps (:,1) cell
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
zRhoPow = powmSafe(zRho,alpha-1,0); % should have alpha-1 >=0, so 0 should be safe.

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