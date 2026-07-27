function K = asymptotic_qubit_bb84(px, f_EC, depolarization, transmittance)
    % ASYMPTOTIC_QUBIT_BB84 Computes the asymptotic BB84 key rate 
    % including basis probabilities, error correction efficiency, 
    % isotropic depolarization, and channel loss (transmittance).
    %
    % Usage:
    %   K = asymptotic_qubit_bb84(px, f_EC, depolarization, transmittance)
    %   transmittance: Probability a pulse reaches Bob and is detected (0 to 1).

    if nargin < 4
        transmittance = 1.0; % Default to no loss if not provided
    end

    % 1. Determine e_x and e_z based on standard isotropic depolarization
    ex = depolarization / 2;
    ez = depolarization / 2;
    
    % 2. Calculate basis probabilities
    pz = 1 - px;
    
    % 3. Calculate binary Shannon entropies
    hx = binary_entropy(ex);
    hz = binary_entropy(ez);
    
    % 4. Calculate secure key rates for each basis stream independently
    % For each basis: Reconcile bit errors (costs f_EC * h_bit), 
    % bound phase errors (costs h_phase).
    rate_x = (px.^2) .* (1 - f_EC * hx - hz);
    rate_z = (pz.^2) .* (1 - f_EC * hz - hx);
    
    % 5. Enforce non-negativity (a compromised basis stream yields 0 key)
    rate_x(rate_x < 0) = 0;
    rate_z(rate_z < 0) = 0;
    
    % 6. Total asymptotic key rate per emitted pulse, scaled by transmittance
    K = transmittance .* (rate_x + rate_z);
end

function H = binary_entropy(e)
    % BINARY_ENTROPY Calculates the binary Shannon entropy
    
    H = zeros(size(e));
    
    % Compute entropy only for valid probabilities strictly between 0 and 1
    valid_idx = (e > 0) & (e < 1);
    e_val = e(valid_idx);
    
    H(valid_idx) = -e_val .* log2(e_val) - (1 - e_val) .* log2(1 - e_val);
end