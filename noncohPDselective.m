function globalPD = noncohPDselective(M, PFA_local, globalThreshold, SNR)
% Computes the global detection probability
% M             : number of local sensors
% PFA_local     : local false alarm probability (same for all sensors)
% PFA_global    : desired global false alarm probability
% SNR           : linear scale SNR

theta = 1./(1 + SNR);

% Define the local threshold
localThreshold = -log(PFA_local);

% Local probability of detection
PD_local = PFA_local.^theta;

% Define the total CDF function F_T(gamma)
globalPD = sum(arrayfun(@(m) binomial_term(M, m, PD_local)*conditional_cdf_gamma_shifted(globalThreshold, m, localThreshold, theta), 0 : M));

end

% ---- Helper: Binomial PMF term ----
function b = binomial_term(n, m, p)
    b = nchoosek(n, m)*p^m*(1 - p)^(n - m);
end

% ---- Helper: Conditional CDF using regularized upper incomplete gamma ----
function cdf_val = conditional_cdf_gamma_shifted(gamma, m, lambda, theta)
    if m == 0
        cdf_val = double(gamma <= 0); % Degenerate at T = 0
        return;
    end

    % Shifted gamma CDF: F_T|M=m(t) = Q(m, (t - m*lambda)/theta)
    if gamma < m*lambda
        cdf_val = 0;
    else
        x = (gamma - m*lambda)./theta;
        cdf_val = gammainc(x, m, 'upper'); % MATLAB's Q(m, x)
    end
end