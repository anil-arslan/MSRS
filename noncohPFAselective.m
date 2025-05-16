function globalThreshold = noncohPFAselective(M, PFA_local, PFA_global)
% Computes the global threshold
% M             : number of local sensors
% PFA_local     : local false alarm probability (same for all sensors)
% PFA_global    : desired global false alarm probability

% Define the local threshold
localThreshold = -log(PFA_local);

% Define the total CDF function F_T(gamma)
F_T = @(t) sum(arrayfun(@(m) binomial_term(M, m, PFA_local)*conditional_cdf_gamma_shifted(t, m, localThreshold), 0 : M));

% Use root-finding to solve F_T(gamma) = 1 - PFA_global
% We search gamma in [0, globalThresholdMax]
globalThresholdMax = M*localThreshold + gammaincinv(PFA_global, M, 'upper'); % Safe upper bound
% it should be bigger than M*localThreshold beacuse if all values are
% passed then the summation will be bigger than that

% t = linspace(0, globalThresholdMax, 1000);
% T = zeros(1, length(t));
% for i = 1 : length(t)
%     T(i) = F_T(t(i)) - (1 - PFA_global);
% end
% figure; loglog(t, abs(T));

globalThreshold = fzero(@(t) F_T(t) - (1 - PFA_global), [0, globalThresholdMax]);

end

% ---- Helper: Binomial PMF term ----
function b = binomial_term(n, m, p)
    b = nchoosek(n, m)*p^m*(1 - p)^(n - m);
end

% ---- Helper: Conditional CDF using regularized upper incomplete gamma ----
function cdf_val = conditional_cdf_gamma_shifted(t, m, lambda)
    if m == 0
        cdf_val = double(t >= 0); % Degenerate at T = 0
        return;
    end

    % Shifted gamma CDF: F_T|M=m(t) = 1 - Q(m, t - m*lambda)
    if t < m*lambda
        cdf_val = 0;
    else
        x = t - m*lambda;
        cdf_val = 1 - gammainc(x, m, 'upper'); % MATLAB's Q(m, x)
    end
end