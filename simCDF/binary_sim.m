clc; clear;
M = 3;
K = 1 : M + 1;
pfaGlobal = 1e-6;
pfaLocal = logspace(0, -6, 1001);
rho = 10;

Pfa = binaryCCDF(M, K, pfaLocal, 0);
Pd = binaryCCDF(M, K, pfaLocal, rho);

q = zeros(M, length(pfaLocal));
for i = 2 : 4
    q(i - 1, :) = (pfaGlobal - Pfa(i, :))./(Pfa(i - 1, :) - Pfa(i, :));
end
q(q > 1 | q < 0) = nan;
pfaNew = q.*Pfa(1 : M, :) + (1 - q).*Pfa((1 : M) + 1, :);
pdNew = q.*Pd(1 : M, :) + (1 - q).*Pd((1 : M) + 1, :);

%%% Visualization
figure(1);
subplot(2, 1, 1);
loglog(pfaLocal, Pfa(1 : M, :));
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of False Alarm');
leg = legend(num2str([1; 2; 3]));
title(leg, 'M-out-of-K');
grid off; grid on; grid minor;

subplot(2, 1, 2);
semilogx(pfaLocal, q);
xlabel('Local Probability of False Alarm');
ylabel('Probability of Threshold Randomization');
leg = legend(num2str([1; 2; 3]));
title(leg, 'M-out-of-K');
grid off; grid on; grid minor;

figure(2);
subplot(2, 1, 1);
semilogx(pfaLocal, Pd(1 : M, :));
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of Detection');
leg = legend(num2str([1; 2; 3]));
title(leg, 'M-out-of-K');
grid off; grid on; grid minor;

subplot(2, 1, 2);
semilogx(pfaLocal, pdNew);
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of Detection');
leg = legend(num2str([1; 2; 3]));
title(leg, 'M-out-of-K');
grid off; grid on; grid minor;

function Pd_global = binaryCCDF(M, K, Pfa_local, rho)
    % Pd_binary returns the global Pd for a k-out-of-M binary combining rule.
    % M           : number of sensors
    % k           : rule parameter (k-out-of-M)
    % Pfa_local   : local false alarm probability
    % rho         : SNR (linear)
    
    % Local detection probability as a function of SNR and Pfa_local
    Pd_local = (Pfa_local).^(1./(1 + rho));
    
    % Global Pd using binomial summation
    Pd_global = zeros(length(K), length(Pfa_local));
    for k = K
        for i = k : M
            Pd_global(k, :) = Pd_global(k, :) + nchoosek(M, i).*Pd_local.^i.*(1 - Pd_local).^(M - i);
        end
    end
end