clc; clear;

kayit = 0;

M = 3;
K = 1 : M + 1;
pfaGlobal = 1e-6;
pfaLocal = logspace(0, -6, 1001);
rho = 10^(8/10);

Pfa = binaryCCDF(M, K, pfaLocal, 0);
Pd = binaryCCDF(M, K, pfaLocal, rho);

q = zeros(M, length(pfaLocal));
for i = 2 : M + 1
    q(i - 1, :) = (pfaGlobal - Pfa(i, :))./(Pfa(i - 1, :) - Pfa(i, :));
end
q(q > 1 | q < 0) = nan;
pfaNew = q.*Pfa(1 : M, :) + (1 - q).*Pfa((1 : M) + 1, :);
pdNew = q.*Pd(1 : M, :) + (1 - q).*Pd((1 : M) + 1, :);

%%% Visualization

fig = figure(1);
loglog(pfaLocal, Pfa(1 : M, :), 'LineStyle', '--'); hold on;
% ax = gca;
% ax.Children(end - 7).Color = [1 0 1];
% ax.Children(end - 8).Color = [0 0.5 0];
Pfa1 = Pfa(1 : M, :);
Pfa1(isnan(q)) = nan;
loglog(pfaLocal, Pfa1, 'LineWidth', 2, 'Color', 'blue');
Pfa2 = Pfa(2 : M + 1, :);
Pfa2(isnan(q)) = nan;
loglog(pfaLocal, Pfa2, 'LineWidth', 2, 'Color', 'red');
ylim([1e-12, 1]);
yline(1e-6, 'LineStyle', '--');
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of False Alarm');
leg = legend(num2str((1 : M)'), 'Location', 'northwest');
title(leg, 'M-out-of-K');
grid off; grid on; grid minor;

if kayit
    figureName = 'binary_combining_pfa_operating';
    savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
end

fig = figure(2);
semilogx(pfaLocal, Pd(1 : M, :), 'LineStyle', '--'); hold on;
% ax = gca;
% ax.Children(end - 7).Color = [1 0 1];
% ax.Children(end - 8).Color = [0 0.5 0];
Pd1 = Pd(1 : M, :);
Pd1(isnan(q)) = nan;
loglog(pfaLocal, Pd1, 'LineWidth', 2, 'Color', 'blue');
Pd2 = Pd(2 : M + 1, :);
Pd2(isnan(q)) = nan;
loglog(pfaLocal, Pd2, 'LineWidth', 2, 'Color', 'red');
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of Detection');
leg = legend(num2str((1 : M)'), 'Location', 'northwest');
title(leg, 'M-out-of-K');
grid off; grid on; grid minor;

if kayit
    figureName = 'binary_combining_pd_operating';
    savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
end

fig = figure(3);
semilogx(pfaLocal, q, 'LineWidth', 2);
% ax = gca;
% ax.Children(end - 7).Color = [1 0 1];
% ax.Children(end - 8).Color = [0 0.5 0];
xlabel('Local Probability of False Alarm');
ylabel('Probability of Threshold Reduction');
leg = legend(num2str((1 : M)'), 'Location', 'best');
title(leg, 'M-out-of-K');
grid off; grid on; grid minor;

if kayit
    figureName = 'probability_of_threshold_randomization';
    savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
end

fig = figure(4);
semilogx(pfaLocal, pdNew, 'LineWidth', 2);
% ax = gca;
% ax.Children(end - 7).Color = [1 0 1];
% ax.Children(end - 8).Color = [0 0.5 0];
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of Detection');
leg = legend(num2str((1 : M)'), 'Location', 'best');
title(leg, 'M-out-of-K');
grid off; grid on; grid minor;

if kayit
    figureName = 'binary_combining_pd_operating_final';
    savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
end

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