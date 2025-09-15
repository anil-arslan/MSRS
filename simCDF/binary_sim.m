clc; clear;

kayit = 1;

M = 9;
pfaGlobal = 1e-6;
pfaLocal = logspace(0, -8, 10001);
rho = 8; % dB

[pfaNew, q_local] = binaryProbability(M, pfaLocal, pfaGlobal, -inf);
pdNew = binaryProbability(M, pfaLocal, pfaGlobal, rho);


%%% Visualization

fig = figure(1);
semilogx(pfaLocal, q_local, 'LineWidth', 2);
ax = gca;
ax.Children(end - 7).Color = [1 0 1];
ax.Children(end - 8).Color = [0 0.5 0];
ax.Children(end - 9).Color = [0 0 0.5];
xlabel('Local Probability of False Alarm');
ylabel('Probability of Threshold Reduction');
leg = legend(num2str((1 : (M + 1))'), 'Location', 'best');
title(leg, 'K-out-of-9');
grid off; grid on; grid minor;

if kayit
    figureName = 'probability_of_threshold_randomization';
    savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
end

fig = figure(2);
semilogx(pfaLocal, pdNew, 'LineWidth', 2);
ax = gca;
ax.Children(end - 7).Color = [1 0 1];
ax.Children(end - 8).Color = [0 0.5 0];
ax.Children(end - 9).Color = [0 0 0.5];
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of Detection');
leg = legend(num2str((1 : (M + 1))'), 'Location', 'best');
title(leg, 'K-out-of-9');
grid off; grid on; grid minor;

if kayit
    figureName = 'binary_combining_pd_operating_final';
    savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
end
