clc; clear;

kayit = 0;

M = 9;
K = 0 : M + 1;
pfaGlobal = 1e-6;
pfaLocal = logspace(0, -8, 1001);
rho = 10^(8/10);

Pfa = binaryCCDF(M, K, pfaLocal, 0);
Pd = binaryCCDF(M, K, pfaLocal, rho);





rulesPfaLocal = zeros(1, M);
for k = 1 : M
    rulesPfaLocal(k) = fzero(@(lambda) binaryCCDF(M, k, lambda, 0) - 1e-6, [0 1]);
end
rulesPfaLocal = [0 rulesPfaLocal 1];

q_local = zeros(M + 1, length(pfaLocal));
for i = 2 : M + 2
    q_local(i - 1, :) = (pfaLocal - rulesPfaLocal(i))./(rulesPfaLocal(i - 1) - rulesPfaLocal(i));
end
q_local(q_local > 1 | q_local < 0) = nan;




Pfa = binaryCCDF(M, K, rulesPfaLocal, 0);
Pd = binaryCCDF(M, K, rulesPfaLocal, rho);



q = zeros(M + 1, length(pfaLocal));
for i = 2 : M + 2
    q(i - 1, :) = (pfaGlobal - q_local(i - 1, :).*Pfa(i, i - 1) - (1 - q_local(i - 1, :)).*Pfa(i, i))./(q_local(i - 1, :).*Pfa(i - 1, i - 1) + (1 - q_local(i - 1, :)).*Pfa(i - 1, i) - q_local(i - 1, :).*Pfa(i, i - 1) - (1 - q_local(i - 1, :)).*Pfa(i, i));
end
q(q > 1 | q < 0) = nan;

idx1 = (1 : M + 2 : (M + 1)*(M + 2)) + (0 : M);
idx2 = (2 : M + 2 : (M + 1)*(M + 2)) + (0 : M);
idx3 = (1 : M + 2 : (M + 1)*(M + 2)) + (0 : M) + (M + 2);
idx4 = (2 : M + 2 : (M + 1)*(M + 2)) + (0 : M) + (M + 2);

pfaNew = q.*q_local.*Pfa(idx1).' + (1 - q).*q_local.*Pfa(idx2).' + ...
q.*(1 - q_local).*Pfa(idx3).' + (1 - q).*(1 - q_local).*Pfa(idx4).';

pdNew = q.*q_local.*Pd(idx1).' + (1 - q).*q_local.*Pd(idx2).' + ...
q.*(1 - q_local).*Pd(idx3).' + (1 - q).*(1 - q_local).*Pd(idx4).';

%%% Visualization

% fig = figure(1);
% loglog(pfaLocal, Pfa(2 : M + 1, :), 'LineStyle', '--'); hold on;
% ax = gca;
% ax.Children(end - 7).Color = [1 0 1];
% ax.Children(end - 8).Color = [0 0.5 0];
% Pfa1 = Pfa(2 : M + 1, :);
% Pfa1(isnan(q(2 : end, :))) = nan;
% loglog(pfaLocal, Pfa1, 'LineWidth', 2, 'Color', 'blue');
% Pfa2 = Pfa(2 : M + 1, :);
% Pfa2(isnan(q(1 : end - 1, :))) = nan;
% loglog(pfaLocal, Pfa2, 'LineWidth', 2, 'Color', 'red');
% ylim([1e-12, 1]);
% yline(1e-6, 'LineStyle', '--');
% xlabel('Local Probability of False Alarm');
% ylabel('Global Probability of False Alarm');
% leg = legend(num2str((1 : M)'), 'Location', 'northwest');
% title(leg, 'K-out-of-9');
% grid off; grid on; grid minor;

if kayit
    figureName = 'binary_combining_pfa_operating';
    savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
end

fig = figure(2);
semilogx(pfaLocal, Pd(2 : M + 1, :), 'LineStyle', '--'); hold on;
ax = gca;
ax.Children(end - 7).Color = [1 0 1];
ax.Children(end - 8).Color = [0 0.5 0];
Pd1 = Pd(2 : M + 1, :);
Pd1(isnan(q(2 : end, :))) = nan;
loglog(pfaLocal, Pd1, 'LineWidth', 2, 'Color', 'blue');
Pd2 = Pd(2 : M + 1, :);
Pd2(isnan(q(1 : end - 1, :))) = nan;
loglog(pfaLocal, Pd2, 'LineWidth', 2, 'Color', 'red');
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of Detection');
leg = legend(num2str((1 : M)'), 'Location', 'northwest');
title(leg, 'K-out-of-9');
grid off; grid on; grid minor;

if kayit
    figureName = 'binary_combining_pd_operating';
    savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
end

fig = figure(3);
semilogx(pfaLocal, q, 'LineWidth', 2);
ax = gca;
ax.Children(end - 7).Color = [1 0 1];
ax.Children(end - 8).Color = [0 0.5 0];
ax.Children(end - 9).Color = [0 0 0.5];
xlabel('Local Probability of False Alarm');
ylabel('Probability of Threshold Reduction');
leg = legend(num2str((1 : M + 1)'), 'Location', 'best');
title(leg, 'K-out-of-9');
grid off; grid on; grid minor;

% if kayit
%     figureName = 'probability_of_threshold_randomization';
%     savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
%     saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
% end

fig = figure(4);
semilogx(pfaLocal, pfaNew, 'LineWidth', 2);
ax = gca;
ax.Children(end - 7).Color = [1 0 1];
ax.Children(end - 8).Color = [0 0.5 0];
ax.Children(end - 9).Color = [0 0 0.5];
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of Detection');
leg = legend(num2str((1 : M + 1)'), 'Location', 'best');
title(leg, 'K-out-of-9');
grid off; grid on; grid minor;

% if kayit
%     figureName = 'binary_combining_pd_operating_final';
%     savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
%     saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
% end





fig = figure(5);
semilogx(pfaLocal, q_local, 'LineWidth', 2);
ax = gca;
ax.Children(end - 7).Color = [1 0 1];
ax.Children(end - 8).Color = [0 0.5 0];
ax.Children(end - 9).Color = [0 0 0.5];
xlabel('Local Probability of False Alarm');
ylabel('Probability of Threshold Reduction');
leg = legend(num2str((1 : M + 1)'), 'Location', 'best');
title(leg, 'K-out-of-9');
grid off; grid on; grid minor;