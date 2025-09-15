% Anıl Arslan 2303980
addpath(genpath('C:/GitRepo/MSRS/library'));
addpath(genpath('C:/GitRepo/MSRS/figureGenerator'));
addpath(genpath('C:/GitRepo/MSRS/simCDF'));
addpath(genpath('/Users/anilarslan/Desktop/MSRS/library'));
addpath(genpath('/Users/anilarslan/Desktop/MSRS/figureGenerator'));
addpath(genpath('/Users/anilarslan/Desktop/MSRS/simCDF'));

%% Centralized Equal SNR, SW 2 PD vs SNR
clc; clear; % close all;
colors = rgb2hex([ ...
    0.0000    0.4470    0.7410 % blue
    0.8500    0.3250    0.0980 % orange
    0.9290    0.6940    0.1250 % yellow
    0.4940    0.1840    0.5560 % magenta
    0.4660    0.6740    0.1880 % green
    0.3010    0.7450    0.9330 % cyan
    0.6350    0.0780    0.1840]);
kayit = 0;
set(groot, "defaultFigurePosition", [680 458 560 420]);


globalPFA = 1e-6;
snrs = linspace(-5, 20, 10001);
M = [1 2 4 8 16 32 64];

fig = figure;
dtc = detector( ...
    "globalPFA", globalPFA, ...
    "numberOfSensors", M, ...
    "SNR_input_dB", snrs, ...
    "localPFA", 1 ...
    );
dtc.setalgorithm( ...
    "signalPhaseModel", "decorrelatedUniform", ...
    "signalAmplitudeModel", "decorrelatedExponential", ...
    "globalFusionRule", "SLC" ...
    );
pd = squeeze(dtc.globalPD);
plot(snrs, pd, 'LineWidth', 2); hold on;

grid on; grid minor;
ylim([0, 1]);
xlabel('Average SNR per Receiving Node (dB)');
ylabel('Global Probability of Detection');
legendStr = num2str(M.');
leg = legend(legendStr, 'Location', 'best');
title(leg, 'M');

if kayit
    figureName = 'analysis_centralized_equal_snr_iid';
    savefig(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end

%% Centralized Different SNR, SW 2 PD vs SNR
clc; clear; % close all;
kayit = 0;
set(groot, "defaultFigurePosition", [680 458 560 420]);

globalPFA = 1e-6;
algorithms = ["SLC", "WSLC"];
% snrs = -5 : .02 : 20;
snrs = linspace(-5, 20, 26);
M = 9;

nmc = 100;
rng(1975);
pd = zeros(length(snrs), length(algorithms), nmc);
for mcID = 1 : nmc
    dtc = detector( ...
        "globalPFA", globalPFA, ...
        "numberOfSensors", M, ...
        "SNR_input_dB", snrs, ...
        "localPFA", 1 ...
        );
    algorithmID = 0;
    dtc.randomSNRwithFixedAverage;
    for algorithm = algorithms
        algorithmID = algorithmID + 1;
        % rng(1975); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 20);
        % rng(1976); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 10);
        % rng(1975); dtc.randomSNRwithFixedAverage;
        dtc.setalgorithm( ...
            "signalPhaseModel", "decorrelatedUniform", ...
            "signalAmplitudeModel", "decorrelatedExponential", ...
            "binaryDetectionConstraint", "fixedGlobal|LocalPFA", ...
            "globalFusionRule", algorithm ...
            );
        pd(:, algorithmID, mcID) = squeeze(dtc.globalPD);
    end
    % pd = squeeze(dtc.globalPD);
    % switch algorithmID
    %     case 1
    %         lineStyle = '-.';
    %     case 2
    %         lineStyle = '-';
    % end
    % plot(snrs, pd, 'LineWidth', 2, 'LineStyle', lineStyle); hold on;
    fprintf('mcID = %d/%d\n', mcID, nmc);
end
save('C:\GitRepo\MSRS\figureGenerator\differentSNR_PD_centralized_SLC.mat', 'pd');
%%
fig1 = figure;
colors = rgb2hex([ ...
    0.0000    0.4470    0.7410 % blue
    0.8500    0.3250    0.0980 % orange
    0.9290    0.6940    0.1250 % yellow
    0.4940    0.1840    0.5560 % magenta
    0.4660    0.6740    0.1880 % green
    0.3010    0.7450    0.9330 % cyan
    0.6350    0.0780    0.1840]);
plot(snrs, nan(length(snrs), 1), '-k'); hold on;
plot(snrs, nan(length(snrs), 1), '-.k');
plot(snrs, nan(length(snrs), 1), '--k');
plot(snrs, mean(pd, 3), 'LineWidth', 2); hold on;
plot(snrs, mean(pd, 3) + std(pd, [], 3), 'LineWidth', 2, 'LineStyle', '-.');
plot(snrs, mean(pd, 3) - std(pd, [], 3), 'LineWidth', 2, 'LineStyle', '--');
ax = gca;
ax.Children(end - 0).Color = [0 0 0];
ax.Children(end - 1).Color = [0 0 0];
ax.Children(end - 2).Color = [0 0 0];

ax.Children(end - 3).Color = colors(1, :);
ax.Children(end - 4).Color = colors(2, :);

ax.Children(end - 5).Color = colors(1, :);
ax.Children(end - 6).Color = colors(2, :);

ax.Children(end - 7).Color = colors(1, :);
ax.Children(end - 8).Color = colors(2, :);
grid on; grid minor;
ylim([0, 1]);
xlabel('Reference Average SNR (dB)');
ylabel('Global Probability of Detection');
legend(["E[P_D^{global}]", "E[P_D^{global}] + std(P_D^{global})", "E[P_D^{global}] - std(P_D^{global})", "SLC", "WSLC"], ...
    'Location', 'best');

if kayit
    figureName = 'analysis_centralized_different_snr_iid';
    savefig(fig1, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig1, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end










%% Decentralized Equal SNR, SW 2 PD vs PFAlocal
clc; clear; % close all;
kayit = 1;
set(groot, "defaultFigurePosition", [680 458 560 420]);
colors = rgb2hex([ ...
    0.0000    0.4470    0.7410 % blue
    0.8500    0.3250    0.0980 % orange
    0.9290    0.6940    0.1250 % yellow
    0.4940    0.1840    0.5560 % magenta
    0.4660    0.6740    0.1880 % green
    0.3010    0.7450    0.9330 % cyan
    0.6350    0.0780    0.1840]);

globalPFA = 1e-6;
pfaLocal = logspace(0, -8, 4001);
algorithms = ["SLC", "BC"];
M = [1 2 3 4 5 7 9];

fig1 = figure;
algorithmID = 0;
for algorithm = algorithms
    algorithmID = algorithmID + 1;
    for sensorID = 1 : length(M)
        dtc = detector( ...
            "globalPFA", globalPFA, ...
            "numberOfSensors", M(sensorID), ...
            "SNR_input_dB", 8, ...
            "localPFA", pfaLocal ...
            );
        dtc.setalgorithm( ...
            "binaryDetectionRule", "notSpecified", ...
            "binaryDetectionConstraint", "fixedGlobal|LocalPFA", ...
            "binaryDetectionConstraint", "timeSharing", ...
            "signalPhaseModel", "decorrelatedUniform", ...
            "signalAmplitudeModel", "decorrelatedExponential", ...
            "globalFusionRule", algorithm ...
            );
        switch algorithmID
            case 1
                lineStyle = '-';
            case 2
                lineStyle = '-.';
        end
        pd = squeeze(dtc.globalPD);
        figure(fig1); semilogx(pfaLocal, pd, 'LineWidth', 2, 'LineStyle', lineStyle); hold on;
        % figure(40); semilogx(pfaLocal, squeeze(dtc.globalRandomizationProbability), 'LineWidth', 2, 'LineStyle', lineStyle, 'color', colors(sensorID)); hold on;
    end
end

figure(fig1);
grid on; grid minor;
ylim([0, 1]);
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of Detection');
legendStr = num2str(M.');
leg = legend(legendStr, 'Location', 'best');
title(leg, 'M');

if kayit
    figureName = 'square_law_combining';
    saveas(fig1, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
    savefig(fig1, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
end

%% Decentralized Different SNR, SW 2 PD vs PFAlocal
clc; clear; % close all;
colors = rgb2hex([ ...
    0.0000    0.4470    0.7410 % blue
    0.8500    0.3250    0.0980 % orange
    0.9290    0.6940    0.1250 % yellow
    0.4940    0.1840    0.5560 % magenta
    0.4660    0.6740    0.1880 % green
    0.3010    0.7450    0.9330 % cyan
    0.6350    0.0780    0.1840]);
kayit = 1;
set(groot, "defaultFigurePosition", [680 458 560 420]);

globalPFA = 1e-6;
pfaLocal = logspace(0, -8, 81);
algorithms = ["SLC", "BC", "WSLC"];
algorithms = ["SLC", "WSLC", "BC"];
M = 9;
snr = 8;

nmc = 100;
rng(1975);
pd = zeros(length(pfaLocal), length(algorithms), nmc);
for mcID = 1 : nmc
    dtc = detector( ...
        "globalPFA", globalPFA, ...
        "numberOfSensors", M, ...
        "SNR_input_dB", snr, ...
        "localPFA", pfaLocal ...
        );
    % rng(1975); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 20);
    % rng(1976); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 10);
    dtc.randomSNRwithFixedAverage;
    %%% Ekle
    % options.distributionType (1, 1) string {ismember(options.distributionType, ["loguniform", "lognormal"])} = "loguniform"
    algorithmID = 0;
    for algorithm = algorithms
        algorithmID = algorithmID + 1;
        dtc.setalgorithm( ...
            "binaryDetectionRule", "notSpecified", ...
            "binaryDetectionConstraint", "fixedGlobal|LocalPFA", ...
            "binaryDetectionConstraint", "timeSharing", ...
            "signalPhaseModel", "decorrelatedUniform", ...
            "signalAmplitudeModel", "decorrelatedExponential", ...
            "globalFusionRule", algorithm ...
            );
        pd(:, algorithmID, mcID) = squeeze(dtc.globalPD);
        % T{algorithmID} = squeeze(dtc.globalThreshold);
        % w{algorithmID} = squeeze(dtc.fusionWeights{1});
    end
    fprintf('mcID = %d/%d\n', mcID, nmc);
end
fig = figure;
semilogx(pfaLocal, nan(length(pfaLocal), 1), '-k'); hold on;
semilogx(pfaLocal, nan(length(pfaLocal), 1), '-.k');
semilogx(pfaLocal, nan(length(pfaLocal), 1), '--k');
semilogx(pfaLocal, mean(pd, 3), 'LineWidth', 2); hold on;
semilogx(pfaLocal, mean(pd, 3) + std(pd, [], 3), 'LineWidth', 2, 'LineStyle', '-.');
semilogx(pfaLocal, mean(pd, 3) - std(pd, [], 3), 'LineWidth', 2, 'LineStyle', '--');
ax = gca;
ax.Children(end - 0).Color = [0 0 0];
ax.Children(end - 1).Color = [0 0 0];
ax.Children(end - 2).Color = [0 0 0];

ax.Children(end - 3).Color = colors(1, :);
ax.Children(end - 4).Color = colors(2, :);
ax.Children(end - 5).Color = colors(3, :);

ax.Children(end - 6).Color = colors(1, :);
ax.Children(end - 7).Color = colors(2, :);
ax.Children(end - 8).Color = colors(3, :);

ax.Children(end - 9).Color = colors(1, :);
ax.Children(end - 10).Color = colors(2, :);
ax.Children(end - 11).Color = colors(3, :);
% ax = gca;
% ax.Children(end).Color = colors(1, :);
% ax.Children(end - 1).Color = colors(2, :);
% ax.Children(end - 2).Color = colors(3, :);
% ax.Children(end - 3).Color = colors(4, :);
% ax.Children(end - 4).Color = colors(1, :);
% ax.Children(end - 5).Color = colors(2, :);
% ax.Children(end - 6).Color = colors(3, :);
% ax.Children(end - 7).Color = colors(4, :);
% ax.Children(end - 8).Color = colors(1, :);
% ax.Children(end - 9).Color = colors(2, :);
% ax.Children(end - 10).Color = colors(3, :);
% ax.Children(end - 11).Color = colors(4, :);
grid on; grid minor;
ylim([0, 1]);
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of Detection');
% legend(["SLC", "BC", "WSLC", "WBC"], 'Location', 'best');
legend(["E[P_D^{global}]", "E[P_D^{global}] + std(P_D^{global})", "E[P_D^{global}] - std(P_D^{global})", "SLC", "WSLC", "BC"], ...
    'Location', 'best');

% figure; semilogx(pfaLocal(2 : end), w{4}(:, 2 : end).', 'LineWidth', 2);
% grid on; grid minor;
% xlim tight;
% xlabel('Local Probability of False Alarm');
% ylabel('Weights');
% legendStr = num2str((1 : M).');
% leg = legend(legendStr, 'Location', 'best');
% title(leg, 'm');

% w{4}(:, 61).'
% w{4}(:, 31).'

% figure; plot(w{3}(2 : end, 1));

% figure; plot(T{4});
if kayit
    figureName = 'analysis_decentralized_different_snr_iid_3';
    % savefig(fig, ['\\Users\anilarslan\Desktop\MSRS\figureGenerator\figures\' figureName '.fig']);
    % saveas(fig, ['\\Users\anilarslan\Desktop\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
    savefig(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end


%%%