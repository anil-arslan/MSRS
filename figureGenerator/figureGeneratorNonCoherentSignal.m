% Anıl Arslan 2303980
addpath(genpath('C:\GitRepo\MSRS/library'));
addpath(genpath('C:\GitRepo\MSRS/figureGenerator'));

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
kayit = 1;
set(groot, "defaultFigurePosition", [680 458 560 420]);

globalPFA = 1e-6;
algorithms = ["SLC", "WSLC"];
snrs = -5 : .02 : 20;
M = 9;

fig1 = figure;
algorithmID = 3;
for algorithm = algorithms
    algorithmID = algorithmID - 1;
    dtc = detector( ...
        "globalPFA", globalPFA, ...
        "numberOfSensors", M, ...
        "SNR_input_dB", snrs, ...
        "localPFA", 1 ...
        );
    % rng(1975); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 20);
    rng(1976); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 10);
    dtc.setalgorithm( ...
        "signalPhaseModel", "decorrelatedUniform", ...
        "signalAmplitudeModel", "decorrelatedExponential", ...
        "globalFusionRule", algorithm ...
        );
    pd = squeeze(dtc.globalPD);
    switch algorithmID
        case 1
            lineStyle = '-.';
        case 2
            lineStyle = '-';
    end
    plot(snrs, pd, 'LineWidth', 2, 'LineStyle', lineStyle); hold on;
end
grid on; grid minor;
ylim([0, 1]);
xlabel('Mean of Average SNRs per Receiving Node (dB)');
ylabel('Global Probability of Detection');
legend(["SLC", "WSLC"], 'Location', 'best');

if kayit
    figureName = 'analysis_centralized_different_snr_iid_2';
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
algorithms = ["BC", "SLC"];
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
    figureName = 'analysis_decentralized_equal_snr_iid';
    saveas(fig1, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
    savefig(fig1, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
end

%% Decentralized Different SNR, SW 2 PD vs PFAlocal
clc; clear; % close all;
kayit = 1;
set(groot, "defaultFigurePosition", [680 458 560 420]);

globalPFA = 1e-6;
pfaLocal = logspace(0, -8, 201);
algorithms = ["SLC", "BC", "WSLC", "CVBC"];
M = 9;
snr = 8;

fig = figure;
for m = M
    algorithmID = 0;
    for algorithm = algorithms
        algorithmID = algorithmID + 1;
        dtc = detector( ...
            "globalPFA", globalPFA, ...
            "numberOfSensors", m, ...
            "SNR_input_dB", snr, ...
            "localPFA", pfaLocal ...
            );
        rng(1975); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 20);
        rng(1976); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 10);
        dtc.setalgorithm( ...
            "binaryDetectionRule", "notSpecified", ...
            "binaryDetectionConstraint", "fixedGlobal|LocalPFA", ...
            "signalPhaseModel", "decorrelatedUniform", ...
            "signalAmplitudeModel", "decorrelatedExponential", ...
            "globalFusionRule", algorithm ...
            );
        pd = squeeze(dtc.globalPD);
        T{algorithmID} = squeeze(dtc.globalThreshold);
        semilogx(pfaLocal, pd, 'LineWidth', 2); hold on;
        w{algorithmID} = squeeze(dtc.fusionWeights{1});
    end
end
grid on; grid minor;
ylim([0, 1]);
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of Detection');
legend(["SLC", "BC", "WSLC", "WBC"], 'Location', 'best');

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
    figureName = 'analysis_decentralized_different_snr_iid_2';
    savefig(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end