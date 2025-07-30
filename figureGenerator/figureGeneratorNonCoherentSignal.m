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
kayit = 0;
set(groot, "defaultFigurePosition", [680 458 560 420]);

globalPFA = 1e-6;
algorithms = ["SLC", "WSLC"];
snrs = -5 : .02 : 15;
M = [9 5];

fig1 = figure;
for m = M
    algorithmID = 0;
    for algorithm = algorithms
        algorithmID = algorithmID + 1;
        dtc = detector( ...
            "globalPFA", globalPFA, ...
            "numberOfSensors", m, ...
            "SNR_input_dB", snrs, ...
            "localPFA", 1 ...
            );
        rng(1975); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 20);
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
        switch m
            case 9
                plot(snrs, pd, 'b', 'LineWidth', 2, 'LineStyle', lineStyle); hold on;
            case 5
                plot(snrs, pd, 'r', 'LineWidth', 2, 'LineStyle', lineStyle); hold on;
        end
    end
end
grid on; grid minor;
ylim([0, 1]);
xlabel('Mean of Average SNRs per Receiving Node (dB)');
ylabel('Global Probability of Detection');
legend(["unweighted SLC - M = 9", "weighted SLC - M = 9", "unweighted SLC - M = 5", "weighted SLC - M = 5"], 'Location', 'best');

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
pfaLocal = logspace(0, -6, 4001);
algorithms = ["BC", "SLC"];
M = 1 : 5;

fig1 = figure;
algorithmID = 0;
for algorithm = algorithms
    algorithmID = algorithmID + 1;
    for sensorID = 1 : length(M)
        dtc = detector( ...
            "globalPFA", globalPFA, ...
            "numberOfSensors", M(sensorID), ...
            "SNR_input_dB", 10, ...
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
        figure(fig1); semilogx(pfaLocal, pd, 'LineWidth', 2, 'LineStyle', lineStyle, 'color', colors(sensorID)); hold on;
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
pfaLocal = logspace(0, -6, 1001);
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
        dtc.setalgorithm( ...
            "binaryDetectionRule", "notSpecified", ...
            "binaryDetectionConstraint", "fixedGlobal|LocalPFA", ...
            "signalPhaseModel", "decorrelatedUniform", ...
            "signalAmplitudeModel", "decorrelatedExponential", ...
            "globalFusionRule", algorithm ...
            );
        pd = squeeze(dtc.globalPD);
        semilogx(pfaLocal, pd, 'LineWidth', 2); hold on;
    end
end
grid on; grid minor;
ylim([0, 1]);
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of Detection');
legend(["SLC", "BC", "weighted SLC", "weighted BC"], 'Location', 'best');

if kayit
    figureName = 'analysis_decentralized_different_snr_iid';
    savefig(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end