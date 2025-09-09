% Anıl Arslan 2303980
addpath(genpath('C:\GitRepo\MSRS/library'));
addpath(genpath('C:\GitRepo\MSRS/figureGenerator'));

%% Binary combining discrete threshold effect
clc; clear; close all;
kayit = 0;
set(groot, "defaultFigurePosition", [680 458 560 420]);

algorithms = ["fixedLocalPFA", "fixedGlobal|LocalPFA"];
signalModels = ["decorrelatedExponential", "deterministic"];
signalModels = "decorrelatedExponential";
pfaLocal = logspace(0, -8, 1001);

dtc = detector( ...
    "globalPFA", 1e-6, ...
    "numberOfSensors", 9, ...
    "SNR_input_dB", 8, ...
    "localPFA", pfaLocal ...
    );
dtc.setalgorithm( ...
    "binaryDetectionRule", "notSpecified", ...
    "signalPhaseModel", "decorrelatedUniform", ...
    "globalFusionRule", "BC" ...
    );

fig1 = figure;
fig2 = figure;
modelID = 0;
for signalModel = signalModels
    modelID = modelID + 1;
    algorithmID = 0;
    for algorithm = algorithms
        algorithmID = algorithmID + 1;
        dtc.setalgorithm("binaryDetectionConstraint", algorithm, "signalAmplitudeModel", signalModel);
        pd = squeeze(dtc.globalPD);
        pfa = squeeze(dtc.globalPFAanalytical);
        T = squeeze(dtc.globalThreshold);
        q = squeeze(dtc.globalRandomizationProbability);
        switch algorithmID
            case 1
                color = 'r';
                lineStyle = '-.';
            case 2
                color = 'b';
                lineStyle = '-';
        end
        if modelID == 1
            figure(fig1); plot(pfaLocal, pd, 'LineWidth', 2, 'Color', color, 'LineStyle', lineStyle); hold on;
            figure(fig2); loglog(pfaLocal, pfa, 'LineWidth', 2, 'Color', color); hold on;
        else
            figure(fig1); plot(pfaLocal, pd, 'LineWidth', 2, 'Color', color, 'LineStyle', lineStyle); hold on;
        end
    end
end
figure(fig2); yyaxis right; semilogx(pfaLocal, T - q - .5, 'LineWidth', 2, 'LineStyle', '-.', 'Color', 'm'); hold on;
ax = gca;
ax.YAxis(2).Color = 'm';

figure(fig1);
grid on; grid minor;
ylim([0, 1]);
xlim([0, 0.3]);
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of Detection');
legendStr = ["fixed thresholding", "randomized thresholding"];
legend(legendStr, 'Location', 'best');

figure(fig2);
grid on; grid minor;
xlabel('Local Probability of False Alarm');
yyaxis left; ylabel('Global Probability of False Alarm'); ylim([1e-14 1e-5])
yyaxis right; ylabel('Effective Global Threshold'); ylim([-0.5 9.5]);
legendStr = ["fixed thresholding", "randomized thresholding", "effective threshold"];
legend(legendStr, 'Location', 'best');

if kayit
    figureName = 'binary_combining_pd_threshold_effect_2';
    savefig(fig1, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig1, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
    % figureName = 'binary_combining_pfa_threshold_effect';
    % savefig(fig2, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    % saveas(fig2, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end

%% Unweighted binary combining
clc; clear; close all;
kayit = 1;
set(groot, "defaultFigurePosition", [680 458 560 420]);
dtc = detector( ...
    "globalPFA", 1e-6, ...
    "numberOfSensors", [1 2 3 4 5 7 9], ...
    "SNR_input_dB", 8, ...
    "localPFA", logspace(0, -8, 10001) ...
    );
dtc.setalgorithm( ...
    "binaryDetectionRule", "notSpecified", ...
    "binaryDetectionConstraint", "fixedGlobal|LocalPFA", ...
    "signalAmplitudeModel", "decorrelatedExponential", ...
    "signalPhaseModel", "decorrelatedUniform", ...
    "globalFusionRule", "BC" ...
    );
% dtc.setbudgetconstraint("constraint", "transmittedPower");
fig = dtc.visualize( ...
    "x_axis", "localPFA", ...
    "y_axis", "globalPD", "dataType", "analytical");
% q = squeeze(dtc.globalRandomizationProbability);
% figure; semilogx(logspace(0, -6, 1001), q(end - 1, :));

if kayit
    figureName = 'binary_combining';
    savefig(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end

%% Weighted binary combining
clc; clear; close all;
kayit = 0;
set(groot, "defaultFigurePosition", [680 458 560 420]);

pfaLocal = logspace(0, -8, 81);
algorithms = ["BC", "CVBC"];
signalModels = ["decorrelatedExponential", "deterministic"];
signalModels = "decorrelatedExponential";

fig = figure;
modelID = 0;
for signalModel = signalModels
    modelID = modelID + 1;
    algorithmID = 3;
    for algorithm = algorithms
        algorithmID = algorithmID - 1;
        dtc = detector( ...
            "globalPFA", 1e-6, ...
            "numberOfSensors", 9, ...
            "SNR_input_dB", 8, ...
            "localPFA", pfaLocal ...
            );
        dtc.setalgorithm( ...
            "binaryDetectionRule", "notSpecified", ...
            "binaryDetectionConstraint", "fixedGlobal|LocalPFA", ...
            "signalAmplitudeModel", signalModel, ...
            "signalPhaseModel", "decorrelatedUniform", ...
            "globalFusionRule", algorithm ...
            );
        % rng(1976); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 10);
        rng(1975); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 20);
        snrr = dtc.SNR_input_dB{1};
        pd = squeeze(dtc.globalPD);
        switch algorithmID
            case 1
                lineStyle = '-.';
            case 2
                lineStyle = '-';
        end
        figure(fig); semilogx(pfaLocal, pd, 'LineWidth', 2, 'LineStyle', lineStyle); hold on;
        w{algorithmID} = squeeze(dtc.fusionWeights{1});
    end
end
figure(fig);
grid on; grid minor;
ylim([0, 1]);
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of Detection');
leg = legend(["BC", "WBC"], 'Location', 'best');

if kayit
    figureName = 'binary_combining_weighted_2';
    savefig(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end

%% Unweighted binary combining ROC
% clc; clear; close all;
% kayit = 0;
% set(groot, "defaultFigurePosition", [680 458 560 420]);
% dtc = detector( ...
%     "globalPFA", logspace(0, -6, 1001), ...
%     "numberOfSensors", 9, ...
%     "SNR_input_dB", linspace(2, 8, 7), ...
%     "localPFA", 1e-3 ...
%     );
% dtc.setalgorithm( ...
%     "binaryDetectionRule", "notSpecified", ...
%     "binaryDetectionConstraint", "fixedGlobal|LocalPFA", ...
%     "signalAmplitudeModel", "decorrelatedExponential", ...
%     "signalPhaseModel", "decorrelatedUniform", ...
%     "globalFusionRule", "BC" ...
%     );
% % dtc.setbudgetconstraint("constraint", "transmittedPower");
% fig = dtc.visualize( ...
%     "x_axis", "globalPFA", ...
%     "y_axis", "globalPD", "dataType", "analytical");
% 
% if kayit
%     figureName = 'binary_combining_roc';
%     savefig(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
%     saveas(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
% end

%% Unweighted square Law Combining
clc; clear; close all;
kayit = 1;
set(groot, "defaultFigurePosition", [680 458 560 420]);
dtc = detector( ...
    "globalPFA", 1e-6, ...
    "numberOfSensors", [1 2 3 4 5 7 9], ...
    "SNR_input_dB", 8, ...
    "localPFA", logspace(0, -8, 10001) ...
    );
dtc.setalgorithm( ...
    "signalAmplitudeModel", "decorrelatedExponential", ...
    "signalPhaseModel", "decorrelatedUniform", ...
    "globalFusionRule", "SLC" ...
    );
fig = dtc.visualize( ...
    "x_axis", "localPFA", ...
    "y_axis", "globalPD", "dataType", "analytical");

if kayit
    figureName = 'square_law_combining';
    savefig(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end

%% Weighted square Law Combining
clc; clear; close all;
kayit = 1;
set(groot, "defaultFigurePosition", [680 458 560 420]);
colors = rgb2hex([ ...
    0.0000    0.4470    0.7410 % blue
    0.8500    0.3250    0.0980 % orange
    0.4940    0.1840    0.5560 % magenta
    0.9290    0.6940    0.1250 % yellow
    0.4660    0.6740    0.1880 % green
    0.3010    0.7450    0.9330 % cyan
    0.6350    0.0780    0.1840]);

pfaLocal = logspace(0, -8, 201);
algorithms = ["SLC", "WSLC"];
snrs = [2 8];

fig = figure;
snrID = 0;
for snr = snrs
    snrID = snrID + 1;
    algorithmID = 3;
    for algorithm = algorithms
        algorithmID = algorithmID - 1;
        dtc = detector( ...
            "globalPFA", 1e-6, ...
            "numberOfSensors", 9, ...
            "SNR_input_dB", snr, ...
            "localPFA", pfaLocal ...
            );
        dtc.setalgorithm( ...
            "signalAmplitudeModel", "decorrelatedExponential", ...
            "signalPhaseModel", "decorrelatedUniform", ...
            "globalFusionRule", algorithm ...
            );
        rng(1975); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 20);
        % rng(1976); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 10);
        snrr = dtc.SNR_input_dB{1};
        pd = squeeze(dtc.globalPD);
        switch algorithmID
            case 1
                lineStyle = '-.';
                figure(fig); semilogx(pfaLocal, pd, 'LineWidth', 2, 'Color', colors(snrID), 'LineStyle', lineStyle); hold on;
            case 2
                lineStyle = '-';
                figure(fig); semilogx(pfaLocal, pd, 'LineWidth', 2, 'Color', colors(snrID + 2), 'LineStyle', lineStyle); hold on;
        end
    end
end
figure(fig);
grid on; grid minor;
ylim([0, 1]);
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of Detection');
leg = legend(["SLC - 2 dB", "WSLC - 2 dB", "SLC - 8 dB", "WSLC - 8 dB"], 'Location', 'best');

if kayit
    figureName = 'square_law_combining_weighted';
    savefig(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end