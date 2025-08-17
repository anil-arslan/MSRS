%% Centralized Equal SNR, SW 1 PD vs SNR
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
algorithms = ["SLC", "EGC"];
snrs = -10 : 1 : 10;
snr = 8;
M = [9 90];
M = 1 : 30;

fig = figure;
% for m = M
    algorithmID = 0;
    for algorithm = algorithms
        algorithmID = algorithmID + 1;
        dtc = detector( ...
            "globalPFA", globalPFA, ...
            "numberOfSensors", 1 : 30, ...
            "SNR_input_dB", snr./(1 : 30), ...
            "localPFA", 1 ...
            );
        % rng(1975); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 20);
        dtc.setalgorithm( ...
            "signalPhaseModel", "correlatedUniform", ...
            "signalAmplitudeModel", "correlatedExponential", ...
            "globalFusionRule", algorithm ...
            );
        % switch algorithm
        %     case "WSLC"
                % dtc.setmontecarlo("numberOfTrials", 1e8); % 1e8 yi gecme
                dtc.setmontecarlo("numberOfTrials", 1e6); % 1e8 yi gecme
                dtc.simulate("printStatus", 1, "simulationData", "globalThreshold", "statistics", 1);
        % end
        % switch algorithm
        %     case "WSLC"
                dtc.setmontecarlo("numberOfTrials", 1e4); % 1e6 yi gecme
                dtc.simulate("printStatus", 1, "simulationData", "globalPD", "statistics", 1);
                pd = squeeze(dtc.globalPDsimulation);
        %     otherwise
        %         pd = squeeze(dtc.globalPD);
        % end
        switch algorithmID
            case 1
                lineStyle = '-.';
            case 2
                lineStyle = '-';
        end
        % switch m
        %     case 9
                plot(M, diag(pd), 'LineWidth', 2, 'LineStyle', lineStyle); hold on;
        %     case 90
        %         plot(snrs, pd, 'r', 'LineWidth', 2, 'LineStyle', lineStyle); hold on;
        % end
    end
% end
grid on; grid minor;
ylim([0, 1]);
xlabel('SNR (dB)');
ylabel('Global Probability of Detection');
legend(["unweighted SLC - M = 9", "weighted SLC - M = 9", "unweighted SLC - M = 90", "weighted SLC - M = 90"], 'Location', 'best');

figure; plot(dtc.SNR_input_dB{1}(:, [1 end]));
xlabel('Receiving Node ID');
ylabel('Average SNR (dB)');
yline(snrs(1), 'k--');
yline(snrs(end), 'k--');
grid on; grid minor;
leg = legend(num2str(snrs([1 end]).'), 'Location', 'best');
title(leg, 'mean whole SNR');

if kayit
    figureName = 'analysis_equal_snr_iid';
    savefig(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end

%% Centralized Different SNR, SW 1 PD vs SNR
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
algorithms = ["SLC", "EGC", "WSLC", "MRC"];
snrs = -10 : 1 : 10;
snr = 8;
M = [9 90];
M = 1 : 30;

fig = figure;
% for m = M
    algorithmID = 0;
    for algorithm = algorithms
        algorithmID = algorithmID + 1;
        dtc = detector( ...
            "globalPFA", globalPFA, ...
            "numberOfSensors", 9, ...
            "SNR_input_dB", snrs, ...
            "localPFA", 1 ...
            );
        % dtc.setbudgetconstraint("constraint", "transmittedPower");
        rng(1975); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 20);
        dtc.setalgorithm( ...
            "signalPhaseModel", "correlatedUniform", ...
            "signalAmplitudeModel", "correlatedExponential", ...
            "globalFusionRule", algorithm ...
            );
        % switch algorithm
        %     case "WSLC"
                % dtc.setmontecarlo("numberOfTrials", 1e8); % 1e8 yi gecme
                dtc.setmontecarlo("numberOfTrials", 1e6); % 1e8 yi gecme
                dtc.simulate("printStatus", 1, "simulationData", "globalThreshold", "statistics", 1);
        % end
        % switch algorithm
        %     case "WSLC"
                dtc.setmontecarlo("numberOfTrials", 1e4); % 1e6 yi gecme
                dtc.simulate("printStatus", 1, "simulationData", "globalPD", "statistics", 1);
                pd = squeeze(dtc.globalPDsimulation);
        %     otherwise
        %         pd = squeeze(dtc.globalPD);
        % end
        % switch algorithmID
        %     case 1
        %         lineStyle = '-.';
        %     case 2
        %         lineStyle = '-';
        % end
        % switch m
        %     case 9
                plot(snrs, pd, 'LineWidth', 2); hold on; drawnow;
        %     case 90
        %         plot(snrs, pd, 'r', 'LineWidth', 2, 'LineStyle', lineStyle); hold on;
        % end
    end
% end
grid on; grid minor;
ylim([0, 1]);
xlabel('SNR (dB)');
ylabel('Global Probability of Detection');
legend(["unweighted SLC - M = 9", "weighted SLC - M = 9", "unweighted SLC - M = 90", "weighted SLC - M = 90"], 'Location', 'best');

figure; plot(dtc.SNR_input_dB{1}(:, [1 end]));
xlabel('Receiving Node ID');
ylabel('Average SNR (dB)');
yline(snrs(1), 'k--');
yline(snrs(end), 'k--');
grid on; grid minor;
leg = legend(num2str(snrs([1 end]).'), 'Location', 'best');
title(leg, 'mean whole SNR');

if kayit
    figureName = 'analysis_equal_snr_iid';
    savefig(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end








%% Decentralized Equal SNR, SW 1 PD vs PFAlocal
clc; clear; % close all;
kayit = 0;
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
pfaLocal = logspace(0, -6, 1001);
algorithms = ["BC", "SLC"];
M = [1 2 4 12 42];

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
            "signalPhaseModel", "correlatedUniform", ...
            "signalAmplitudeModel", "correlatedExponential", ...
            "globalFusionRule", algorithm ...
            );
        switch algorithmID
            case 1
                lineStyle = '-';
            case 2
                lineStyle = '-.';
        end
        switch algorithm
            case "EGC"
                dtc.setmontecarlo("numberOfTrials", 1e6); % 1e8 yi gecme
                dtc.simulate("printStatus", 1, "simulationData", "globalThreshold", "statistics", 1);
                T = squeeze(dtc.globalThresholdSimulation);
        end
        dtc.setmontecarlo("numberOfTrials", 1e5); % 1e6 yi gecme
        dtc.simulate("printStatus", 1, "simulationData", "globalPD", "statistics", 1);
    
        pd = squeeze(dtc.globalPDsimulation);
        
        figure(fig1); semilogx(pfaLocal, pd, 'LineWidth', 2, 'LineStyle', lineStyle, 'Color', colors(sensorID)); hold on;
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
    figureName = 'analysis_decentralized_equal_snr_coh';
    savefig(fig1, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig1, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end


%% Decentralized Different SNR, SW 1 PD vs PFAlocal
clc; clear; % close all;
kayit = 1;
set(groot, "defaultFigurePosition", [680 458 560 420]);

globalPFA = 1e-6;
pfaLocal = logspace(0, -6, 101);
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
            "signalAmplitudeModel", "correlatedExponential", ...
            "globalFusionRule", algorithm ...
            );
        % switch algorithm
        %     case "WSLC"
                % dtc.setmontecarlo("numberOfTrials", 1e8); % 1e8 yi gecme
        %         dtc.setmontecarlo("numberOfTrials", 1e7); % 1e8 yi gecme
        %         dtc.simulate("printStatus", 1, "simulationData", "globalThreshold", "statistics", 1);
        % % end
        % % switch algorithm
        % %     case "WSLC"
        %         dtc.setmontecarlo("numberOfTrials", 1e4); % 1e6 yi gecme
        %         dtc.simulate("printStatus", 1, "simulationData", "globalPD", "statistics", 1);
        %         pd = squeeze(dtc.globalPDsimulation);
        %     otherwise

        dtc.setmontecarlo("numberOfTrials", 1e5); % 1e6 yi gecme
        dtc.simulate("printStatus", 1, "simulationData", "globalPD", "statistics", 1);
    
        pd = squeeze(dtc.globalPDsimulation);
        semilogx(pfaLocal, pd, 'LineWidth', 2); hold on;
    end
end
grid on; grid minor;
ylim([0, 1]);
xlabel('Local Probability of False Alarm');
ylabel('Global Probability of Detection');
% legend(["unweighted SLC - M = 9", "weighted SLC - M = 9", "unweighted SLC - M = 90", "weighted SLC - M = 90"], 'Location', 'best');
legend(["SLC", "BC", "weighted SLC", "weighted BC"], 'Location', 'best');

if kayit
    figureName = 'analysis_decentralized_different_snr_coh';
    savefig(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end











%% Binary Fusion Equal SNR, SW 1 and SW 2 PD vs PFAlocal
clc; clear; % close all;
kayit = 0;
set(groot, "defaultFigurePosition", [680 458 560 420]);
colors = rgb2hex([ ...
    0.0000    0.4470    0.7410 % blue
    0.8500    0.3250    0.0980 % orange
    0.9290    0.6940    0.1250 % yellow
    0.4940    0.1840    0.5560 % magenta
    0.4660    0.6740    0.1880 % green
    0.3010    0.7450    0.9330 % cyan
    0.6350    0.0780    0.1840]);

globalPFA = 1e-4;
pfaLocal = logspace(0, -8, 101);
signalModels = ["correlatedUniform", "correlatedExponential"; "decorrelatedUniform", "decorrelatedExponential"];
M = [2 7 12];

fig1 = figure;
modelID = 0;
for signalModel = signalModels.'
    modelID = modelID + 1;
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
            "signalPhaseModel", signalModel(1), ...
            "signalAmplitudeModel", signalModel(2), ...
            "globalFusionRule", "BC" ...
            );
        switch modelID
            case 1
                lineStyle = '-';
            case 2
                lineStyle = '-.';
        end
        dtc.setmontecarlo("numberOfTrials", 1e5); % 1e6 yi gecme
        dtc.simulate("printStatus", 1, "simulationData", "globalPD", "statistics", 1);
    
        pd = squeeze(dtc.globalPDsimulation);
        
        figure(fig1); semilogx(pfaLocal, pd, 'LineWidth', 2, 'LineStyle', lineStyle, 'Color', colors(sensorID)); hold on;
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
    figureName = 'analysis_bc_equal_snr';
    savefig(fig1, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig1, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end