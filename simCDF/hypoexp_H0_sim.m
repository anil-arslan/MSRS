clc; clear; close all;
kayit = 0;
set(groot, "defaultFigurePosition", [680 458 560 420]);

dtc = detector( ...
    "globalPFA", 1e-6, ...
    "numberOfSensors", 9, ...
    "SNR_input_dB", 8, ...
    "localPFA", 1e-3 ...
    );
dtc.setalgorithm( ...
    "signalAmplitudeModel", "decorrelatedExponential", ...
    "signalPhaseModel", "decorrelatedUniform", ...
    "globalFusionRule", "WSLC" ...
    );
rng(1975); dtc.randomSNRwithFixedAverage("rangeSNR_dB", 20);
ccdf = dtc.totalCDF(9, 9, dtc.localThreshold{1}, dtc.localPFA{1}, dtc.fusionWeights{1}, dtc.SNR_input_dB{1});
T = logspace(0, 2, 101);
pfa = zeros(1, length(T));
for i = 1 : length(T)
    pfa(i) = ccdf(T(i));
end
figure; loglog(T, pfa);
grid off; grid on; grid minor;
xlabel('Global Threshold');
ylabel('Global Probability of False Alarm');




if kayit
    figureName = 'square_law_combining_weighted';
    savefig(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.fig']);
    saveas(fig, ['C:\GitRepo\MSRS\figureGenerator\figures\' figureName '.eps'], 'epsc');
end