%%% Anıl ARSLAN
% 2303980

%%% LRT
% EGC: Equal Gain Combining (LRT under equal SNR)                               Implemented
% MRC: Maximal Ratio Combining (LRT under different SNR)                        Implemented

%%% GLRT
% LDC: Log-Divergence Combining (GLRT w/ joint MLE)                             Implemented
% GLDC: Generalized Log-Divergence Combining (GLRT w/ independent MLE)          No closed form expression

%%% Binary
% BC: Binary Combining (LRT under equal SNR)                                    Implemented
% CVBC: Chair-Varshney Binary Combining (LRT under different SNR)               Implemented

%%%
% NEGC: Normalized Equal Gain Combining                                         Implemented
% MRCbiased: Maximal Ratio Combining without unbias term                        Implemented
% SC: Selective Combining (Partial Selective Combining K = 1)                   Implemented
% PSC: Partial Selective Combining (Combine largest K samples)


% implemented = ROC curve under different SNR


% clc;
% clear;
% close all;
addpath(genpath([pwd '/library']));

averageSNR_dB = 10; % dB
snr = flipud(linspace(45, 0, 10).');
snrs = cell(1, 10);
for i = 1 : 10
    snrs{i} = 10*log10(i*10.^(0.1*snr(1 : i))./sum(10.^(0.1*snr(1 : i)))) + averageSNR_dB;
end

detector = detectorNoncoherent( ...
    "globalPFA", 1e-6, ...
    "numberOfSensors", 1 : 10, ...
    "numberOfSensors", 10 , ...
    "SNR_input_dB", 6, ...
    ..."SNR_input_dB", 10, ...
    ..."SNR_input_dB", linspace(0, 20, 20), ...
    ..."localPFA", [logspace(0, -8, 9), 0.5] ...
    ..."localPFA", logspace(-1, -8, 8) ...
    "localPFA", logspace(0, -8, 9) ...
    ..."localPFA", 1 ...
    ..."localPFA", .1 ...
    );
% detector.setbudgetconstraint("constraint", ["dataRate", "transmittedPower"]);
detector.setbudgetconstraint("constraint", "dataRate");
detector.setalgorithm( ...
    "globalFusionRule", "BC", ...
    "numberOfBinaryDetectionsRule", "fixedGlobalPFA");
% detector.setuserdefinednumberofdetections("userDefinedNumberOfDetections", 5);
detector.visualize( ...
    "x_axis", ["numberOfSensors", "SNR"], ...
    "x_axis", "numberOfSensors", ...
    "y_axis", "globalPD", "dataType", "analytical");
detector.visualize( ...
    "x_axis", ["numberOfSensors", "SNR"], ...
    "x_axis", "numberOfSensors", ...
    "x_axis", ["localPFA", "globalThreshold"], ...
    "y_axis", "globalPFA", "dataType", "analytical");



%% PD simulator
detector.setmontecarlo("numberOfTrials", 1e5);
detector.simulate("printStatus", 1, "simulationData", "globalPD");
detector.visualize( ...
    "x_axis", ["globalPFA", "numberOfSensors", "SNR", "globalThreshold"], ...
    "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);
detector.visualize( ...
    "x_axis", ["localPFA", "localThreshold", "localPD"], ...
    "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);




%% High Monte Carlo PFA simulator
clc;
clear;

detector = detectorNoncoherent( ...
    "globalPFA", 1e-5, ...
    "numberOfSensors", 1 : 10, ...
    "SNR_input_dB", 14, ...
    "SNR_input_dB", flipud(linspace(14, 0, 10).'), ... % meaningful for SNR dependent weighted fusion rules e.g. MRC
    ..."localPFA", 1 ...
    "localPFA", logspace(0, -8, 5) ...
    ..."localPFA", logspace(0, -2, 9) ...
    );
% detector.setconstraint;
detector.setmontecarlo("numberOfTrials", 2e7);
detector.setalgorithm("globalFusionRule", "MRCbiased");
detector.simulate("printStatus", 1, "simulationData", "globalPFA");

% Visualization
% close all;
detector.visualize( ...
    "x_axis", "numberOfSensors", ...
    "y_axis", "globalPFA", "dataType", ["analytical", "empirical"]);
% detector.visualize( ...
%     "x_axis", ["globalPFA", "numberOfSensors", "SNR", "globalThreshold", "localPFA", "localThreshold", "localPD"], ...
%     "y_axis", "globalPFA", "dataType", ["analytical", "empirical"]);





%%

%%%% TO DO
    % mean/variance formulas, E[M] both simulation and formula
    % average output SNR
    % what is mean p_d, conditioned on given SNR?
    % binary integration: 1 bit, others: float; 32/16/8 bit (Power = nonnegative real number, so also unsigned)
    % independent yapma (swerling 0)
    % binary integration da local PFA fixleme global PFA ayarlamak icin kullan

    % tum integer thresholdlar icin unique local PFA cikacaktir.
    % ayni global thresholdda dusuk local PFA dusuk global PFA demek
    % yuksek global threshold --> dusuk local PFA yani ters oranti ve monotonic
    % mesela m out of 5in sabit local PFAde ki global PFAi m out of 4 deki duruma gore daha yuksek olacaktir.
    % global PFAi dusurmek icin m out of 5 de daha dusuk local PFA kullanmak gerekir
    % yani yuksek global threshold (5) daha dusuk local PFA ister
    % sabit global PFA icin her integer threshold degerinde unique local PFA vardir.

    % average snr in sabit oldugu random snr uret
    % discrete threshold u dB olarak cizme


%%%% Questions
    % optimum rule under local thresholding? GLRT LRT
    % is EGC optimum?
    % SNR/sqrt(M) [since there will M^2 samples]
    % Kullback-Leibler divergence x - log(x) - 1
    % censoring










%% Calisma, order statistics
M = 11;
K = 4;

F = @(t) 1 - sum(arrayfun(@(j) (-1).^j*nchoosek(M - K, j)*gammainc(t, K + j, 'lower'), 0 : (M - K)));

x = linspace(0, 30, 1000);

for i = 1 : length(x)
    % c(i) = F(x(i))*K*nchoosek(M, K);
    c(i) = F(x(i));
end

figure; plot(x, c);
