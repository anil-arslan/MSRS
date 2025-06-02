%%% Anıl ARSLAN
% 2303980

%%% LRT
% EGC: Equal Gain Combining (LRT under equal SNR)                               Implemented
% MRC: Maximal Ratio Combining (LRT under different SNR)                        Implemented
% MRCbiased: 

%%% GLRT
% LDC: Log-Divergence Combining (GLRT w/ joint MLE)                             Implemented
% GLDC: Generalized Log-Divergence Combining (GLRT w/ independent MLE)

%%% Binary
% BC: Binary Combining (LRT under equal SNR)                                    Implemented
% CVBC: Chair-Varshney Binary Combining (LRT under different SNR)               Implemented

%%%
% NEGC: Normalized Equal Gain Combining                                         Implemented
% SC: Selective Combining (Partial Selective Combining K = 1)                   Implemented
% PSC: Partial Selective Combining (Combine largest K samples)


% implemented = ROC curve under different SNR


% clc;
% clear;
% close all;
addpath(genpath([pwd '/library']));

detector = detectorNoncoherent( ...
    "globalPFA", 1e-6, ...
    "numberOfSensors", 1 : 15, ...
    "SNR_input_dB", flipud(linspace(10, 0, 15).'), ...
    "SNR_input_dB", 6, ...
    ..."localPFA", [logspace(0, -8, 9), 0.5] ...
    ..."localPFA", logspace(-1, -8, 8) ...
    "localPFA", logspace(0, -8, 9) ...
    ..."localPFA", 1 ...
    ..."localPFA", .1 ...
    );
% detector.setbudgetconstraint("constraint", ["dataRate", "transmittedPower"]);
detector.setalgorithm( ...
    "globalFusionRule", "CVBC", ...
    "numberOfBinaryDetectionsRule", "fixedGlobalPFA");
% detector.setuserdefinednumberofdetections("userDefinedNumberOfDetections", 5);
detector.visualize( ...
    "x_axis", ["numberOfSensors", "SNR"], ...
    "x_axis", "numberOfSensors", ...
    "y_axis", "globalPD", "dataType", "analytical");
detector.visualize( ...
    "x_axis", ["numberOfSensors", "SNR"], ...
    "x_axis", "numberOfSensors", ...
    "y_axis", "globalPFA", "dataType", "analytical");

%% Simulation
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
    "SNR_input_dB", flipud(linspace(14, 0, 10).'), ... % meaningful for MRC
    ..."localPFA", 1 ...
    "localPFA", logspace(0, -8, 5) ...
    ..."localPFA", logspace(0, -2, 9) ...
    );
% detector.setconstraint;
detector.setmontecarlo("numberOfTrials", 1e7);
detector.setalgorithm("globalFusionRule", "CVBC");
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
    % mean/variance formulas
    % average output SNR
    % what is mean p_d, conditioned on given SNR?
    % binary integration: 1 bit, others: float; 32/16/8 bit (Power = nonnegative real number, so also unsigned)


%%%% Questions
    % optimum rule under local thresholding? GLRT LRT
    % is EGC optimum?
    % SNR/sqrt(M) [since there will M^2 samples]
    % Kullback-Leibler divergence x - log(x) - 1










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
