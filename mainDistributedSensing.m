%%% Anıl ARSLAN 2303980
clc; % clear; % close all;
addpath(genpath([pwd '/library']));

detector = detectorNoncoherent( ...
    "globalPFA", 1e-6, ...
    "numberOfSensors", 1 : 10, ...
    "SNR_input_dB", 0 : 1 : 10, ...
    "localPFA", 1e-1);
detector.setmontecarlo("numberOfTrials", 1e4);
detector.simulate("printStatus", 1, "simulationData", ["globalPD", "globalPFA"]);

%%
detector.visualize("x_axis", "globalPFA", "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);
detector.visualize("x_axis", "numberOfSensors", "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);
detector.visualize("x_axis", "SNR", "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);
detector.visualize("x_axis", "globalThreshold", "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);

% detector.visualize("x_axis", "globalPFA", "y_axis", "globalPFA", "dataType", ["analytical", "empirical"]);
% detector.visualize("x_axis", "numberOfSensors", "y_axis", "globalPFA", "dataType", ["analytical", "empirical"]);
% detector.visualize("x_axis", "SNR", "y_axis", "globalPFA", "dataType", ["analytical", "empirical"]);


%%%% TO DO
% weighting
% different SNR

% 0.75/M
% SNR/M ile SNR vs Pd
% 1/M
% Binary Integration, 1, 0 1 bit, biz 32 bit iz ??