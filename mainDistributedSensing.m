%%% Anıl ARSLAN 2303980
clc; % clear; % close all;
addpath(genpath([pwd '/library']));

%%% lejant yanlis
detector = detectorNoncoherent( ...
    "globalPFA", [1e-6 1e-4], ...
    "numberOfSensors", 1 : 20, ...
    "SNR_input_dB", 15, ...
    "localPFA", 1e-1);
detector.setconstraint;
detector.visualize("x_axis", "localPFA", "y_axis", "globalPD", "dataType", "analytical");



%%% cizim cakildi
detector = detectorNoncoherent( ...
    "globalPFA", [1e-6 1e-4], ...
    "numberOfSensors", 1 : 20, ...
    "SNR_input_dB", [15 20], ...
    "localPFA", 1e-1);
% detector.visualize("x_axis", "localPFA", "y_axis", "globalPD", "dataType", "analytical");
detector.setconstraint;
detector.visualize("x_axis", "localPFA", "y_axis", "globalPD", "dataType", "analytical");



%%
% Simulation
detector.setmontecarlo("numberOfTrials", 1e4);
detector.simulate("printStatus", 1, "simulationData", ["globalPD", "globalPFA"]);

% detector.setmontecarlo("numberOfTrials", 1e5);
% detector.simulate("printStatus", 1, "simulationData", "globalPFA");

%%
% Visualization
% detector.visualize("x_axis", "globalPFA", "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);
% detector.visualize("x_axis", "numberOfSensors", "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);
% detector.visualize("x_axis", "SNR", "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);
% detector.visualize("x_axis", "globalThreshold", "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);
detector.visualize("x_axis", "localPFA", "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);
detector.visualize("x_axis", "localThreshold", "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);
detector.visualize("x_axis", "localPD", "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);

% detector.visualize("x_axis", "globalPFA", "y_axis", "globalPFA", "dataType", ["analytical", "empirical"]);
% detector.visualize("x_axis", "numberOfSensors", "y_axis", "globalPFA", "dataType", ["analytical", "empirical"]);
% detector.visualize("x_axis", "SNR", "y_axis", "globalPFA", "dataType", ["analytical", "empirical"]);
% detector.visualize("x_axis", "localPFA", "y_axis", "globalPFA", "dataType", ["analytical", "empirical"]);


%%%% TO DO
% weighting
% different SNR

% 0.75/M
% SNR/M ile SNR vs Pd
% 1/M
% Binary Integration, 1, 0 1 bit, biz 32 bit iz ??