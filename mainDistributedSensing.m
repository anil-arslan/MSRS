%%% Anıl ARSLAN 2303980
clc; % clear;
addpath(genpath([pwd '/library']));

detector = detectorNoncoherent( ...
    "globalPFA", 1e-6, ...
    "numberOfSensors", 1 : 20, ...
    "SNR_input_dB", 14, ...
    "localPFA", logspace(0, -8, 9) ...
    );
detector.setconstraint;
detector.setmontecarlo("numberOfTrials", 1e4);
detector.setalgorithm("globalFusionRule", "EGN");
detector.simulate("printStatus", 1, "simulationData", "globalPD");
%%% High Monte Carlo PFA simulator
% detector.setmontecarlo("numberOfTrials", 1e5);
% detector.simulate("printStatus", 1, "simulationData", "globalPFA");


%%% Visualization
% close all;
detector.visualize( ...
    "x_axis", ["globalPFA", "numberOfSensors", "SNR", "globalThreshold"], ...
    "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);
% detector.visualize( ...
%     "x_axis", ["localPFA", "localThreshold", "localPD"], ...
%     "y_axis", "globalPD", "dataType", "analytical");



%%%% TO DO
% EGN vs Normalized EGN, closed form
% mean/variance formullerini objeye koyalim

% weighting
% different SNR

% Binary Integration, 1, 0 1 bit, biz 32/16/8 bit iz ?? (Power = nonnegative real number)

% SNR/sqrt(M) [cunku M^2 kadar samples olacak]