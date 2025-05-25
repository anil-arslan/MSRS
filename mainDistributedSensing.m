%%% Anıl ARSLAN 2303980
clc;
% clear;
% close all;
addpath(genpath([pwd '/library']));

detector = detectorNoncoherent( ...
    "globalPFA", 1e-6, ...
    "numberOfSensors", 1 : 30, ...
    "SNR_input_dB", flipud(linspace(0, 10, 30).'), ...
    "localPFA", logspace(0, -8, 9) ...
    );
% detector.setbudgetconstraint("constraint", ["dataRate", "transmittedPower"]);
detector.setalgorithm("globalFusionRule", "EGC", "numberOfSensorsPSC", 5);
detector.visualize( ...
    "x_axis", ["numberOfSensors", "SNR"], ...
    "x_axis", "numberOfSensors", ...
    "y_axis", "globalPD", "dataType", "analytical");

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
    "numberOfSensors", 1 : 15, ...
    "SNR_input_dB", 14, ...
    "localPFA", logspace(0, -8, 5), ...
    "localPFA", logspace(0, -2, 9) ...
    );
% detector.setconstraint;
detector.setmontecarlo("numberOfTrials", 1e7);
detector.setalgorithm("globalFusionRule", "PSC");
detector.simulate("printStatus", 1, "simulationData", "globalPFA");

% Visualization
% close all;
detector.visualize( ...
    "x_axis", ["globalPFA", "numberOfSensors", "SNR", "globalThreshold", "localPFA", "localThreshold", "localPD"], ...
    "y_axis", "globalPFA", "dataType", ["analytical", "empirical"]);

%%

%%%% TO DO
% mean/variance formullerini objeye koyalim
% average output SNR

% what is mean p_d, conditioned on given SNR?

% weighting, MRC
% different SNR

% Binary Integration, 1, 0 1 bit, biz 32/16/8 bit iz ?? (Power = nonnegative real number)

% SNR/sqrt(M) [cunku M^2 kadar samples olacak]

% Local thresholding de optimum rule nedir?
% EGC optimum mudur?







%% Calisma
M = 11;
K = 4;

F = @(t) 1 - sum(arrayfun(@(j) (-1).^j*nchoosek(M - K, j)*gammainc(t, K + j, 'lower'), 0 : (M - K)));

x = linspace(0, 30, 1000);

for i = 1 : length(x)
    % c(i) = F(x(i))*K*nchoosek(M, K);
    c(i) = F(x(i));
end

figure; plot(x, c);
