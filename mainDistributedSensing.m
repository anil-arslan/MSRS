%%% Anıl ARSLAN
% 2303980

% widely seperated receivers | spatially i.i.d circularly symmetric complex gaussian signal model
    %%% LRT
        % SLC: Square Law Combining (LRT under equal SNR)                               Implemented
        % WSLC: Weighted Sqaure Law Combining (LRT under different SNR)                 Implemented
    
    %%% GLRT
        % LDC: Log-Divergence Combining (GLRT w/ joint MLE)                             Implemented
        % GLDC: Generalized Log-Divergence Combining (GLRT w/ independent MLE)          *Unsuccesful: no closed form expression
    
    %%% Binary
        % BC: Binary Combining (LRT under equal SNR)                                    Implemented
        % CVBC: Chair-Varshney Binary Combining (LRT under different SNR)               Implemented
    
    %%% Other
        % NSLC: Normalized Square Law Combining                                         Implemented
        % WSLCbiased: Weighted Sqaure Law Combining without unbias term                 Implemented
        % SC: Selective Combining (Partial Selective Combining K = 1)                   Implemented
        % PSC: Partial Selective Combining (Combine largest K samples)                  *Unsuccesful: successful for K = 1 (SC), 2 and M (SLED)

% phase asynchronous closely spaced receivers | spatially i.i.d uniform phase, fully correlated amplitude

% phase synchronous closely spaced receivers | spatially fully correlated circularly symmetric complex gaussian signal model
    %%% LRT
        % EGC: Equal Gain Combining


% clc;
% clear;
% close all;
addpath(genpath([pwd '/library']));

detector = detectorNoncoherent( ...
    "globalPFA", 1e-6, ...
    "numberOfSensors", 1 : 15, ...
    ..."numberOfSensors", 10 , ...
    "SNR_input_dB", 6, ...
    ..."SNR_input_dB", 10, ...
    ..."SNR_input_dB", linspace(0, 20, 20), ...
    ..."localPFA", [logspace(0, -8, 9), 0.5] ...
    ..."localPFA", logspace(-1, -8, 8) ...
    "localPFA", logspace(0, -8, 9) ...
    ..."localPFA", 1 ...
    ..."localPFA", .1 ...
    );
detector.setalgorithm( ...
    "numberOfBinaryDetectionsRule", "notSpecified", ...
    ..."numberOfBinaryDetectionsRule", "userDefined", ...
    ..."userDefinedNumberOfDetections", 2, ...
    "binaryDetectionConstraint", "fixedGlobalPFA", ...
    "signalAmplitudeModel", "correlatedExponential", ...
    ..."signalAmplitudeModel", "decorrelatedExponential", ...
    ..."signalPhaseModel", "correlatedUniform", ...
    "signalPhaseModel", "decorrelatedUniform", ...
    "globalFusionRule", "EGC" ...
    );
% detector.randomSNRwithFixedAverage("rangeSNR_dB", 30);
% detector.randomLocalPFAwithFixedAverage("averageLocalPFA", 0.01);
% detector.setbudgetconstraint("constraint", ["dataRate", "transmittedPower"]);
% detector.setbudgetconstraint("constraint", "dataRate");
% detector.visualize( ...
%     "x_axis", ["numberOfSensors", "SNR"], ...
%     "x_axis", "numberOfSensors", ...
%     "y_axis", "globalPD", "dataType", "analytical");
detector.visualize( ...
    "x_axis", ["numberOfSensors", "SNR"], ...
    "x_axis", "numberOfSensors", ...
    ..."x_axis", ["localPFA", "globalThreshold"], ...
    "y_axis", "globalPFA", "dataType", "analytical");



%% PD simulator
detector.setmontecarlo("numberOfTrials", 1e5);
detector.simulate("printStatus", 1, "simulationData", "globalPD", "statistics", 1);
detector.visualize( ...
    "x_axis", ["globalPFA", "numberOfSensors", "SNR", "globalThreshold"], ...
    "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);
detector.visualize( ...
    "x_axis", ["localPFA", "localThreshold", "localPD"], ...
    "y_axis", "globalPD", "dataType", ["analytical", "empirical"]);

% detector.visualize( ...
%     "x_axis", ["globalPFA", "numberOfSensors", "SNR", "globalThreshold", "localPFA", "localThreshold", "localPD"], ...
%     "y_axis", "numberOfActiveSensorsExpectationUnderSignal", "dataType", ["analytical", "empirical"]);
% detector.visualize( ...
%     "x_axis", ["globalPFA", "numberOfSensors", "SNR", "globalThreshold", "localPFA", "localThreshold", "localPD"], ...
%     "y_axis", "numberOfActiveSensorsStandardDeviationUnderSignal", "dataType", ["analytical", "empirical"]);


%% High Monte Carlo PFA simulator
clc;
clear;

detector = detectorNoncoherent( ...
    "globalPFA", 1e-5, ...
    "numberOfSensors", 1 : 10, ...
    "SNR_input_dB", 6, ...
    ..."SNR_input_dB", flipud(linspace(14, 0, 10).'), ... % meaningful for SNR dependent weighted fusion rules e.g. WSLC
    ..."localPFA", 1 ...
    "localPFA", logspace(0, -8, 5) ...
    ..."localPFA", logspace(0, -2, 9) ...
    );
% detector.setconstraint;
detector.setmontecarlo("numberOfTrials", 2e7);
detector.setalgorithm( ...
    "numberOfBinaryDetectionsRule", "notSpecified", ...
    ..."numberOfBinaryDetectionsRule", "userDefined", ...
    ..."userDefinedNumberOfDetections", 2, ...
    "binaryDetectionConstraint", "fixedLocalPFA", ...
    "globalFusionRule", "EGC" ...
    );
detector.simulate("printStatus", 1, "simulationData", "globalPFA", "statistics", 1);

% Visualization
% close all;
detector.visualize( ...
    "x_axis", "numberOfSensors", ...
    "y_axis", "globalPFA", "dataType", ["analytical", "empirical"]);
% detector.visualize( ...
%     "x_axis", ["globalPFA", "numberOfSensors", "SNR", "globalThreshold", "localPFA", "localThreshold", "localPD"], ...
%     "y_axis", "globalPFA", "dataType", ["analytical", "empirical"]);



% detector.visualize( ...
%     "x_axis", ["globalPFA", "numberOfSensors", "SNR", "globalThreshold", "localPFA", "localThreshold", "localPD"], ...
%     "y_axis", "numberOfActiveSensorsExpectationUnderNoise", "dataType", ["analytical", "empirical"]);
% detector.visualize( ...
%     "x_axis", ["globalPFA", "numberOfSensors", "SNR", "globalThreshold", "localPFA", "localThreshold", "localPD"], ...
%     "y_axis", "numberOfActiveSensorsStandardDeviationUnderNoise", "dataType", ["analytical", "empirical"]);




%%

%%%% TO DO

    % average output SNR, deflection coefficient? energy loss at fixed p_d p_fa?
    % Mean p_d, random variance? mixture?
    % i.i.d yapma (comparison with spatial models like swerling 0 & 1 to see the advantage of spatial diversity)
    % algorithm & signal model comparison script

%%%% Questions
    % binary integration: 1 bit, others: float; 32/16/8 bit (Power = nonnegative real number, so also unsigned)
    % SNR/sqrt(M) [since there will M^2 samples]
    % censoring
    % simulated RCS histogram is exponential, what about histogram of a true data?




%%% DONE
    % is SLC optimum?
    % optimum rule under local thresholding? GLRT LRT
    % binary integration exact global PFA using local PFA
    % binary integration exact local PFA using global PFA
    % mean and std of #active sensors
    % different local PFA case | censoring
    % random SNR, random local PFA



%%% Swerling

% Swerling 0, scalar
% Swerling 1, complex gaussian zero mean with some variance (random phase, exponential power)
% Swerling 2, same but independent from pulse to pulse

% Swerling 3, complex gaussian nonzero mean and variance
% Swerling 4, same but independent from pulse to pulse


%%
% Parameters
T = 4;                          % Truncation threshold
N = 1e5;                        % Number of samples
s_vals = linspace(0, 40, 500);  % S-axis for empirical CDF

% Generate truncated Rayleigh(1) samples using raylrnd
A1 = [];
A2 = [];

while numel(A1) < N
    samples = raylrnd(1/sqrt(2), 2*N, 1);     % Oversample for efficiency
    samples = samples(samples >= sqrt(T));    % Apply truncation
    A1 = [A1; samples(1 : min(N - numel(A1), numel(samples)))];
end

while numel(A2) < N
    samples = raylrnd(1/sqrt(2), 2*N, 1);
    samples = samples(samples >= sqrt(T));
    A2 = [A2; samples(1 : min(N - numel(A2), numel(samples)))];
end

% Compute S = |Z1 + Z2|^2
S = abs(A1.*exp(1j*2*pi*rand(N, 1)) + A2.*exp(1j*2*pi*rand(N, 1))).^2;

% Empirical CDF
F_empirical = arrayfun(@(s) mean(S <= s), s_vals);

% Plot
figure; hold on;
plot(s_vals, F_empirical, 'LineWidth', 2);
plot(s_vals, gammainc(s_vals/2, 1, 'lower'), 'LineWidth', 2);
xlabel('s = |Z_1 + Z_2|^2');
ylabel('Empirical CDF');
title(sprintf('Monte Carlo CDF using raylrnd (T = %.2f)', T));
grid on;
