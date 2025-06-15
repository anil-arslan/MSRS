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

detector = detectorNoncoherent( ...
    "globalPFA", 1e-6, ...
    "numberOfSensors", 1 : 15, ...
    ..."numberOfSensors", 10 , ...
    "SNR_input_dB", 6, ...
    ..."SNR_input_dB", 10, ...
    ..."SNR_input_dB", linspace(0, 20, 20), ...
    ..."localPFA", [logspace(0, -8, 9), 0.5] ...
    ..."localPFA", logspace(-1, -8, 8) ...
    ..."localPFA", logspace(0, -8, 9) ...
    "localPFA", 1 ...
    ..."localPFA", .1 ...
    );
detector.setalgorithm( ...
    "numberOfBinaryDetectionsRule", "notSpecified", ...
    ..."numberOfBinaryDetectionsRule", "userDefined", ...
    ..."userDefinedNumberOfDetections", 2, ...
    "binaryDetectionConstraint", "fixedGlobalPFA", ...
    "signalAmplitudeModel", "exponential", ...
    ..."signalAmplitudeModel", "independentExponential", ...
    "signalPhaseModel", "independentUniform", ...
    ..."signalPhaseModel", "uniform", ...
    "globalFusionRule", "EGC" ...
    );
% detector.setuserdefinednumberofdetections();
% detector.randomSNRwithFixedAverage("rangeSNR_dB", 30);
% detector.setbudgetconstraint("constraint", ["dataRate", "transmittedPower"]);
% detector.setbudgetconstraint("constraint", "dataRate");
detector.visualize( ...
    "x_axis", ["numberOfSensors", "SNR"], ...
    "x_axis", "numberOfSensors", ...
    "y_axis", "globalPD", "dataType", "analytical");
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
    "SNR_input_dB", 14, ...
    "SNR_input_dB", flipud(linspace(14, 0, 10).'), ... % meaningful for SNR dependent weighted fusion rules e.g. MRC
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
    "globalFusionRule", "BC" ...
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

    % average output SNR, deflection coefficient? energy loss at fixed pd pfa?
    % Mean p_d, random variance? mixture?
    % independent yapma (comparison with swerling 0 ve swerling 1 to see the advantage of spatial diversity)
    % different local PFA case
    % algorithm comparison dimension

%%%% Questions
    % binary integration: 1 bit, others: float; 32/16/8 bit (Power = nonnegative real number, so also unsigned)
    % SNR/sqrt(M) [since there will M^2 samples]
    % censoring




%%% DONE
    % is EGC optimum?
    % optimum rule under local thresholding? GLRT LRT
    % binary integration exact global PFA using local PFA



%%% Swerling

% Swerling 0, scalar
% Swerling 1, complex gaussian zero mean with some variance (random phase, exponential power)
% Swerling 2, same but independent from pulse to pulse
% Swerling 3, complex gaussian nonzero mean and variance
% Swerling 4, same but independent from pulse to pulse


%%
clc; clear;

% Parameters
M = 1;                  % Number of samples
sigma = 1;            % Rayleigh scale parameter
Nmc = 5e5;              % Monte Carlo simulation size
t_vals = linspace(0, 20, 300);  % Evaluation grid

% Simulation
A = raylrnd(sigma/sqrt(2), [1, Nmc]);
phi = 2*pi*rand(M, Nmc);
signal = A.*exp(1j*phi);
noise = (randn(M, Nmc) + 1j * randn(M, Nmc))/sqrt(2);
x = signal + noise;
T_sim = sum(abs(x).^2, 1);
cdf_sim = arrayfun(@(t) mean(T_sim <= t), t_vals);

% Theoretical CDF via integral
cdf_th = zeros(size(t_vals));

for k = 1 : length(t_vals)
    integrand = @(a) ncx2cdf(t_vals(k), 2*M, a).*exppdf(a, sigma);
    integrand = @(a) ncx2cdf(t_vals(k), 2*M, a.^2).*raylpdf(a, sigma/sqrt(2));
    cdf_th(k) = integral(integrand, 0, inf);
end


cdf_exp = 1 - exp(-t_vals / (1 + sigma^2));  % Analytical exponential CDF

% Plot
figure; hold on; grid on;
plot(t_vals, cdf_sim, 'b', 'LineWidth', 2);
plot(t_vals, cdf_exp, 'g--', 'LineWidth', 2);
plot(t_vals, cdf_th, 'r--', 'LineWidth', 2);
xlabel('t'); ylabel('CDF');
legend('Empirical CDF', 'Exponential CDF', 'Theoretical CDF');
title(sprintf('Validation of Energy Detector CDF, M = %d', M));

%%
clc; clear;
compare_pdf_with_montecarlo_final_final;

function compare_pdf_with_montecarlo_final_final()
    lambda = 3;
    t_vals = linspace(0.01, 15, 300);
    dt = t_vals(2) - t_vals(1);

    % Numerical integral (robust)
    pdf_vals = magnitude_square_pdf_safe(t_vals, lambda);

    % Monte Carlo Simulation
    N = 1e6;
    A = exprnd(1/lambda, N, 1);
    Z = (randn(N,1) + 1i*randn(N,1)) / sqrt(2);
    Y = sqrt(A) + Z;
    T_samples = abs(Y).^2;

    % Histogram
    [counts, ~] = histcounts(T_samples, [t_vals, t_vals(end)+dt]);
    mc_pdf = counts / (N * dt);

    % Plot
    figure;
    plot(t_vals, pdf_vals, 'b-', 'LineWidth', 2); hold on;
    stairs(t_vals, mc_pdf, 'r--', 'LineWidth', 1.5);
    legend('Numerical Integration (safe)', 'Monte Carlo');
    xlabel('t = |Y|^2'); ylabel('PDF');
    title('PDF of |sqrt{A} + Z|^2 (Fully Stable)');
    grid on;
end

function ft = magnitude_square_pdf_safe(t_vec, lambda)
    ft = zeros(size(t_vec));
    for idx = 1:length(t_vec)
        t = t_vec(idx);
        sqrt_t = sqrt(t);
        integrand = @(u) stable_integrand(u, sqrt_t, lambda);
        val = integral(integrand, 0, Inf, 'RelTol',1e-8,'AbsTol',1e-12);
        ft(idx) = lambda * exp(-t) * val*2;
    end
end

function val = stable_integrand(u, sqrt_t, lambda)
    x = 2 * u * sqrt_t;

    % Use exact Bessel where safe, asymptotic elsewhere
    val = zeros(size(x));
    safe = x < 700;
    val(safe) = besseli(0, x(safe));

    % Asymptotic for large x: I_0(x) ≈ exp(x) / sqrt(2πx)
    large = ~safe;
    val(large) = exp(x(large)) ./ sqrt(2 * pi * x(large));

    % Final expression
    val = u .* exp(-(1 + lambda) * u.^2) .* val;
    
    % Ensure no NaNs or Infs
    val(~isfinite(val)) = 0;
end



%%
% Parameters
lambda = 2;
t_vals = linspace(0, 10, 200);
F_analytic = magnitude_square_cdf_ncx2_fixed(t_vals, lambda);

% Monte Carlo
N = 1e6;
A = exprnd(1/lambda, N, 1);
Z = (randn(N,1) + 1i*randn(N,1))/sqrt(2); % CN(0,1)
T = abs(sqrt(A) + Z).^2;
F_mc = arrayfun(@(t) mean(T <= t), t_vals);

% Plot
plot(t_vals, F_analytic, 'b-', 'LineWidth', 2); hold on;
plot(t_vals, F_mc, 'r--', 'LineWidth', 1.5);
legend('Corrected Analytical (ncx2cdf)', 'Monte Carlo');
xlabel('t'); ylabel('CDF');
title('CDF of |sqrt{A} + Z|^2');
grid on;


function F = magnitude_square_cdf_ncx2_fixed(t_vals, lambda)
    % Corrected: Rescale argument and noncentrality for MATLAB's ncx2cdf
    F = zeros(size(t_vals));
    for idx = 1:length(t_vals)
        t = t_vals(idx);
        integrand = @(a) ncx2cdf(2*t, 2, 2*a) .* lambda .* exp(-lambda * a);
        F(idx) = integral(integrand, 0, Inf, 'RelTol',1e-8, 'AbsTol',1e-12);
    end
end


%%

lambda = .11;
M = 5;
t_vals = linspace(0, 100, 300);

F_series = magnitude_square_sum_cdf_series_factored(t_vals, lambda, M, 100);
F_integral = magnitude_square_sum_cdf(t_vals, lambda, M);  % from previous function

% Monte Carlo
N = 1e6;
A = exprnd(1/lambda, N, M);
Z = (randn(N, M) + 1i*randn(N, M)) / sqrt(2);
T = sum(abs(sqrt(A) + Z).^2, 2);
F_mc = arrayfun(@(t) mean(T <= t), t_vals);

% Plot
plot(t_vals, F_series, 'b-', 'LineWidth', 2); hold on;
plot(t_vals, F_integral, 'k-.', 'LineWidth', 2);
plot(t_vals, F_mc, 'r--', 'LineWidth', 1.2);
legend('Series Expansion', 'Numerical Integral', 'Monte Carlo');
xlabel('t'); ylabel('CDF');
title(['CDF of sum_{i=1}^{', num2str(M), '} |sqrt{A_i} + Z_i|^2']);
grid on;



function F = magnitude_square_sum_cdf(t_vals, lambda, M)
    % Computes CDF of T = sum_{i=1}^M |sqrt{A_i} + Z_i|^2
    % where A_i ~ Exp(lambda), Z_i ~ CN(0,1)

    F = zeros(size(t_vals));

    % Define the gamma density for sum S = sum a_i ~ Gamma(M, lambda)

    for idx = 1:length(t_vals)
        t = t_vals(idx);

        % Integrand: ncx2cdf(2t, 2M, 2s) * gamma_pdf(s)
        integrand = @(s) ncx2cdf(2*t, 2*M, 2*s) .* gampdf(s, M, 1/lambda);

        % Numerical integration over s
        F(idx) = integral(integrand, 0, Inf, 'RelTol',1e-8, 'AbsTol',1e-12);
    end
end

function F = magnitude_square_sum_cdf_series_factored(t_vals, lambda, M, K_max)
    % Computes the CDF of T = sum_{i=1}^M |sqrt{A_i} + Z_i|^2
    % using the cleaned-up series expansion
    % Inputs:
    %   t_vals  - vector of t values
    %   lambda  - exponential rate parameter
    %   M       - number of independent terms
    %   K_max   - number of terms in the summation

    if nargin < 4
        K_max = 100;  % Default number of terms if not specified
    end

    F = zeros(size(t_vals));  % Initialize output
    pre_factor = (lambda^M) / ((lambda + 1)^M * factorial(M - 1));  % Outside the sum

    for idx = 1:length(t_vals)
        t = t_vals(idx);
        sum_k = 0;

        for k = 0:K_max
            % Regularized incomplete gamma multiplied by Gamma gives gamma()
            gam = gammainc(t, M + k, 'lower') * gamma(M + k);
            term = gam / (factorial(k) * (lambda + 1)^k);
            sum_k = sum_k + term;
        end

        F(idx) = pre_factor * sum_k;
    end
end


