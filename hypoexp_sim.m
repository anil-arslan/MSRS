%%% Anıl ARSLAN
% 2303980
clc; clear; close all;

%%% parameters
weights = [0.4, 0.4, 0.4, 0.5 0.5, 0.45];
weights = [0.5, 0.5, 0.5];
nSamples = 100000;
tValues = linspace(0, 15, 500);

lambdas = 1./weights;
samples = zeros(nSamples, 1);
for i = 1 : length(lambdas)
    samples = samples + exprnd(1 / lambdas(i), nSamples, 1);
end

empiricalCCDF = arrayfun(@(t) mean(samples >= t), tValues);
analyticCCDF = arrayfun(@(t) hypoexponentialCCDF(t, lambdas), tValues);

%%% visualization
figure;
plot(tValues, empiricalCCDF, 'b', 'LineWidth', 2); hold on;
plot(tValues, analyticCCDF, 'r--', 'LineWidth', 2);
plot(tValues, gammainc(tValues.*unique(lambdas), 3, 'upper'), 'g-', 'LineWidth', 2);
legend('Empirical CDF (Simulation)', 'Analytical CDF (Matrix Exponential)');
xlabel('t'); ylabel('CDF');
title('Validation of Hypoexponential CDF (with repeated weights)');
grid on;

function probability = hypoexponentialCCDF(t, lambdas)
    % hypoexponentialCCDF - Computes the complementary CDF of a hypoexponential distribution
    % with potentially repeated weights using the matrix exponential method.
    mustBeNonnegative(lambdas);
    N = length(lambdas);
    
    Q = diag(-lambdas); % Construct generator matrix Q
    Q(N + 1 : N + 1 : end) = lambdas(1 : end - 1);

    alpha = [1 zeros(1, N - 1)]; % Initial distribution (starting in state 1)
    e = ones(N, 1); % Exit vector (absorbing state after last)
    
    probability = alpha*expm(Q*t)*e;
end