M_list = 2:2:20;
SNR_dB = 0:2:20;
PFA_local = 1e-2;
PFA_global = 1e-6;

[M_vals, PD_matrix, eta_vals] = sweep_PD_vs_M(M_list, SNR_dB, PFA_local, PFA_global);

% Plotting
figure;
hold on;
for i = 1:length(SNR_dB)
    plot(M_vals, PD_matrix(:, i), 'DisplayName', sprintf('SNR = %d dB', SNR_dB(i)));
end
grid on;
xlabel('Number of Receivers (M)');
ylabel('Probability of Detection');
title('Selective Noncoherent Detector: $P_D$ vs. $M$');
legend show;

function [M_vals, PD_matrix, eta_vals] = sweep_PD_vs_M(M_list, SNR_dB, PFA_local, PFA_global)
    % Inputs:
    %   M_list      - array of receiver counts
    %   SNR_dB      - array of SNRs in dB
    %   PFA_local   - local false alarm probability
    %   PFA_global  - desired global false alarm probability
    % Outputs:
    %   M_vals      - receiver counts (same as M_list)
    %   PD_matrix   - PD(i,j): i = M index, j = SNR index
    %   eta_vals    - global threshold for each M
    
    num_M = length(M_list);
    num_SNR = length(SNR_dB);
    PD_matrix = zeros(num_M, num_SNR);
    eta_vals = zeros(num_M, 1);
    
    for i = 1:num_M
        M = M_list(i);
        [eta, PD] = selective_PD_with_threshold(M, SNR_dB, PFA_local, PFA_global);
        eta_vals(i) = eta;
        PD_matrix(i, :) = PD;
    end
    M_vals = M_list;
end

function [eta_global, PD] = selective_PD_with_threshold(M, SNR_dB, PFA_local, PFA_global)
    % Inputs:
    %   M            - number of receivers
    %   SNR_dB       - array of SNR values in dB
    %   PFA_local    - local false alarm probability
    %   PFA_global   - target global false alarm probability
    % Outputs:
    %   eta_global   - global threshold to achieve PFA_global
    %   PD           - analytical probability of detection at each SNR

    % Step 1: Compute local threshold
    eta_m = -log(PFA_local);

    % Step 2: Compute global threshold eta such that P_FA = PFA_global
    func = @(eta) total_PFA(eta, M, PFA_local) - PFA_global;
    eta_global = fzero(func, [0.01, 100]);  % root finding

    % Step 3: Compute P_D for each SNR
    SNR = 10.^(SNR_dB / 10);
    PD = zeros(size(SNR));
    
    for idx = 1:length(SNR)
        snr = SNR(idx);
        theta = 1 / (1 + snr);  % scale
        pd_local = exp(-eta_m / (1 + snr));  % local detection prob
        
        pd_total = 0;
        for k = 1:M
            binom_weight = nchoosek(M, k) * pd_local^k * (1 - pd_local)^(M - k);
            gamma_cdf = gammainc(eta_global / theta, k, 'lower');
            tail_prob = 1 - gamma_cdf;
            pd_total = pd_total + binom_weight * tail_prob;
        end
        PD(idx) = pd_total;
    end
end

function PFA = total_PFA(eta, M, PFA_local)
    % Computes global PFA for given eta (under H0)
    PFA = 0;
    for k = 1:M
        binom_weight = nchoosek(M, k) * PFA_local^k * (1 - PFA_local)^(M - k);
        gamma_cdf = gammainc(eta, k, 'lower');  % under H0: scale = 1
        tail_prob = 1 - gamma_cdf;
        PFA = PFA + binom_weight * tail_prob;
    end
end