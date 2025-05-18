%%% Anıl ARSLAN 2303980
clc;

M = 1 : 1 : 10;
SNR_in_dB = 0 : 10; % dB
PFA_local = 1;
PFA_global = 1e-6;
nmc = 1e4;

numSNR = length(SNR_in_dB);
numRX = length(M);

maxNrx = max(M);
globalThresholdClosedForm = zeros(maxNrx, 1);
globalPDanalytical= zeros(maxNrx, numSNR);
for i = 1 : maxNrx
    globalThresholdClosedForm(i) = noncohPFAselective(i, PFA_local/i, PFA_global);
    for k = 1 : numSNR
        globalPDanalytical(i, k) = noncohPDselective(i, PFA_local/i, globalThresholdClosedForm(i), 10.^(.1*SNR_in_dB(k)));
    end
end

localThreshold = -log(PFA_local./M); % to fix data rate
globalThresholdsChoosen = globalThresholdClosedForm;

%% Simulation
PDselectiveProcessing = zeros(numRX, numSNR);
PFAselectiveProcessing = zeros(numRX, 1);

for i = 1 : numRX
    fprintf('%d/%d\n', i, numRX)

    %%% signal model
        n = (randn(M(i), 1, nmc) + 1j*randn(M(i), 1, nmc))/sqrt(2);
        x = (randn(M(i), 1, nmc) + 1j*randn(M(i), 1, nmc)).*10.^(.05*SNR_in_dB)/sqrt(2);
        s = x + n; % [M by 1]

    %%% selective processing PD simulation
        ind = abs(s).^2 > localThreshold(i);
        m = sum(ind, 1);
        LH1 = zeros(1, numSNR, nmc);
        for q = 1 : numSNR
            for mc = find(m(:, q, :) ~= 0).'
                LH1(1, q, mc) = sum(abs(s(ind(:, q, mc), q, mc)).^2);
            end
        end
        PDselectiveProcessing(i, :) = sum(LH1 > globalThresholdsChoosen(M(i)), 3)/nmc;

    %%% selective processing PFA simulation
        ind = abs(n).^2 > localThreshold(i);
        m = sum(ind, 1);
        LH0 = zeros(1, 1, nmc);
        for mc = find(m ~= 0).'
            LH0(1, 1, mc) = sum(abs(n(ind(:, 1, mc), 1, mc)).^2, 1);
        end
        PFAselectiveProcessing(i) = sum(LH0 > globalThresholdsChoosen(M(i)), 3)/nmc;
end

%%
figure; plot(M, PDselectiveProcessing, '-b');
hold on; plot(M, globalPDanalytical, '--r');
xlabel('# RX'); xlim([1, max(M)]);
ylabel('P_D'); ylim([0, 1]);
title('P_D vs #receivers');
grid on; grid minor;