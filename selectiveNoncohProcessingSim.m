%%% Anıl ARSLAN 2303980
clc;

M = [1 2 3 4 5 7 9 11 14 18];
M = 1 : 55;
SNR_in_dB = 10*log10(10.^(.1*20)./M); % dB
PFA_local = 1;
PFA_global = 1e-6;
nmc = 1e4;

onlyPFAsim = false;

numSNR = length(SNR_in_dB);
numRX = length(M);

maxNrx = 55;
globalThresholds = zeros(maxNrx, 1);
globalThresholdClosedForm = zeros(maxNrx, 1);
globalThresholdFullSummation = zeros(maxNrx, 1);
for i = 1 : maxNrx
    % globalThresholds(i) = noncohSimPFAselective(i, PFA_local./i, PFA_global, 1e6);
    globalThresholdClosedForm(i) = noncohPFAselective(i, PFA_local/i, PFA_global);
    globalThresholdFullSummation(i) = gammaincinv(PFA_global, i, 'upper');
end

localThreshold = -log(PFA_local./M); % to fix data rate

globalThresholdsChoosen = globalThresholdClosedForm;

%% Analytical ROC curve

SNR = 0 : .25 : 5; % dB
globalPDclosedForm = zeros(maxNrx, length(SNR));
for k = 1 : length(SNR)
    for i = 1 : maxNrx
        globalPDclosedForm(i, k) = noncohPDselective(i, PFA_local/i, globalThresholdsChoosen(i), 10.^(1*SNR(k)));
    end
end

figure; plot(1 : maxNrx, globalPDclosedForm);
xlabel('# RX'); xlim([1, max(M)]);
ylabel('P_{D}'); ylim([0, 1]);
leg = legend(num2str(SNR.'), 'Location', 'best');
title(leg, 'SNR_{in}');
title('Global probability od detection vs #receivers');
grid on; grid minor;

%% Simulation
PDselectiveProcessing = zeros(numRX, numSNR);
PDfullProcessing = zeros(numRX, numSNR);
PFAselectiveProcessing = zeros(numRX, 1);
PFAfullProcessing = zeros(numRX, 1);

for i = 1 : numRX
    fprintf('%d/%d\n', i, numRX)

    %%% signal model
        n = (randn(M(i), 1, nmc) + 1j*randn(M(i), 1, nmc))/sqrt(2);
    if ~onlyPFAsim
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

    %%% full processing PD simulation
        LH1 = sum(abs(s).^2, 1);
        PDfullProcessing(i, :) = sum(LH1 > globalThresholdFullSummation(M(i)), 3)/nmc;
    end

    %%% selective processing PFA simulation
        ind = abs(n).^2 > localThreshold(i);
        m = sum(ind, 1);
        LH0 = zeros(1, 1, nmc);
        for mc = find(m ~= 0).'
            LH0(1, 1, mc) = sum(abs(n(ind(:, 1, mc), 1, mc)).^2, 1);
        end
        PFAselectiveProcessing(i) = sum(LH0 > globalThresholdsChoosen(M(i)), 3)/nmc;

    %%% full processing PFA simulation
        LH0 = sum(abs(n).^2, 1);
        PFAfullProcessing(i) = sum(LH0 > globalThresholdFullSummation(M(i)), 3)/nmc;
end

%%
figure; plot(M, diag(PDselectiveProcessing), '-b');
hold on; plot(M, diag(PDfullProcessing), '--b');
xlabel('# RX'); xlim([1, max(M)]);
ylabel('P_D'); ylim([0, 1]);
title('P_D vs #receivers');
legend('selective', 'full');
grid on; grid minor;

figure; semilogy(M, PFAselectiveProcessing);
hold on; semilogy(M, PFAfullProcessing);
xlabel('# RX'); xlim([1, max(M)]);
ylabel('P_{FA}'); ylim([0, 1]);
title('P_{FA} vs #receivers');
legend('selective', 'full');
grid on; grid minor;

figure; semilogy(M, globalThresholdsChoosen(M));
hold on; semilogy(M, globalThresholdFullSummation(M));
xlabel('# RX'); xlim([1, max(M)]);
ylabel('threshold');
title('Global threshold vs #receivers');
legend('selective', 'full');
grid on; grid minor;

%%%% TO DO
% weighting
% different SNR

% 0.75/M
% SNR/M ile SNR vs Pd
% 1/M
% Binary Integration, 1, 0 1 bit, biz 32 bit iz ??