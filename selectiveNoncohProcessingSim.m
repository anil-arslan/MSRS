%%% Anıl ARSLAN 2303980
clc;

M = [1 2 3 4 5 7 9 11 14 18];
SNR_in_dB = 0 : .25 : 5; % dB
PFA_local = 0.4;
PFA_global = 1e-5;
nmc = 1e6;

numSNR = length(SNR_in_dB);
numRX = length(M);

maxNrx = 20;
globalThresholds = zeros(maxNrx, 1);
globalThresholdClosedForm = zeros(maxNrx, 1);
globalThresholdFullSummation = zeros(maxNrx, 1);
for i = 1 : size(globalThresholds, 1)
    globalThresholds(i) = noncohSimPFAselective(i, PFA_local./i, PFA_global, 1e6);
    % globalThresholdClosedForm(i) = noncohPFAselective(i, PFA_local, PFA_global);
    globalThresholdFullSummation(i) = gammaincinv(PFA_global, i, 'upper');
end

localThreshold = -log(PFA_local./M); % to fix data rate

%%
PDselectiveProcessing = zeros(numRX, numSNR);
PDfullProcessing = zeros(numRX, numSNR);
PFAselectiveProcessing = zeros(numRX, 1);
PFAfullProcessing = zeros(numRX, 1);

for i = 1 : numRX
    fprintf('%d/%d\n', i, numRX)

    %%% signal model
        n = (randn(M(i), 1, nmc) + 1j*randn(M(i), 1, nmc))/sqrt(2);
        x = (randn(M(i), 1, nmc) + 1j*randn(M(i), 1, nmc)).*10.^(.5*SNR_in_dB)/sqrt(2);
        s = x + n; % [M by 1]

    %%% selective processing PD simulation
        ind = abs(s).^2 > localThreshold(i);
        m = sum(ind, 1);
        LH1 = zeros(1, numSNR, nmc);
        globalThreshold = zeros(1, numSNR, nmc);
        for q = 1 : numSNR
            for mc = find(m(:, q, :) ~= 0).'
                globalThreshold(1, q, mc) = globalThresholds(m(:, q, mc)); % correctly works
                LH1(1, q, mc) = sum(abs(s(ind(:, q, mc), q, mc)).^2);
            end
        end
        PDselectiveProcessing(i, :) = sum(LH1 > globalThreshold, 3)/nmc;

    %%% full processing PD simulation
        globalThreshold = globalThresholdFullSummation(M(i));
        LH1 = sum(abs(s).^2, 1);
        PDfullProcessing(i, :) = sum(LH1 > globalThreshold, 3)/nmc;

    %%% selective processing PFA simulation
        ind = abs(n).^2 > localThreshold(i);
        m = sum(ind, 1);
        LH0 = zeros(1, 1, nmc);
        globalThreshold = zeros(1, 1, nmc);
        for mc = find(m ~= 0).'
            % globalThreshold(1, 1, mc) = globalThresholdClosedForm(m(mc)); % does not fix PFA
            % globalThreshold(1, 1, mc) = globalThresholdFullSummation(m(mc)); % does not fix PFA
            globalThreshold(1, 1, mc) = globalThresholds(m(mc)); % correctly works
            LH0(1, 1, mc) = sum(abs(n(ind(:, 1, mc), 1, mc)).^2, 1);
        end
        PFAselectiveProcessing(i) = sum(LH0 > globalThreshold, 3)/nmc;

    %%% full processing PFA simulation
        globalThreshold = globalThresholdFullSummation(M(i));
        LH0 = sum(abs(n).^2, 1);
        PFAfullProcessing(i) = sum(LH0 > globalThreshold, 3)/nmc;
end

%%
figure; plot(M, PDselectiveProcessing);
hold on; plot(M, PDfullProcessing, '--');
xlabel('# RX'); xlim([1, max(M)]);
ylabel('P_D'); ylim([0, 1]);
title('P_D vs #receivers');
legend(num2str(SNR_in_dB.'));
grid on; grid minor;

figure; semilogy(M, PFAselectiveProcessing);
hold on; semilogy(M, PFAfullProcessing);
xlabel('# RX'); xlim([1, max(M)]);
ylabel('P_{FA}'); ylim([0, 1]);
title('P_{FA} vs #receivers');
legend('selective', 'full');
grid on; grid minor;

figure; semilogy(M, globalThresholds(M));
hold on; semilogy(M, globalThresholdFullSummation(M));
xlabel('# RX'); xlim([1, max(M)]);
ylabel('threshold');
title('Global threshold vs #receivers');
legend('selective', 'full');
grid on; grid minor;