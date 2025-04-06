clc;
N = 1 : 10;
PFA = logspace(-8, 0, 81);
T = zeros(length(PFA), length(N));
for k = 1 : length(N)
    T(:, k) = 10*log10(gammaincinv(PFA, N(k), 'upper'));
end

figure; semilogx(PFA, T);
grid on; grid minor;
ylabel('threshold (dB)'); xlabel('PFA');