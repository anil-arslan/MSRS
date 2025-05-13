function globalThreshold = noncohSimPFAselective(M, PFA_local, PFA_global, num_trials)
    % Empirical threshold computation for selective noncoherent detection
    % Inputs:
    %   M            - number of receivers
    %   PFA_local    - local false alarm rate (used to threshold T_m)
    %   PFA_global   - desired global false alarm rate
    %   num_trials   - number of Monte Carlo trials
    % Output:
    %   eta_empirical - empirically derived threshold achieving PFA_global

    % Local threshold
    localThreshold = -log(PFA_local);

    % Preallocate test statistic under H0
    T_sel_H0 = zeros(num_trials, 1);

    for i = 1:num_trials
        Tm = exprnd(1, M, 1);              % T_m ~ Exp(1) under H0
        % delta = Tm > localThreshold;                % apply local threshold
        T_sel_H0(i) = sum(Tm + localThreshold);      % sum only passed values
    end

    % Empirical threshold at (1 - P_FA) percentile
    sorted_T = sort(T_sel_H0);
    idx = max(1, round((1 - PFA_global) * num_trials));
    globalThreshold = sorted_T(idx);
end
