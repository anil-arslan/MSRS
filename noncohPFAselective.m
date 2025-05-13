function globalThreshold = noncohPFAselective(M, PFA_local, PFA_global)
    % Closed-form threshold estimation for selective noncoherent detection
    % Inputs:
    %   M           - number of receivers
    %   PFA_local   - local false alarm probability (per receiver)
    %   PFA_global  - desired global false alarm probability
    % Output:
    %   eta         - global threshold

    localThreshold = -log(PFA_local);
    p = PFA_local;

    % Helper: Q_k(x) = upper tail of Erlang(k,1), shifted by k*eta_m
    function q = Qk(k, x)
        if x < 0
            q = 1;
        else
            q = sum((x .^ (0:k-1)) ./ factorial(0:k-1)) * exp(-x);
        end
    end

    % Global PFA expression
    function pfa = globalPFA(eta)
        pfa = 0;
        for k = 1 : M
            pk = nchoosek(M, k) * p^k * (1 - p)^(M - k);
            x = eta - k * localThreshold;
            pfa = pfa + pk * Qk(k, x);
        end
    end

    % Solve for eta such that globalPFA(eta) == PFA_global
    eta_lb = M * localThreshold;
    eta_ub = eta_lb + 20;
    globalThreshold = fzero(@(x) globalPFA(x) - PFA_global, [eta_lb, eta_ub]);
end
