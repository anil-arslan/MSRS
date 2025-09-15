function Pd_global = binaryCCDF(M, K, Pfa_local, rho)
    % Pd_binary returns the global Pd for a k-out-of-M binary combining rule.
    % M           : number of sensors
    % k           : rule parameter (k-out-of-M)
    % Pfa_local   : local false alarm probability
    % rho         : SNR (linear)
    
    % Local detection probability as a function of SNR and Pfa_local
    Pd_local = (Pfa_local).^(1./(1 + rho));
    
    % Global Pd using binomial summation
    Pd_global = zeros(length(K), length(Pfa_local));
    for j = 1 : length(K)
        k = K(j);
        for i = k : M
            Pd_global(j, :) = Pd_global(j, :) + nchoosek(M, i).*Pd_local.^i.*(1 - Pd_local).^(M - i);
        end
    end
end