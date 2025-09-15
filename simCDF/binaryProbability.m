function [prob, q_local] = binaryProbability(M, pfaLocal, pfaGlobal, rho)
    K = 1 : M;

    rulesPfaLocal = zeros(1, M);
    for k = 1 : M
        rulesPfaLocal(k) = fzero(@(pfaLoc) binaryCCDF(M, k, pfaLoc, 0) - pfaGlobal, [0 1]);
    end
    pRules = diag(binaryCCDF(M, K, rulesPfaLocal, 10^(rho/10)));
    
    q_local = zeros(M - 1, length(pfaLocal));
    for i = 2 : M
        q_local(i - 1, :) = (pfaLocal - rulesPfaLocal(i))./(rulesPfaLocal(i - 1) - rulesPfaLocal(i));
    end
    q_local(q_local > 1 | q_local < 0) = nan;

    lowLocalPFAregion = find(pfaLocal <= rulesPfaLocal(1), 1, 'first') : length(pfaLocal);
    highLocalPFAregion = 1 : find(pfaLocal > rulesPfaLocal(end), 1, 'last');
    
    % pfaLocall = q_local.*rulesPfaLocal(1 : (M - 1)).' + (1 - q_local).*rulesPfaLocal(2 : M).';
    % figure; loglog(pfaLocal, pfaLocall.');
    
    prob = q_local.*pRules(1 : (M - 1)) +  (1 - q_local).*pRules(2 : M);

    ORrule = binaryCCDF(M, 1, pfaLocal(lowLocalPFAregion), 0);
    ANDrule = binaryCCDF(M, M, pfaLocal(highLocalPFAregion), 0);
    q_lowLocalPFA = (pfaGlobal - ORrule)./(1 - ORrule);
    q_highLocalPFA = pfaGlobal./ANDrule;

    ORrule = binaryCCDF(M, 1, pfaLocal(lowLocalPFAregion), 10^(rho/10));
    ANDrule = binaryCCDF(M, M, pfaLocal(highLocalPFAregion), 10^(rho/10));
    probLowLocalPFAregion = nan(1, length(pfaLocal));
    probHighLocalPFAregion = nan(1, length(pfaLocal));
    
    probLowLocalPFAregion(lowLocalPFAregion) = q_lowLocalPFA + (1 - q_lowLocalPFA).*ORrule;
    probHighLocalPFAregion(highLocalPFAregion) = q_highLocalPFA.*ANDrule;



    prob = [probLowLocalPFAregion; prob; probHighLocalPFAregion];

    q_lowLocalPFAvec = nan(1, length(pfaLocal));
    q_highLocalPFAvec = nan(1, length(pfaLocal));

    q_lowLocalPFAvec(lowLocalPFAregion) = q_lowLocalPFA;
    q_highLocalPFAvec(highLocalPFAregion) = q_highLocalPFA;
    q_local = [q_lowLocalPFAvec; q_local; q_highLocalPFAvec];
end