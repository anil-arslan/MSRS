classdef detector < handle
    %detector Summary of this class goes here
    %   Detailed explanation goes here

    %%% Model parameters

    properties (SetAccess = private, GetAccess = public)
        seed (1, 1) double {mustBeNonnegative, mustBeInteger, mustBeInRange(seed, 0, 4294967295)} = 0
        numberOfTrials (1, 1) double {mustBeNonnegative, mustBeInteger} = 1e3
        numberOfLoops (1, 1) double {mustBeNonnegative, mustBeInteger} = 1
        SNR_input_dB (1, :) cell = {4} % [1 x M] cell of [M x Nsnr] matrices
        numberOfSensors (1, :) double {mustBeNonnegative, mustBeInteger} = 10 % [1 x M] vector
    end

    %%% Fusion parameters

    properties (SetAccess = private, GetAccess = public)
        numberOfRequiredDetections double {mustBeNonnegative, mustBeInteger} = 1 % [1 x M] vector
        globalFusionRule (1, 1) string {mustBeMember(globalFusionRule, ["SLC", "WSLC", "BC", "CVBC", "LLC", "WLLC", "NSLC", "SC", "PSC", "EGC", "MRC"])} = "SLC"
        %%% LRT: widely seperated receivers | spatially i.i.d circularly symmetric complex gaussian signal model
            % SLC: Square Law Combining (LRT under equal SNR)
            % WSLC: Weighted Sqaure Law Combining (LRT under different SNR)

        %%% Binary
            % BC: Binary Combining (LRT under equal SNR)
            % CVBC: Chair-Varshney Binary Combining (LRT under different SNR)

        %%% Other
            % LLC: Linear Law Combining
            % WLLC: Weighted Linear Law Combining
            % NSLC: Normalized Square Law Combining
            % SC: Selective Combining (Partial Selective Combining K = 1)
            % PSC: Partial Selective Combining (Combine largest K samples)

        %%% LRT: phase asynchronous closely spaced receivers | spatially i.i.d uniform phase, fully correlated amplitude

        %%% LRT: phase synchronous closely spaced receivers | spatially fully correlated circularly symmetric complex gaussian signal model
            % EGC: Equal Gain Combining
            % MRC: Matched Filter

        signalAmplitudeModel (1, 1) string {mustBeMember(signalAmplitudeModel, ["decorrelatedExponential", "correlatedExponential", "deterministic"])} = "decorrelatedExponential"
        signalPhaseModel (1, 1) string {mustBeMember(signalPhaseModel, ["decorrelatedUniform", "correlatedUniform"])} = "decorrelatedUniform"
        binaryDetectionRule (1, 1) string {mustBeMember(binaryDetectionRule, ["and", "or", "majority", "notSpecified", "userDefined"])} = "notSpecified"
        binaryDetectionPFAtype (1, 1) string {mustBeMember(binaryDetectionPFAtype, ["fixedGlobalPFA", "fixedLocalPFA", "fixedGlobal|LocalPFA", "timeSharing"])} = "fixedLocalPFA"

        fusionWeights (1, :) cell % [1 x M] cell of [M x Nsnr x 1 x Nlocal] matrices
    end

    properties (SetAccess = private, GetAccess = public)
        globalPFA (:, 1) double {mustBeNonnegative, mustBeInRange(globalPFA, 0, 1)} = 1e-6 % [NpfaGlobal x 1] vector
        globalPFAsimulation double {mustBeNonnegative, mustBeInRange(globalPFAsimulation, 0, 1)} = [] % [NpfaGlobal x Nsnr x M x Nlocal] matrix
    end

    properties (Dependent)
        globalThreshold double {mustBeNonnegative} % [NpfaGlobal x Nsnr x M x Nlocal] matrix
    end

    properties (SetAccess = private, GetAccess = public)
        globalNeighbourThresholdCVBC double {mustBeNonnegative} = [] % [NpfaGlobal x Nsnr x M x Nlocal] matrix
        globalThresholdSimulation double {mustBeNonnegative} = [] % [NpfaGlobal x Nsnr x M x Nlocal] matrix
    end

    properties (Dependent)
        globalPFAanalytical double {mustBeNonnegative} % [NpfaGlobal x Nsnr x M x Nlocal] matrix
        globalPD double {mustBeNonnegative, mustBeInRange(globalPD, 0, 1)} % [NpfaGlobal x Nsnr x M x Nlocal] matrix
    end

    properties (SetAccess = private, GetAccess = public)
        % global probability of random threshold
        globalRandomizationProbability double {mustBeNonnegative, mustBeInRange(globalRandomizationProbability, 0, 1)} = [] % [NpfaGlobal x Nsnr x M x Nlocal] matrix
        globalPDsimulation double {mustBeNonnegative, mustBeInRange(globalPDsimulation, 0, 1)} = [] % [NpfaGlobal x Nsnr x M x Nlocal] matrix
    end

    %%% Local detection parameters

    % properties (SetAccess = private, GetAccess = public)
    properties (SetAccess = public, GetAccess = public)
        localPFA (1, :) cell = {0.1} % [1 x M] cell of [M x 1 x 1 x Nlocal] matrices
        localPFAsimulation (1, :) cell = {} % [1 x M] cell of [M x 1 x 1 x Nlocal] matrices
    end

    properties (Dependent)
        localThreshold (1, :) cell % [1 x M] cell of [M x 1 x 1 x Nlocal] matrices
        localPD (1, :) cell % [1 x M] cell of [M x Nsnr x 1 x Nlocal] matrices
    end

    properties (SetAccess = private, GetAccess = public)
        localPDsimulation (1, :) cell = {} % [1 x M] cell of [M x Nsnr x 1 x Nlocal] matrices
    end

    %%% Random variable statistics

    properties (Dependent)
        numberOfActiveSensorsExpectationUnderNoise double {mustBeNonnegative} % [1 x 1 x M x Nlocal] matrix
        numberOfActiveSensorsStandardDeviationUnderNoise double {mustBeNonnegative} % [1 x 1 x M x Nlocal] matrix
        numberOfActiveSensorsExpectationUnderSignal double {mustBeNonnegative} % [1 x Nsnr x M x Nlocal] matrix
        numberOfActiveSensorsStandardDeviationUnderSignal double {mustBeNonnegative} % [1 x Nsnr x M x Nlocal] matrix
    end

    properties (SetAccess = private, GetAccess = public)
        numberOfActiveSensorsExpectationUnderNoiseSimulation double {mustBeNonnegative} % [1 x 1 x M x Nlocal] matrix
        numberOfActiveSensorsStandardDeviationUnderNoiseSimulation double {mustBeNonnegative} % [1 x 1 x M x Nlocal] matrix
        numberOfActiveSensorsExpectationUnderSignalSimulation double {mustBeNonnegative} % [1 x Nsnr x M x Nlocal] matrix
        numberOfActiveSensorsStandardDeviationUnderSignalSimulation double {mustBeNonnegative} % [1 x Nsnr x M x Nlocal] matrix
    end

    %%% Data generation

    properties (Dependent, Hidden)
        noise
        signal
    end

    methods
        function obj = detector(options)
            %detector Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                options.globalPFA (:, 1) double {mustBeNonnegative, mustBeInRange(options.globalPFA, 0, 1)} = 1e-6
                options.numberOfSensors (1, :) double {mustBeNonnegative, mustBeInteger} = 10
                options.SNR_input_dB = {4} % [1 x M] cell of [M x Nsnr] matrices
                options.localPFA = {0.1} % [1 x M] cell of [M x Nlocal] matrices
            end
            obj.globalPFA = options.globalPFA;
            obj.numberOfSensors = options.numberOfSensors;
            obj.numberOfRequiredDetections = ceil((obj.numberOfSensors + 1)/2);
            numberOfScansPFAglobal = length(obj.globalPFA);
            numberOfScans = length(obj.numberOfSensors);
            obj.setSNR(options.SNR_input_dB)
            numberOfScansSNR = unique((cellfun(@(c) size(c, 2), obj.SNR_input_dB)));
            obj.localPFA = cell(1, numberOfScans); % [1 x M] cell
            if ~iscell(options.localPFA)
                mustBeNonnegative(options.localPFA);
                mustBeInRange(options.localPFA, 0, 1)
                options.localPFA = {options.localPFA};
            else
                cellfun(@(c) mustBeNonnegative(c), options.localPFA);
                cellfun(@(c) mustBeInRange(c, 0, 1), options.localPFA);
                mustBeVector(options.localPFA);
            end
            if isscalar(options.localPFA)
                options.localPFA = cell2mat(options.localPFA);
                for scanID = 1 : numberOfScans
                    if size(options.localPFA, 1) == 1 % same SNR across all sensors
                        obj.localPFA{scanID} = repmat(options.localPFA, obj.numberOfSensors(scanID), 1);
                    else
                        obj.localPFA{scanID} = options.localPFA(1 : obj.numberOfSensors(scanID), :);
                    end
                end
            else
                mustBeScalarOrEmpty(unique((cellfun(@(c) size(c, 2), options.localPFA)))); % to ensure same Nlocal
                for scanID = 1 : numberOfScans
                    localPFAmatrix = options.localPFA{scanID};
                    firstSizeLocalPFAmatrix = size(localPFAmatrix, 1);
                    assert(firstSizeLocalPFAmatrix == obj.numberOfSensors(scanID) || firstSizeLocalPFAmatrix == 1);
                    if firstSizeLocalPFAmatrix == 1 % same SNR across all sensors
                        localPFAmatrix = repmat(localPFAmatrix, obj.numberOfSensors(scanID), 1);
                    end
                    obj.localPFA{scanID} = localPFAmatrix;
                end
            end
            obj.localPFA = cellfun(@(c) permute(c, [1 3 4 2]), obj.localPFA, 'UniformOutput', false);
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            obj.globalPFAsimulation = zeros(numberOfScansPFAglobal, numberOfScansSNR, numberOfScans, numberOfScansLocalPFA);
            obj.globalPDsimulation = zeros(numberOfScansPFAglobal, numberOfScansSNR, numberOfScans, numberOfScansLocalPFA);
            obj.localPFAsimulation = cell(1, numberOfScans);
            obj.localPDsimulation = cell(1, numberOfScans);
            for scanID = 1 : numberOfScans
                obj.localPFAsimulation{scanID} = zeros(numberOfScans, 1, 1, numberOfScansLocalPFA); % [1 x M] cell of [M x 1 x 1 x Nlocal] matrices
                obj.localPDsimulation{scanID} = zeros(numberOfScans, numberOfScansSNR, 1, numberOfScansLocalPFA); % [1 x M] cell of [M x Nsnr x 1 x Nlocal] matrices
            end
            obj.numberOfActiveSensorsExpectationUnderNoiseSimulation = zeros(1, 1, numberOfScans, numberOfScansLocalPFA); % [1 x 1 x M x Nlocal] matrix
            obj.numberOfActiveSensorsStandardDeviationUnderNoiseSimulation = zeros(1, 1, numberOfScans, numberOfScansLocalPFA); % [1 x 1 x M x Nlocal] matrix
            obj.numberOfActiveSensorsExpectationUnderSignalSimulation = zeros(1, numberOfScansSNR, numberOfScans, numberOfScansLocalPFA); % [1 x Nsnr x M x Nlocal] matrix
            obj.numberOfActiveSensorsStandardDeviationUnderSignalSimulation = zeros(1, numberOfScansSNR, numberOfScans, numberOfScansLocalPFA); % [1 x Nsnr x M x Nlocal] matrix
        end

        %%% Fusion rule

        function K = get.numberOfRequiredDetections(obj)
            % [1 x M] vector
            switch obj.binaryDetectionRule
                case "and"
                    K = obj.numberOfSensors;
                case "or"
                    K = ones(1, length(obj.numberOfSensors));
                case "majority"
                    K = ceil((obj.numberOfSensors + 1)/2);
                case "notSpecified"
                    K = obj.numberOfSensors;
                case "userDefined"
                    K = obj.numberOfRequiredDetections;
            end
        end

        %%% Fusion parameters

        function w = get.fusionWeights(obj)
            % [1 x M] cell of [M x Nsnr x 1 x Nlocal] matrices
            numberOfScansSNR = unique(cellfun(@(c) size(c, 2), obj.SNR_input_dB));
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            w = cell(1, numberOfScans);
            switch obj.globalFusionRule
                case "WSLC"
                    for scanID = 1 : numberOfScans
                        switch obj.signalAmplitudeModel
                            case "decorrelatedExponential"
                                snr = 10.^(.1*obj.SNR_input_dB{scanID});
                                w{scanID} = repmat(snr./(1 + snr), [1 1 1 numberOfScansLocalPFA]);
                            case "correlatedExponential"
                                w{scanID} = repmat(10.^(.1*obj.SNR_input_dB{scanID}), [1 1 1 numberOfScansLocalPFA]);
                                w{scanID} = w{scanID}./sum(w{scanID}, 1);
                        end
                    end
                case "CVBC"
                    for scanID = 1 : numberOfScans
                        pfaLocal = obj.localPFA{scanID};
                        pdLocal = obj.localPD{scanID}; % requires knowledge of SNR
                        w{scanID} = log(pdLocal.*(1 - pfaLocal)./pfaLocal./(1 - pdLocal));
                        % lets make the sum of weights equal to #sensors so that CVBC and BC has the same maximum test statistics value
                        w{scanID} = w{scanID}./sum(w{scanID}, 1)*obj.numberOfSensors(scanID);
                        w{scanID}(isnan(w{scanID})) = 1;
                    end
                case "WLLC"
                    for scanID = 1 : numberOfScans
                        w{scanID} = repmat(10.^(.05*obj.SNR_input_dB{scanID}), [1 1 1 numberOfScansLocalPFA]);
                    end
                case "MRC"
                    for scanID = 1 : numberOfScans
                        w{scanID} = repmat(10.^(.05*obj.SNR_input_dB{scanID}), [1 1 1 numberOfScansLocalPFA]);
                    end
                otherwise
                    for scanID = 1 : numberOfScans
                        w{scanID} = ones(obj.numberOfSensors(scanID), numberOfScansSNR, 1, numberOfScansLocalPFA);
                    end
            end
        end

        function Tglobal = get.globalThreshold(obj)
            % [NpfaGlobal x Nsnr x M x Nlocal]
            simulateThreshold = {"WSLC", "SLC", "CVBC", "BC"};
            simulateThreshold = {"WSLC", "SLC", "CVBC"};
            simulateThreshold = {};

            averageLocalPFA = permute(cell2mat(cellfun(@(c) mean(c, 1), obj.localPFA, 'UniformOutput', false).'), [1 4 2 3]);
            uniqueLocalPFA = unique(averageLocalPFA, 'rows', 'stable');
            switch obj.globalFusionRule
                case {"EGC", "MRC"} % coherent processing
                    if ~isscalar(uniqueLocalPFA) || uniqueLocalPFA ~= 1
                        % for localPFA is not equal to 1 we do not have a closed form expression for the distribution of noise
                        Tglobal = obj.globalThresholdSimulation;
                        return;
                    end
                case [{"LLC", "WLLC"}, simulateThreshold]
                    % sum of i.i.d. Rayleigh random variables do not have closed form expression
                    Tglobal = obj.globalThresholdSimulation;
                    return;
            end
            Tlocal = obj.localThreshold;
            weights = obj.fusionWeights;
            switch obj.globalFusionRule
                case {"WSLC", "CVBC", "WLLC", "MRC"} % weighted
                    numberOfScansSNR = unique(cellfun(@(c) size(c, 2), obj.SNR_input_dB));
                otherwise
                    numberOfScansSNR = 1;
            end
            numberOfScansLocalPFA = unique(cellfun(@(c) size(c, 4), obj.localPFA));
            numberOfScans = length(obj.numberOfSensors);
            Tglobal = zeros(length(obj.globalPFA), numberOfScansSNR, length(obj.numberOfSensors), numberOfScansLocalPFA);
            for scanID = 1 : numberOfScans
                M = obj.numberOfSensors(scanID);
                m = obj.numberOfRequiredDetections(scanID);
                if strcmp(obj.globalFusionRule, "CVBC")
                    allSubsets = dec2bin(0 : 2^M - 1) - '0';
                    numberOfSubsets = size(allSubsets, 1);
                    pmf = zeros(numberOfSubsets, 1);
                    weightSet = zeros(numberOfSubsets, 1);
                end
                for PFAID = 1 : length(obj.globalPFA)
                    pfaGlobal = obj.globalPFA(PFAID);
                    for localPFAID = 1 : numberOfScansLocalPFA
                        pfaLocal = obj.localPFA{scanID}(:, :, :, localPFAID);
                        lambda = Tlocal{scanID}(:, :, :, localPFAID);
                        for snrID = 1 : numberOfScansSNR
                            w = weights{scanID}(:, snrID, 1, localPFAID);
                            F = obj.totalCDF(M, m, lambda, pfaLocal, w);
                            switch obj.globalFusionRule
                            %%% LRT: widely seperated receivers | spatially i.i.d circularly symmetric complex gaussian signal model
                                % SLC: Square Law Combining (LRT under equal SNR)
                                % WSLC: Weighted Sqaure Law Combining (LRT under different SNR)
                        
                            %%% Binary
                                % BC: Binary Combining (LRT under equal SNR)
                                % CVBC: Chair-Varshney Binary Combining (LRT under different SNR)
                        
                            %%% Other
                                % LLC: Linear Law Combining
                                % WLLC: Weighted Linear Law Combining
                                % NSLC: Normalized Square Law Combining
                                % SC: Selective Combining (Partial Selective Combining K = 1)
                                % PSC: Partial Selective Combining (Combine largest K samples)

                            %%% LRT: phase synchronous closely spaced receivers | spatially fully correlated circularly symmetric complex gaussian signal model
                                % EGC: Equal Gain Combining
                                % MRC: Matched Filter

                                case {"SLC", "WSLC", "NSLC", "SC", "PSC"}
                                    thresholdNSC = gammaincinv(pfaGlobal, M, 'upper'); % unweighted non-selective power combiner threshold, conventional
                                    thresholdSearchSpace = [0, sum(w.*lambda) + 2*thresholdNSC];
                                case {"BC", "CVBC"} %%% Binary
                                    thresholdSearchSpace = [0, sum(w)];
                            end
                            switch obj.globalFusionRule
                                case {"SLC", "WSLC", "NSLC", "PSC"} % continous float threshold
                                    Tglobal(PFAID, snrID, scanID, localPFAID) = fzero(@(t) F(t) - pfaGlobal, thresholdSearchSpace);
                                case {"BC", "CVBC"} % discrete thresholds
                                    if all(abs(w - 1) < 1e-12, 'all') || any(strcmp(obj.globalFusionRule, "BC")) % discrete uniform integer thresholds
                                        % - 0.5 prevents boundary problems
                                        switch obj.binaryDetectionRule
                                            case "notSpecified"
                                                pfaGlobalK = zeros(1, M + 2);
                                                for K = -1 : M
                                                    pfaGlobalK(K + 2) = F(K, pfaLocal);
                                                end
                                                switch obj.binaryDetectionPFAtype
                                                    case "fixedGlobalPFA"
                                                        % since global PFA is fixed with local PFA already
                                                        [~, threshold] = min(abs(pfaGlobalK - pfaGlobal));
                                                        threshold = threshold - 1.5;
                                                    case {"fixedLocalPFA", "fixedGlobal|LocalPFA"}
                                                        threshold = find(pfaGlobalK < pfaGlobal, 1, 'first') - 1.5;
                                                        if strcmp(obj.binaryDetectionPFAtype, "fixedGlobal|LocalPFA")
                                                            % T has less global PFA than desired
                                                            % T - 1 has higher global PFA than desired
                                                            pfaGlobalLower = pfaGlobalK(ceil(threshold + 1));
                                                            pfaGlobalUpper = pfaGlobalK(max(ceil(threshold), 1));
                                                            if pfaGlobalUpper ~= pfaGlobalLower
                                                                obj.globalRandomizationProbability(PFAID, snrID, scanID, localPFAID) = (pfaGlobal - pfaGlobalLower)/(pfaGlobalUpper - pfaGlobalLower);
                                                            else
                                                                obj.globalRandomizationProbability(PFAID, snrID, scanID, localPFAID) = 0;
                                                            end
                                                        end
                                                    case "timeSharing"
                                                        [~, qLocal] = binaryProbability(M, pfaLocal(1), pfaGlobal, -inf);
                                                        threshold = find(~isnan(qLocal)) + 0.5;
                                                        obj.globalRandomizationProbability(PFAID, snrID, scanID, localPFAID) = qLocal(~isnan(qLocal));
                                                end
                                            otherwise
                                                threshold = m - 0.5;
                                        end
                                    else % discrete but nonuniform float thresholds
                                        % for subSetID = 1 : numberOfSubsets
                                        %     indices = logical(allSubsets(subSetID, :));
                                        %     pmf(subSetID) = prod(pfaLocal(indices))*prod(1 - pfaLocal(~indices));
                                        %     weightSet(subSetID) = sum(w(indices));
                                        % end
                                        % [weightSet, sortIndices] = sort(weightSet);
                                        % ccdf = 1 - cumsum(pmf(sortIndices));

                                        ccdf = CVBCCCDF(pfaLocal);
                                        thresholdIndex = find(ccdf - eps < pfaGlobal, 1, 'first');
                                        if thresholdIndex ~= numberOfSubsets
                                            threshold = mean(weightSet([thresholdIndex, thresholdIndex + 1]));
                                        else
                                            threshold = weightSet(thresholdIndex) + 0.5;
                                        end

                                        switch obj.binaryDetectionPFAtype
                                            case "fixedGlobal|LocalPFA"
                                                % T has less global PFA than desired
                                                % T - 1 has higher global PFA than desired
                                                pfaGlobalLower = ccdf(thresholdIndex);
                                                pfaGlobalUpper = ccdf(max(ceil(thresholdIndex - 1), 1));
                                                if pfaGlobalUpper ~= pfaGlobalLower
                                                    obj.globalRandomizationProbability(PFAID, snrID, scanID, localPFAID) = max((pfaGlobal - pfaGlobalLower)/(pfaGlobalUpper - pfaGlobalLower), 0);
                                                else
                                                    obj.globalRandomizationProbability(PFAID, snrID, scanID, localPFAID) = 0;
                                                end
                                                neighbourThresholdIndex = thresholdIndex - 1;
                                                if neighbourThresholdIndex ~= 0
                                                    neighbourThreshold = mean(weightSet([neighbourThresholdIndex, neighbourThresholdIndex + 1]));
                                                else
                                                    neighbourThreshold = 0;
                                                end
                                            case "timeSharing"
                                                rulesPfaLocal = zeros(1, numberOfSubsets);
                                                for subSetID = 1 : numberOfSubsets
                                                    rulesPfaLocal(subSetID) = fzero(@(pfaLoc) CVBCCCDF(pfaLoc) - pfaGlobal, [0 1]);
                                                end
                                                pRules = diag(CVBCCCDF(M, 1 : numberOfSubsets, rulesPfaLocal, 10^(rho/10)));
                                                qLocal = zeros(M - 1, length(pfaLocal));
                                                for i = 2 : M
                                                    qLocal(i - 1, :) = (pfaLocal - rulesPfaLocal(i))./(rulesPfaLocal(i - 1) - rulesPfaLocal(i));
                                                end
                                                qLocal(qLocal > 1 | qLocal < 0) = nan;
                                        end
                                        obj.globalNeighbourThresholdCVBC(PFAID, snrID, scanID, localPFAID) = neighbourThreshold;
                                    end
                                    Tglobal(PFAID, snrID, scanID, localPFAID) = threshold;
                                case "SC" % inverse exists
                                    Tglobal(PFAID, snrID, scanID, localPFAID) = max(unique(lambda), -log(1 - (1 - pfaGlobal).^(1/M)));
                                case {"EGC", "MRC"}
                                    Tglobal(PFAID, snrID, scanID, localPFAID) = -sum(w.^2).*log(pfaGlobal);
                            end
                        end
                    end
                end
            end
            function ccdfOut = CVBCCCDF(localPFA)
                for subSetID = 1 : numberOfSubsets
                    indices = logical(allSubsets(subSetID, :));
                    pmf(subSetID) = prod(localPFA(indices))*prod(1 - localPFA(~indices));
                    weightSet(subSetID) = sum(w(indices));
                end
                [weightSet, sortIndices] = sort(weightSet);
                ccdfOut = 1 - cumsum(pmf(sortIndices));
            end
        end
        
        function pfa = get.globalPFAanalytical(obj)
            % [NpfaGlobal x Nsnr x M x Nlocal]
            Tglobal = obj.globalThreshold;
            if isempty(Tglobal)
                pfa = [];
                fprintf('simulate global threshold\n');
                return;
            end
            Tlocal = obj.localThreshold;
            pFAlocal = obj.localPFA;
            weights = obj.fusionWeights;
            numberOfScansSNR = unique(cellfun(@(c) size(c, 2), obj.SNR_input_dB));
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            pfa = zeros(length(obj.globalPFA), numberOfScansSNR, length(obj.numberOfSensors), numberOfScansLocalPFA);
            for scanID = 1 : numberOfScans
                M = obj.numberOfSensors(scanID);
                m = obj.numberOfRequiredDetections(scanID);
                for PFAID = 1 : length(obj.globalPFA)
                    rulesPfaLocal = zeros(1, M);
                    for k = 1 : M
                        rulesPfaLocal(k) = fzero(@(pfaLoc) binaryCCDF(M, k, pfaLoc, 0) - obj.globalPFA(PFAID), [0 1]);
                    end
                    pRules = diag(binaryCCDF(M, 1 : M, rulesPfaLocal, 0));
                    for localPFAID = 1 : numberOfScansLocalPFA
                        if any(strcmp(obj.globalFusionRule, "CVBC")) && strcmp(obj.binaryDetectionPFAtype, "fixedGlobal|LocalPFA")
                            gammaNeighbour = obj.globalNeighbourThresholdCVBC(PFAID, 1, scanID, localPFAID);
                        end
                        lambda = Tlocal{scanID}(:, :, :, localPFAID);
                        pfaLocal = pFAlocal{scanID}(:, :, :, localPFAID);
                        for snrID = 1 : numberOfScansSNR
                            gamma = Tglobal(PFAID, snrID, scanID, localPFAID);
                            w = weights{scanID}(:, snrID, 1, localPFAID);
                            F = obj.totalCDF(M, m, lambda, pfaLocal, w);
                            pRulesOR = binaryCCDF(M, 1, pfaLocal(1), 0);
                            pRulesAND = binaryCCDF(M, M, pfaLocal(1), 0);
                            if any(strcmp(obj.globalFusionRule, ["BC", "CVBC"]))
                                if all(abs(w - 1) < 1e-12, 'all') % discrete uniform integer thresholds
                                    switch obj.binaryDetectionPFAtype
                                        case "fixedGlobal|LocalPFA"
                                            q = obj.globalRandomizationProbability(PFAID, 1, scanID, localPFAID);
                                            pfa(PFAID, snrID, scanID, localPFAID) = (1 - q)*F(gamma, pfaLocal) + q*F(gamma - 1, pfaLocal);
                                        case "timeSharing"
                                            q = obj.globalRandomizationProbability(PFAID, 1, scanID, localPFAID);
                                            if gamma > M + 1
                                                pfa(PFAID, snrID, scanID, localPFAID) = q*0 + q*pRulesAND;
                                            elseif gamma < 2
                                                pfa(PFAID, snrID, scanID, localPFAID) = q*1 + (1 - q).*pRulesOR;
                                            else
                                                pfa(PFAID, snrID, scanID, localPFAID) = (1 - q)*pRules(floor(gamma)) + q*pRules(floor(gamma) - 1);
                                            end
                                        otherwise
                                            pfa(PFAID, snrID, scanID, localPFAID) = F(gamma, pfaLocal);
                                    end
                                else
                                    if ~strcmp(obj.binaryDetectionPFAtype, "fixedGlobal|LocalPFA")
                                        pfa(PFAID, snrID, scanID, localPFAID) = F(gamma, pfaLocal);
                                    else
                                        q = obj.globalRandomizationProbability(PFAID, snrID, scanID, localPFAID);
                                        pfa(PFAID, snrID, scanID, localPFAID) = (1 - q)*F(gamma, pfaLocal) + q*F(gammaNeighbour, pfaLocal);
                                    end
                                end
                            else
                                pfa(PFAID, snrID, scanID, localPFAID) = F(gamma);
                            end
                        end
                        % evaluating @(gamma - 1) is important for binary fusion with integer threshold
                    end
                end
            end
        end

        function ROC = get.globalPD(obj)
            % [NpfaGlobal x Nsnr x M x Nlocal]
            Tglobal = obj.globalThreshold;
            if isempty(Tglobal)
                ROC = [];
                fprintf('simulate global threshold\n');
                return;
            end
            Tlocal = obj.localThreshold;
            pFAlocal = obj.localPFA;
            pDlocal = obj.localPD;
            weights = obj.fusionWeights;
            numberOfScansSNR = unique(cellfun(@(c) size(c, 2), obj.SNR_input_dB));
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            ROC = zeros(length(obj.globalPFA), numberOfScansSNR, numberOfScans, numberOfScansLocalPFA);
            for scanID = 1 : numberOfScans
                M = obj.numberOfSensors(scanID);
                m = obj.numberOfRequiredDetections(scanID);
                snr = 10.^(.1*obj.SNR_input_dB{scanID});
                for PFAID = 1 : length(obj.globalPFA)
                    rulesPfaLocal = zeros(1, M);
                    for k = 1 : M
                        rulesPfaLocal(k) = fzero(@(pfaLoc) binaryCCDF(M, k, pfaLoc, 0) - obj.globalPFA(PFAID), [0 1]);
                    end
                    for localPFAID = 1 : numberOfScansLocalPFA
                        if any(strcmp(obj.globalFusionRule, "CVBC")) && strcmp(obj.binaryDetectionPFAtype, "fixedGlobal|LocalPFA")
                            gammaNeighbour = obj.globalNeighbourThresholdCVBC(PFAID, 1, scanID, localPFAID);
                        end
                        lambda = Tlocal{scanID}(:, :, :, localPFAID);
                        pfaLocal = pFAlocal{scanID}(:, :, :, localPFAID);
                        for snrID = 1 : numberOfScansSNR
                            if size(Tglobal, 2) == 1
                                gamma = Tglobal(PFAID, :, scanID, localPFAID);
                            else
                                gamma = Tglobal(PFAID, snrID, scanID, localPFAID);
                            end
                            pdLocal = pDlocal{scanID}(:, snrID, 1, localPFAID);
                            w = weights{scanID}(:, snrID, 1, localPFAID);
                            F = obj.totalCDF(M, m, lambda, pdLocal, w, snr(:, snrID));
                            pRules = diag(binaryCCDF(M, 1 : M, rulesPfaLocal, snr(1, snrID)));
                            pRulesOR = binaryCCDF(M, 1, pfaLocal(1), snr(1, snrID));
                            pRulesAND = binaryCCDF(M, M, pfaLocal(1), snr(1, snrID));
                            if any(strcmp(obj.globalFusionRule, ["BC", "CVBC"]))
                                if all(abs(w - 1) < 1e-12, 'all') % discrete uniform integer thresholds
                                    switch obj.binaryDetectionPFAtype
                                        case "fixedGlobal|LocalPFA"
                                            q = obj.globalRandomizationProbability(PFAID, 1, scanID, localPFAID);
                                            ROC(PFAID, snrID, scanID, localPFAID) = (1 - q)*F(gamma, pdLocal) + q*F(gamma - 1, pdLocal);
                                        case "timeSharing"
                                            q = obj.globalRandomizationProbability(PFAID, 1, scanID, localPFAID);
                                            if gamma > M + 1
                                                ROC(PFAID, snrID, scanID, localPFAID) = q*0 + q*pRulesAND;
                                            elseif gamma < 2
                                                ROC(PFAID, snrID, scanID, localPFAID) = q*1 + (1 - q).*pRulesOR;
                                            else
                                                ROC(PFAID, snrID, scanID, localPFAID) = (1 - q)*pRules(floor(gamma)) + q*pRules(floor(gamma) - 1);
                                            end
                                        otherwise
                                            ROC(PFAID, snrID, scanID, localPFAID) = F(gamma, pdLocal);
                                    end
                                else
                                    if ~strcmp(obj.binaryDetectionPFAtype, "fixedGlobal|LocalPFA")
                                        ROC(PFAID, snrID, scanID, localPFAID) = F(gamma, pdLocal);
                                    else
                                        q = obj.globalRandomizationProbability(PFAID, snrID, scanID, localPFAID);
                                        ROC(PFAID, snrID, scanID, localPFAID) = (1 - q)*F(gamma, pdLocal) + q*F(gammaNeighbour, pdLocal);
                                    end
                                end
                            else
                                ROC(PFAID, snrID, scanID, localPFAID) = F(gamma);
                            end
                            if ~mod(snrID, 10000)
                                fprintf('%d/%d\n', snrID, numberOfScansSNR);
                            end
                            % evaluating @(gamma - 1) is important for binary fusion with integer threshold
                        end
                        if ~mod(localPFAID, 10)
                            fprintf('%d/%d\n', localPFAID, numberOfScansLocalPFA);
                        end
                    end
                end
            end
        end

        %%% Local detection parameters

        function T = get.localThreshold(obj)
            % [1 x M] cell of [M x 1 x 1 x Nlocal] matrices
            T = cellfun(@(c) -log(c), obj.localPFA, 'UniformOutput', false);
        end

        function ROC = get.localPD(obj)
            % [1 x M] cell of [M x Nsnr x 1 x Nlocal] matrices
            switch obj.signalAmplitudeModel
                case "deterministic"
                    theta = cellfun(@(c) 10.^(.05*c), obj.SNR_input_dB, 'UniformOutput', false); % [1 x M] cell of [M x Nsnr] matrices
                    ROC = cellfun(@(c, t) marcumq(repmat(sqrt(2)*t, [1, 1, 1, size(c, 4)]), repmat(sqrt(-2*log(c)), [1, size(t, 2)]), 1), obj.localPFA, theta, 'UniformOutput', false);
                otherwise
                    theta = cellfun(@(c) 1./(1 + 10.^(.1*c)), obj.SNR_input_dB, 'UniformOutput', false); % [1 x M] cell of [M x Nsnr] matrices
                    ROC = cellfun(@(c, t) c.^t, obj.localPFA, theta, 'UniformOutput', false);
            end
        end

        %%% Random variable statistics

        function M = get.numberOfActiveSensorsExpectationUnderNoise(obj)
            % [1 x 1 x M x Nlocal] matrices
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            M = reshape(cell2mat(cellfun(@(c) sum(c, 1), obj.localPFA, 'UniformOutput', false)), [1, 1, numberOfScans, numberOfScansLocalPFA]);
        end

        function Mstd = get.numberOfActiveSensorsStandardDeviationUnderNoise(obj)
            % [1 x 1 x M x Nlocal] matrices
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            Mstd = reshape(cell2mat(cellfun(@(c) sqrt(sum(c.*(1 - c), 1)), obj.localPFA, 'UniformOutput', false)), [1, 1, numberOfScans, numberOfScansLocalPFA]);
        end

        function M = get.numberOfActiveSensorsExpectationUnderSignal(obj)
            % [1 x Nsnr x M x Nlocal] matrices
            numberOfScansSNR = unique(cellfun(@(c) size(c, 2), obj.SNR_input_dB));
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            M = reshape(cell2mat(cellfun(@(c) sum(c, 1), obj.localPD, 'UniformOutput', false)), [1, numberOfScansSNR, numberOfScans, numberOfScansLocalPFA]);
        end

        function Mstd = get.numberOfActiveSensorsStandardDeviationUnderSignal(obj)
            % [1 x Nsnr x M x Nlocal] matrices
            numberOfScansSNR = unique(cellfun(@(c) size(c, 2), obj.SNR_input_dB));
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            Mstd = reshape(cell2mat(cellfun(@(c) sqrt(sum(c.*(1 - c), 1)), obj.localPD, 'UniformOutput', false)), [1, numberOfScansSNR, numberOfScans, numberOfScansLocalPFA]);
        end

        %%% Data generation

        function n = get.noise(obj)
            % [M x 1 x Nmc]
            % independent and identically distributed across samples
            % zero mean, unit variance circularly symmetric complex gaussian noise
            n = (randn(max(obj.numberOfSensors), 1, obj.numberOfTrials, 'single') + 1j*randn(max(obj.numberOfSensors), 1, obj.numberOfTrials, 'single'))/sqrt(2);
        end

        function s = get.signal(obj)
            % [1 x M] cell of [M x Nsnr] matrices
            switch obj.signalPhaseModel
                case "decorrelatedUniform"
                    switch obj.signalAmplitudeModel
                        case "decorrelatedExponential"
                            % swerling-2 like spatially independent unit power signal
                            % meaningful for widely seperated receivers, "FLUCTUATING target"
                            complexUnitPowerSignal = (randn(max(obj.numberOfSensors), 1, obj.numberOfTrials, 'single') + 1j*randn(max(obj.numberOfSensors), 1, obj.numberOfTrials, 'single'))/sqrt(2);
                        case "correlatedExponential"
                            % swerling-1 like spatially fully correlated unit power signal
                            % meaningful for "phase ASYNCHRONOUS" closely spaced receivers, "FLUCTUATING target"
                            complexUnitPowerSignal = single(raylrnd(1/sqrt(2), 1, 1, obj.numberOfTrials)).*exp(1j*2*pi*rand(max(obj.numberOfSensors), 1, obj.numberOfTrials, 'single'));
                        case "deterministic"
                            % swerling-0 like spatially non-fluctuating unit power signal
                            % meaningful for "phase ASYNCHRONOUS" closely spaced receivers, "NON-FLUCTUATING target"
                            complexUnitPowerSignal = exp(1j*2*pi*rand(max(obj.numberOfSensors), 1, obj.numberOfTrials, 'single'));
                    end
                case "correlatedUniform"
                    % Coherent processing is possible
                    switch obj.signalAmplitudeModel
                        case "decorrelatedExponential"
                            error('not meaningful');
                        case "correlatedExponential"
                            % meaningful for "phase SYNCHRONIZED" closely spaced receivers, "FLUCTUATING target"
                            complexUnitPowerSignal = repmat((randn(1, 1, obj.numberOfTrials, 'single') + 1j*randn(1, 1, obj.numberOfTrials, 'single'))/sqrt(2), max(obj.numberOfSensors), 1);
                        case "deterministic"
                            % meaningful for "phase SYNCHRONIZED" closely spaced receivers, "NON-FLUCTUATING target"
                            complexUnitPowerSignal = exp(1j*2*pi*rand(1, 1, obj.numberOfTrials, 'single'));
                    end
            end
            numberOfScans = length(obj.numberOfSensors);
            s = cell(1, numberOfScans);
            for scanID = 1 : numberOfScans
                s{scanID} = 10.^(.05*obj.SNR_input_dB{scanID}).*complexUnitPowerSignal(1 : obj.numberOfSensors(scanID), 1, :);
            end
        end

        %%% set functions

        function setalgorithm(obj, options)
            arguments
                obj
                options.globalFusionRule (1, 1) string {mustBeMember(options.globalFusionRule, ["SLC", "WSLC", "BC", "CVBC", "LLC", "WLLC", "NSLC", "SC", "PSC", "EGC", "MRC"])} = obj.globalFusionRule
                %%% LRT: widely seperated receivers | spatially i.i.d circularly symmetric complex gaussian signal model
                    % SLC: Square Law Combining (LRT under equal SNR)
                    % WSLC: Weighted Sqaure Law Combining (LRT under different SNR)
        
                %%% Binary
                    % BC: Binary Combining (LRT under equal SNR)
                    % CVBC: Chair-Varshney Binary Combining (LRT under different SNR)
        
                %%% Other
                    % LLC: Linear Law Combining
                    % WLLC: Weighted Linear Law Combining
                    % NSLC: Normalized Square Law Combining
                    % SC: Selective Combining (Partial Selective Combining K = 1)
                    % PSC: Partial Selective Combining (Combine largest K samples)

                %%% LRT: phase synchronous closely spaced receivers | spatially fully correlated circularly symmetric complex gaussian signal model
                    % EGC: Equal Gain Combining
                    % MRC: Matched Filter

                options.signalAmplitudeModel (1, 1) string {mustBeMember(options.signalAmplitudeModel, ["decorrelatedExponential", "correlatedExponential", "deterministic"])} = obj.signalAmplitudeModel
                options.signalPhaseModel (1, 1) string {mustBeMember(options.signalPhaseModel, ["decorrelatedUniform", "correlatedUniform"])} = obj.signalPhaseModel
                options.binaryDetectionRule (1, 1) string {mustBeMember(options.binaryDetectionRule, ["and", "or", "majority", "notSpecified", "userDefined"])} = obj.binaryDetectionRule
                options.binaryDetectionConstraint (1, 1) string {mustBeMember(options.binaryDetectionConstraint, ["fixedGlobalPFA", "fixedLocalPFA", "fixedGlobal|LocalPFA", "timeSharing"])} = obj.binaryDetectionPFAtype
                options.numberOfRequiredDetections (1, :) double {mustBeNonnegative, mustBeInteger} = obj.numberOfRequiredDetections % [1 x M]
            end
            obj.globalFusionRule = options.globalFusionRule;
            if strcmp(options.signalPhaseModel, "correlatedUniform") && strcmp(options.signalAmplitudeModel, "decorrelatedExponential")
                error('expecting correlated phase and decorrelated amplitude at the same time is not meaningful');
            end
            obj.signalAmplitudeModel = options.signalAmplitudeModel;
            obj.signalPhaseModel = options.signalPhaseModel;
            obj.binaryDetectionRule = options.binaryDetectionRule;
            obj.binaryDetectionPFAtype = options.binaryDetectionConstraint;
            if isscalar(options.numberOfRequiredDetections)
                obj.numberOfRequiredDetections = repmat(options.numberOfRequiredDetections, 1, length(obj.numberOfSensors));
            elseif length(options.numberOfRequiredDetections) == length(obj.numberOfSensors)
                obj.numberOfRequiredDetections = options.numberOfRequiredDetections;
            else
                error('Length of "#detections for binary Combining algorithm" vector must be 1 or %d', length(obj.numberOfSensors));
            end
            if strcmp(obj.globalFusionRule, "BC") && strcmp(obj.binaryDetectionPFAtype, "fixedGlobalPFA")
                numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
                numberOfScans = length(obj.numberOfSensors);
                for scanID = 1 : numberOfScans
                    M = obj.numberOfSensors(scanID);
                    for PFAID = 1 : length(obj.globalPFA)
                        pfaGlobal = obj.globalPFA(PFAID);
                        for localPFAID = 1 : numberOfScansLocalPFA
                            F = obj.totalCDF(M, 0, zeros(M, 1), ones(M, 1), ones(M, 1), zeros(M, 1), zeros(M, 1));
                            pfaLocalK = zeros(M, 1);
                            for K = 0 : M - 1
                                pfaLocalK(K + 1) = fzero(@(p) F(K, p*ones(M, 1)) - pfaGlobal, [0, 1]);
                            end
                            % local PFA monotonically increasing with K
                            switch obj.binaryDetectionRule
                                case "notSpecified"
                                    pfaLocal = obj.localPFA{scanID}(:, :, :, localPFAID);
                                    k = find(pfaLocalK <= pfaLocal, 1, 'first'); % find the greatest threshold possible
                                    if isempty(k)
                                        % to have fixed global PFA, local PFA must be greater than what is desired
                                        k = 1;
                                    end
                                otherwise
                                    k = min(obj.numberOfRequiredDetections(scanID), M);
                            end
                            obj.localPFA{scanID}(:, :, :, localPFAID) = pfaLocalK(k);
                        end
                    end
                end
            end
        end

        function setmontecarlo(obj, options)
            arguments
                obj
                options.numberOfTrials (1, 1) double {mustBeNonnegative, mustBeInteger} = obj.numberOfTrials
                options.numberOfLoops (1, 1) double {mustBeNonnegative, mustBeInteger} = obj.numberOfLoops
                options.seed (1, 1) double {mustBeNonnegative, mustBeInteger, mustBeInRange(options.seed, 0, 4294967295)} = obj.seed
            end
            obj.numberOfTrials = options.numberOfTrials;
            obj.numberOfLoops = options.numberOfLoops;
            obj.seed = options.seed;
        end

        function setbudgetconstraint(obj, options)
            arguments
                obj
                options.constraint (1, :) string {mustBeMember(options.constraint, ["dataRate", "transmittedPower"])} = ["dataRate", "transmittedPower"]
                options.numberOfTransmitters (1, :) double {mustBeNonnegative} = sqrt(obj.numberOfSensors)
            end
            if any(strcmp(options.constraint, "transmittedPower"))
                for scanID = 2 : length(obj.numberOfSensors)
                    % mustBeLessThanOrEqual(options.numberOfTransmitters(scanID), sqrt(obj.numberOfSensors(scanID)))
                    % obj.SNR_input_dB{scanID} = repmat(10*log10(10.^(.1*obj.SNR_input_dB{1})*options.numberOfTransmitters(1)/options.numberOfTransmitters(scanID)), obj.numberOfSensors(scanID), 1);
                    obj.SNR_input_dB{scanID} = repmat(10*log10(10.^(.1*obj.SNR_input_dB{1})*obj.numberOfSensors(1)/obj.numberOfSensors(scanID)), obj.numberOfSensors(scanID), 1);
                end
            end
            if any(strcmp(options.constraint, "dataRate"))
                for scanID = 2 : length(obj.numberOfSensors)
                    obj.localPFA{scanID} = repmat(obj.localPFA{1}*obj.numberOfSensors(1)/obj.numberOfSensors(scanID), obj.numberOfSensors(scanID), 1);
                end
            end
        end

        function setSNR(obj, SNR_input_dB)
            arguments
                obj
                SNR_input_dB = {4} % [1 x M] cell of [M x Nsnr] matrices
            end
            numberOfScans = length(obj.numberOfSensors);
            obj.SNR_input_dB = cell(1, numberOfScans); % [1 x M] cell
            if ~iscell(SNR_input_dB)
                SNR_input_dB = {SNR_input_dB};
            else
                mustBeVector(SNR_input_dB);
            end
            if isscalar(SNR_input_dB)
                SNR_input_dB = cell2mat(SNR_input_dB);
                for scanID = 1 : numberOfScans
                    if size(SNR_input_dB, 1) == 1 % same SNR across all sensors
                        obj.SNR_input_dB{scanID} = repmat(SNR_input_dB, obj.numberOfSensors(scanID), 1);
                    else
                        obj.SNR_input_dB{scanID} = SNR_input_dB(1 : obj.numberOfSensors(scanID), :);
                    end
                end
            else
                mustBeScalarOrEmpty(unique((cellfun(@(c) size(c, 2), SNR_input_dB)))); % to ensure same Nsnr
                for scanID = 1 : numberOfScans
                    SNRmatrix = SNR_input_dB{scanID};
                    firstSizeSNRmatrix = size(SNRmatrix, 1);
                    assert(firstSizeSNRmatrix == obj.numberOfSensors(scanID) || firstSizeSNRmatrix == 1);
                    if firstSizeSNRmatrix == 1 % same SNR across all sensors
                        SNRmatrix = repmat(SNRmatrix, obj.numberOfSensors(scanID), 1);
                    end
                    obj.SNR_input_dB{scanID} = SNRmatrix;
                end
            end
        end

        %%% simulation

        function simulate(obj, options)
            arguments
                obj
                options.simulationData (1, :) string {mustBeMember(options.simulationData, ["globalPFA", "globalPD", "globalThreshold"])} = ["globalPFA", "globalPD"]
                options.statistics (1, 1) logical {mustBeMember(options.statistics, [0, 1])} = false
                options.printStatus (1, 1) logical {mustBeMember(options.printStatus, [0, 1])} = true
            end
            rng(obj.seed);
            TglobalAll = obj.globalThreshold; % [NpfaGlobal x Nsnr x M x Nlocal]
            Tlocal = obj.localThreshold; % [1 x M] cell of [M x 1 x 1 x Nlocal] matrices
            numberOfSensorsMax = max(obj.numberOfSensors);
            numberOfScansSNR = unique(cellfun(@(c) size(c, 2), obj.SNR_input_dB));
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            if isempty(TglobalAll) || any(strcmp(options.simulationData, "globalThreshold"))
                thresholdSimulation = true;
                if isempty(obj.globalThresholdSimulation)
                    obj.globalThresholdSimulation = zeros(length(obj.globalPFA), numberOfScansSNR, numberOfScans, numberOfScansLocalPFA);
                end
            else
                thresholdSimulation = false;
            end
            for monteCarloID = 1 : obj.numberOfLoops
                if options.printStatus && ~mod(obj.numberOfLoops, 10)
                    fprintf('#loop = %d/%d\n', monteCarloID, obj.numberOfLoops);
                end
                n = obj.noise; % [M x 1 x Nmc] matrix
                if any(strcmp(options.simulationData, "globalPD"))
                    s = obj.signal; % [M x 1] cell of [1 x Nsnr x Nmc] matrices
                    x = cellfun(@(s) s + n(1 : size(s, 1), 1, :), s, 'UniformOutput', false); % [M x Nsnr x Nmc] matrix
                end
                if any(strcmp(obj.globalFusionRule, ["BC", "CVBC"])) && strcmp(obj.binaryDetectionPFAtype, "fixedGlobal|LocalPFA")
                    uniform01 = rand(1, 1, obj.numberOfTrials, 1, 'single');
                end
                for scanID = 1 : numberOfScans
                    if options.printStatus
                        fprintf('#sensors = %d/%d\n', obj.numberOfSensors(scanID), numberOfSensorsMax);
                    end
                    if thresholdSimulation
                        Tglobal = [];
                    else
                        Tglobal = TglobalAll(:, :, scanID, :);
                        if any(strcmp(obj.globalFusionRule, ["BC", "CVBC"])) && strcmp(obj.binaryDetectionPFAtype, "fixedGlobal|LocalPFA")
                            switch obj.globalFusionRule
                                case "BC"
                                    Tglobal = Tglobal - double(uniform01 < obj.globalRandomizationProbability(:, :, scanID, :));
                                case "CVBC"
                                    gammaNeighbour = obj.globalNeighbourThresholdCVBC(:, 1, scanID, :);
                                    useNeighbour = double(uniform01 < obj.globalRandomizationProbability(:, :, scanID, :));
                                    Tglobal = Tglobal + useNeighbour.*(gammaNeighbour - Tglobal);
                            end
                        end
                    end
                    % no signal is present
                    if any(strcmp(options.simulationData, "globalPFA")) || thresholdSimulation
                        localTestStatisticsH0 = abs(n(1 : obj.numberOfSensors(scanID), 1, :)).^2; % [M x 1 x Nmc]
                        indicatorFunction = localTestStatisticsH0 > Tlocal{scanID}; % [M x 1 x Nmc x Nlocal]
                        switch obj.globalFusionRule
                            case {"EGC", "MRC"} % coherent combiners
                                globalTestStatisticsH0 = fusioncenter(n(1 : obj.numberOfSensors(scanID), 1, :));
                            otherwise % power combiners
                                globalTestStatisticsH0 = fusioncenter(localTestStatisticsH0);
                        end
                        if monteCarloID == 1
                            obj.localPFAsimulation{scanID} = mean(indicatorFunction, 3)/obj.numberOfLoops; % [M x 1 x 1 x Nlocal]
                        else
                            obj.localPFAsimulation{scanID} = obj.localPFAsimulation{scanID} + mean(indicatorFunction, 3)/obj.numberOfLoops; % [M x 1 x 1 x Nlocal]
                        end
                        if thresholdSimulation
                            for localPFAID = 1 : numberOfScansLocalPFA
                                switch obj.globalFusionRule
                                    case {"WSLC", "CVBC", "MRC"} % weighted
                                        for snrID = 1 : numberOfScansSNR
                                            sortedTestStatisticsH0 = sort(globalTestStatisticsH0(:, snrID, :, localPFAID)); % [1 x Nsnr x Nmc x Nlocal]
                                            index = ceil((1 - obj.globalPFA)*obj.numberOfTrials);
                                            obj.globalThresholdSimulation(:, snrID, scanID, localPFAID) = sortedTestStatisticsH0(index);
                                            obj.globalThresholdSimulation(all(obj.globalThresholdSimulation(:, snrID, scanID, localPFAID) == 0, 2) & index > 1, snrID, scanID, localPFAID) = inf;
                                        end
                                    otherwise
                                        sortedTestStatisticsH0 = sort(globalTestStatisticsH0(:, :, :, localPFAID)); % [1 x 1 x Nmc x Nlocal]
                                        index = ceil((1 - obj.globalPFA)*obj.numberOfTrials);

                                        if obj.numberOfLoops ~= 1
                                            % not implemented
                                            % sortedTestStatisticsH0(end)
                                        end

                                        obj.globalThresholdSimulation(:, :, scanID, localPFAID) = repmat(squeeze(sortedTestStatisticsH0(index)), 1, numberOfScansSNR);
                                        obj.globalThresholdSimulation(all(obj.globalThresholdSimulation(:, :, scanID, localPFAID) == 0, 2) & index > 1, :, scanID, localPFAID) = inf;
                                end
                            end
                            Tglobal = obj.globalThresholdSimulation(:, :, scanID, :);
                        end
                        if monteCarloID == 1
                            obj.globalPFAsimulation(:, :, scanID, :) = mean(globalTestStatisticsH0 >= Tglobal, 3)/obj.numberOfLoops;
                        else
                            obj.globalPFAsimulation(:, :, scanID, :) = obj.globalPFAsimulation(:, :, scanID, :) + mean(globalTestStatisticsH0 >= Tglobal, 3)/obj.numberOfLoops; % [M x 1 x 1 x Nlocal]
                        end
                        % equality is important for binary fusion with integer threshold
                        if options.statistics && monteCarloID == obj.numberOfLoops
                            obj.numberOfActiveSensorsExpectationUnderNoiseSimulation(:, :, scanID, :) = sum(obj.localPFAsimulation{scanID}, 1); % [M x 1 x 1 x Nlocal] matrices
                            obj.numberOfActiveSensorsStandardDeviationUnderNoiseSimulation(:, :, scanID, :) = sqrt(sum(var(indicatorFunction, [], 3), 1)); % [M x 1 x 1 x Nlocal] matrices
                        end
                    end
                    % signal + noise is present
                    if any(strcmp(options.simulationData, "globalPD"))
                        localTestStatisticsH1 = abs(x{scanID}).^2; % [M x Nsr x Nmc]
                        indicatorFunction = localTestStatisticsH1 > Tlocal{scanID}; % [M x Nsr x Nmc x Nlocal]
                        switch obj.globalFusionRule
                            case {"EGC", "MRC"} % coherent combiners
                                globalTestStatisticsH1 = fusioncenter(x{scanID});
                            otherwise % power combiners
                                globalTestStatisticsH1 = fusioncenter(localTestStatisticsH1);
                        end
                        if monteCarloID == 1
                            obj.localPDsimulation{scanID} = mean(indicatorFunction, 3)/obj.numberOfLoops; % [M x Nsnr x 1 x Nlocal]
                            obj.globalPDsimulation(:, :, scanID, :) = mean(globalTestStatisticsH1 >= Tglobal, 3);
                        else
                            obj.localPDsimulation{scanID} = obj.localPDsimulation{scanID} + mean(indicatorFunction, 3)/obj.numberOfLoops; % [M x Nsnr x 1 x Nlocal]
                            obj.globalPDsimulation(:, :, scanID, :) = obj.globalPDsimulation(:, :, scanID, :) + mean(globalTestStatisticsH1 >= Tglobal, 3)/obj.numberOfLoops;
                        end
                        % equality is important for binary fusion with integer threshold
                        if options.statistics && monteCarloID == obj.numberOfLoops
                            obj.numberOfActiveSensorsExpectationUnderSignalSimulation(:, :, scanID, :) = sum(obj.localPDsimulation{scanID}, 1); % [M x Nsnr x 1 x Nlocal] matrices
                            obj.numberOfActiveSensorsStandardDeviationUnderSignalSimulation(:, :, scanID, :) = sqrt(sum(var(indicatorFunction, [], 3), 1)); % [M x Nsnr x 1 x Nlocal] matrices
                        end
                    end
                end
            end

            function globalTestStatistics = fusioncenter(localTestStatistics)
                localData = indicatorFunction.*localTestStatistics; % [M x Nsnr x Nmc x Nlocal]
                switch obj.globalFusionRule
                %%% LRT: widely seperated receivers | spatially i.i.d circularly symmetric complex gaussian signal model
                    case "SLC" % Square Law Combining (LRT under equal SNR)
                        globalTestStatistics = sum(localData, 1); % [1 x Nsnr x Nmc x Nlocal]
                    case "WSLC" % Weighted Sqaure Law Combining (LRT under different SNR)
                        weightError = false;
                        if weightError
                            % not fully implemented
                            rangedB = 20;
                            delta_dB = -rangedB/2 + rangedB*rand(obj.numberOfSensors(scanID), numberOfScansSNR, obj.numberOfTrials, numberOfScansLocalPFA);  % uniform in [-3, 3] dB
                            scaling_factor = 10.^(delta_dB/10);
                        end
                        globalTestStatistics = sum(obj.fusionWeights{scanID}.*localData, 1); % [1 x Nsnr x Nmc x Nlocal]
                %%% Binary
                    case "BC" % Binary Combining (LRT under equal SNR)
                        globalTestStatistics = sum(indicatorFunction, 1); % [1 x Nsnr x Nmc x Nlocal]
                    case "CVBC" % Chair-Varshney Binary Combining (LRT under different SNR)
                        globalTestStatistics = sum(obj.fusionWeights{scanID}.*indicatorFunction, 1); % [1 x Nsnr x Nmc x Nlocal]
                %%% Other
                    case "LLC" % Linear Law Combining
                        globalTestStatistics = sum(sqrt(localData), 1); % [1 x Nsnr x Nmc x Nlocal]
                    case "WLLC" % Weighted Linear Law Combining
                        globalTestStatistics = sum(obj.fusionWeights{scanID}.*sqrt(localData), 1); % [1 x Nsnr x Nmc x Nlocal]
                    case "NSLC" % Normalized Square Law Combining
                        globalTestStatistics = sum(localData, 1)./sum(indicatorFunction, 1); % [1 x Nsnr x Nmc x Nlocal]
                    case "SC" % Selective Combining (Partial Selective Combining K = 1)
                        globalTestStatistics = max(localTestStatistics, [], 1); % [1 x Nsnr x Nmc x Nlocal]
                    case "PSC" % Partial Selective Combining (Combine largest K samples)
                        localData = sort(localData, 1, 'descend');
                        globalTestStatistics = sum(localData(1 : min(obj.numberOfRequiredDetections(scanID), obj.numberOfSensors(scanID)), :, :, :), 1); % [1 x Nsnr x Nmc x Nlocal]
                %%% LRT: phase synchronous closely spaced receivers | spatially fully correlated circularly symmetric complex gaussian signal model
                    case "EGC" % Equal Gain Combining
                        globalTestStatistics = abs(sum(localData, 1)).^2; % [1 x Nsnr x Nmc x Nlocal]
                    case "MRC" % Matched Filter
                        globalTestStatistics = abs(sum(obj.fusionWeights{scanID}.*localData, 1)).^2;
                end
            end
        end

        %%% utility functions

        function randomSNRwithFixedAverage(obj, options)
            arguments
                obj
                options.averageSNR_dB double = 10*log10(cell2mat(cellfun(@(c) mean(10.^(.1*c), 1), obj.SNR_input_dB, 'UniformOutput', false).'))
                options.rangeSNR_dB (1, 1) double = 10
            end
            snrSamples = 10.^(.1*options.rangeSNR_dB*rand(max(obj.numberOfSensors), 1)); % between 0 dB to "rangeSNR_dB" dB
            for scanID = 1 : length(obj.numberOfSensors)
                M = obj.numberOfSensors(scanID);
                snr = snrSamples(1 : M)./mean(snrSamples(1 : M)); % mean of snr is now 0 dB 
                obj.SNR_input_dB{scanID} = 10*log10(snr) + options.averageSNR_dB(scanID, :);
            end
        end

        function randomLocalPFAwithFixedAverage(obj, options)
            arguments
                obj
                options.averageLocalPFA (1, 1) double = mean(cell2mat(cellfun(@(c) mean(c, 1), obj.localPFA, 'UniformOutput', false).'), 'all')
            end
            localPFAsamples = rand(max(obj.numberOfSensors), 1); % between 0 to 1
            for scanID = 1 : length(obj.numberOfSensors)
                M = obj.numberOfSensors(scanID);
                obj.localPFA{scanID} = localPFAsamples(1 : M)./mean(localPFAsamples(1 : M))*options.averageLocalPFA;
            end
        end

        function probability = totalCDF(obj, numberOfSamples, numberOfRequiredDetections, localThreshold, localProbability, fusionWeights, snr)
            if nargin < 7
                noiseEquation = true;
                snr = zeros(numberOfSamples, 1);
            else
                noiseEquation = false;
            end
            localThreshold = localThreshold(:);
            %%% ensure row vectors
            localProbability = localProbability(:).';
            fusionWeights = fusionWeights(:).';
            snr = snr(:).';

            uniqueLocalThreshold = unique(localThreshold);
            uniqueFusionWeights = unique(fusionWeights);
            uniqueSNR = unique(snr);
            uniqueRandomVariable = isscalar(uniqueSNR) && isscalar(uniqueLocalThreshold);

            switch obj.globalFusionRule
            %%% LRT: widely seperated receivers | spatially i.i.d circularly symmetric complex gaussian signal model
                case "SLC" % Square Law Combining (LRT under equal SNR)
                    if uniqueRandomVariable
                        if noiseEquation
                            probability = @(t) summationCCDF(t);
                        else
                            switch obj.signalPhaseModel
                                case "decorrelatedUniform"
                                    switch obj.signalAmplitudeModel
                                        case "decorrelatedExponential"
                                            % swerling-2 like spatially independent unit power signal
                                            % meaningful for widely seperated receivers
                                            probability = @(t) summationCCDF(t);
                                        case "correlatedExponential" % only for localPFA = 1
                                            % meaningful for phase asynchronous closely spaced receivers
                                            probability = @(t) summationCCDFcorrelatedAmplitude(t);
                                        case "deterministic"
                                            probability = @(t) summationCCDFdeterministic(t);
                                    end
                                case "correlatedUniform"
                                    switch obj.signalAmplitudeModel
                                        case "correlatedExponential" % only for localPFA = 1
                                            % swerling-1 like spatially fully correlated unit power signal
                                            % meaningful for phase synchronized closely spaced receivers
                                            probability = @(t) summationCCDFcorrelatedAmplitude(t);
                                        case "deterministic"
                                            error('not implemented');
                                    end
                            end
                        end
                    else
                        if noiseEquation
                            probability = @(t) generalizedSummationCCDF(t);
                        else
                            switch obj.signalPhaseModel
                                case "decorrelatedUniform"
                                    switch obj.signalAmplitudeModel
                                        case "decorrelatedExponential"
                                            % swerling-2 like spatially independent unit power signal
                                            % meaningful for widely seperated receivers
                                            probability = @(t) generalizedSummationCCDF(t);
                                        case "correlatedExponential" % only for localPFA = 1
                                            % meaningful for phase asynchronous closely spaced receivers
                                            error('not implemented');
                                        case "deterministic"
                                            error('not implemented');
                                    end
                                case "correlatedUniform"
                                    switch obj.signalAmplitudeModel
                                        case "correlatedExponential" % only for localPFA = 1
                                            % swerling-1 like spatially fully correlated unit power signal
                                            % meaningful for phase synchronized closely spaced receivers
                                            error('not implemented');
                                        case "deterministic"
                                            error('not implemented');
                                    end
                            end
                        end
                    end
                case "WSLC" % Weighted Sqaure Law Combining (LRT under different SNR)
                    probability = @(t) generalizedSummationCCDF(t);
            %%% Binary
                case "BC" % Binary Combining (LRT under equal SNR)
                    probability = @(t, p) summationBinaryCCDF(t, p);
                case "CVBC" % Chair-Varshney Binary Combining (LRT under different SNR)
                    if isscalar(uniqueFusionWeights) % integer threshold
                        probability = @(t, p) summationBinaryCCDF(t, p);
                    else
                        probability = @(t, p) generalizedSummationBinaryCCDF(t, p);
                    end
            %%% Others
                case "LLC" % Linear Law Combining
                    probability = @(t) nan;
                case "WLLC" % Weighted Linear Law Combining
                    probability = @(t) nan;
                case "NSLC" % Normalized Square Law Combining
                    if uniqueRandomVariable
                        probability = @(t) summationCCDF(t);
                    else
                        probability = @(t) generalizedSummationCCDF(t);
                    end
                case "SC" % Selective Combining (Partial Selective Combining K = 1)
                    probability = @(t) orderStatisticsCCDF(t);
                case "PSC" % Partial Selective Combining (Combine largest K samples)
                    probability = @(t) orderStatisticsCCDF(t);
                case {"EGC", "MRC"} % Equal Gain Combining
                    % for localPFA = 1 closed form exist
                    if uniqueRandomVariable && ~strcmp(obj.globalFusionRule, "MRC")
                        if noiseEquation
                            probability = @(t) cohSummationCCDFnonCohSignal(t); % only for localPFA = 1
                        else
                            switch obj.signalPhaseModel
                                case "decorrelatedUniform"
                                    switch obj.signalAmplitudeModel
                                        case "decorrelatedExponential"
                                            % swerling-2 like spatially independent unit power signal
                                            % meaningful for widely seperated receivers
                                            probability = @(t) cohSummationCCDFnonCohSignal(t); % only for localPFA = 1
                                        case "correlatedExponential"
                                            % meaningful for phase asynchronous closely spaced receivers
                                            probability = @(t) nan;
                                        case "deterministic"
                                            error('not implemented');
                                    end
                                case "correlatedUniform"
                                    switch obj.signalAmplitudeModel
                                        case "correlatedExponential"
                                            % swerling-1 like spatially fully correlated unit power signal
                                            % meaningful for phase synchronized closely spaced receivers
                                            probability = @(t) cohSummationCCDFcohSignal(t); % only for localPFA = 1
                                        case "deterministic"
                                            error('not implemented');
                                    end
                            end
                        end
                    else
                        if noiseEquation
                            probability = @(t) generalizedCohSummationCCDFnonCohSignal(t); % only for localPFA = 1
                        else
                            switch obj.signalPhaseModel
                                case "decorrelatedUniform"
                                    switch obj.signalAmplitudeModel
                                        case "decorrelatedExponential"
                                            % swerling-2 like spatially independent unit power signal
                                            % meaningful for widely seperated receivers
                                            probability = @(t) generalizedCohSummationCCDFnonCohSignal(t); % only for localPFA = 1
                                        case "correlatedExponential"
                                            % meaningful for phase asynchronous closely spaced receivers
                                            probability = @(t) nan;
                                        case "deterministic"
                                            error('not implemented');
                                    end
                                case "correlatedUniform"
                                    switch obj.signalAmplitudeModel
                                        case "correlatedExponential"
                                            % swerling-1 like spatially fully correlated unit power signal
                                            % meaningful for phase synchronized closely spaced receivers
                                            probability = @(t) generalizedCohSummationCCDFcohSignal(t); % only for localPFA = 1
                                        case "deterministic"
                                            error('not implemented');
                                    end
                            end
                        end
                    end
            end

            function probability = summationBinaryCCDF(t, p)
                p = p(:).'; % make sure of row vector
                pmf = poissonbinomialpdf(0 : min(floor(t), numberOfSamples), p);
                probability = 1 - sum(pmf);
            end

            function probability = generalizedSummationBinaryCCDF(t, p)
                if numberOfSamples <= 15
                    allSubsets = dec2bin(0 : 2^numberOfSamples - 1) - '0';
                    pmf = zeros(size(allSubsets, 1), 1);
                    weights = zeros(size(allSubsets, 1), 1);
                    for subSetID = 1 : size(allSubsets, 1)
                        indices = logical(allSubsets(subSetID, :));
                        pmf(subSetID) = prod(p(indices))*prod(1 - p(~indices));
                        weights(subSetID) = sum(fusionWeights(indices));
                    end
                    idx = weights <= t;
                    probability = 1 - sum(pmf(idx));
                else
                    error('#samples is huge for the calculation');
                end
            end

            function probability = cohSummationCCDFnonCohSignal(t)
                pmf = poissonbinomialpdf(0 : numberOfSamples, localProbability);
                probability = 0;
                for k = 0 : numberOfSamples
                    subSetProbability = pmf(k + 1);
                    if k == 0
                        conditionalCDF = double(t <= 0);
                    else
                        shift = k*uniqueLocalThreshold;
                        if any(shift)
                            %%% cannot implemented
                            conditionalCDF = nan;
                        else
                            switch obj.globalFusionRule
                                case "EGC"
                                    rate = 1/(k*(1 + uniqueSNR));
                                    conditionalCDF = gammaCCDF(t, rate, 0, 1);
                            end
                        end
                    end
                    probability = probability + subSetProbability*conditionalCDF;
                end
            end

            function probability = generalizedCohSummationCCDFnonCohSignal(t)
                if numberOfSamples <= 15
                    allSubsets = dec2bin(0 : 2^numberOfSamples - 1) - '0';
                    probability = 0;
                    for subSetID = 1 : size(allSubsets, 1)
                        indices = logical(allSubsets(subSetID, :));
                        subSetProbability = prod(localProbability(indices))*prod(1 - localProbability(~indices));
                        if subSetID == 1
                            conditionalCDF = double(t <= 0);
                        else
                            shifts = sum(localThreshold(indices));
                            if any(shifts)
                                %%% cannot implemented
                                conditionalCDF = nan;
                            else
                                K = nnz(indices);
                                switch obj.globalFusionRule
                                    case "EGC"
                                        rate = 1/(K + sum(snr(indices)));
                                        conditionalCDF = gammaCCDF(t, rate, 0, 1);
                                    case "MRC"
                                        rate = 1/(sum(fusionWeights(indices).^2) + sum(fusionWeights(indices).^2.*snr(indices)));
                                        conditionalCDF = gammaCCDF(t, rate, 0, 1);
                                end
                            end
                        end
                        probability = probability + subSetProbability*conditionalCDF;
                    end
                else
                    error('#samples is huge for the calculation');
                end
            end

            function probability = cohSummationCCDFcohSignal(t)
                pmf = poissonbinomialpdf(0 : numberOfSamples, localProbability);
                probability = 0;
                for k = 0 : numberOfSamples
                    subSetProbability = pmf(k + 1);
                    if k == 0
                        conditionalCDF = double(t <= 0);
                    else
                        shift = k*uniqueLocalThreshold;
                        if any(shift)
                            %%% cannot implemented
                            conditionalCDF = nan;
                        else
                            switch obj.globalFusionRule
                                case "EGC"
                                    rate = 1/(k*(1 + k*uniqueSNR));
                                    conditionalCDF = gammaCCDF(t, rate, 0, 1);
                            end
                        end
                    end
                    probability = probability + subSetProbability*conditionalCDF;
                end
            end

            function probability = generalizedCohSummationCCDFcohSignal(t)
                if numberOfSamples <= 15
                    allSubsets = dec2bin(0 : 2^numberOfSamples - 1) - '0';
                    probability = 0;
                    for subSetID = 1 : size(allSubsets, 1)
                        indices = logical(allSubsets(subSetID, :));
                        subSetProbability = prod(localProbability(indices))*prod(1 - localProbability(~indices));
                        if subSetID == 1
                            conditionalCDF = double(t <= 0);
                        else
                            shifts = sum(localThreshold(indices));
                            if any(shifts)
                                %%% cannot implemented
                                conditionalCDF = nan;
                            else
                                K = nnz(indices);
                                switch obj.globalFusionRule
                                    case "EGC"
                                        rate = 1/(K + sum(sqrt(snr(indices))).^2);
                                        conditionalCDF = gammaCCDF(t, rate, 0, 1);
                                    case "MRC"
                                        rate = 1/(sum(fusionWeights(indices).^2) + sum(fusionWeights(indices).*sqrt(snr(indices))).^2);
                                        conditionalCDF = gammaCCDF(t, rate, 0, 1);
                                end
                            end
                        end
                        probability = probability + subSetProbability*conditionalCDF;
                    end
                else
                    error('#samples is huge for the calculation');
                end
            end

            function probability = summationCCDF(t)
                pmf = poissonbinomialpdf(0 : numberOfSamples, localProbability);
                probability = 0;
                for k = 0 : numberOfSamples
                    subSetProbability = pmf(k + 1);
                    if k == 0
                        conditionalCDF = double(t <= 0);
                    else
                        switch obj.globalFusionRule
                            case "SLC" % Square Law Combining (LRT under equal SNR)
                                rate = 1/(1 + uniqueSNR);
                                shift = k*uniqueLocalThreshold;
                                conditionalCDF = gammaCCDF(t, rate, shift, k);
                            case "NSLC" % Normalized Square Law Combining
                                rate = k/(1 + uniqueSNR);
                                shift = uniqueLocalThreshold;
                                conditionalCDF = gammaCCDF(t, rate, shift, k);
                        end
                    end
                    probability = probability + subSetProbability*conditionalCDF;
                end
            end

            function probability = summationCCDFcorrelatedAmplitude(t)
                pmf = poissonbinomialpdf(0 : numberOfSamples, localProbability);
                probability = 0;
                for k = 0 : numberOfSamples
                    subSetProbability = pmf(k + 1);
                    if k == 0
                        conditionalCDF = double(t <= 0);
                    else
                        switch obj.globalFusionRule
                            case "SLC" % Square Law Combining (LRT under equal SNR)
                                rate = 1/(1 + k*uniqueSNR);
                                shift = k*uniqueLocalThreshold;
                                conditionalCDF = gammaCCDFcorrelatedAmplitude(t, rate, shift, k);
                            case "NSLC" % Normalized Square Law Combining
                                % rate = k/(1 + uniqueSNR);
                                % shift = uniqueLocalThreshold;
                                % conditionalCDF = gammaCCDF(t, rate, shift, k);
                        end
                    end
                    probability = probability + subSetProbability*conditionalCDF;
                end
            end

            function probability = summationCCDFdeterministic(t)
                pmf = poissonbinomialpdf(0 : numberOfSamples, localProbability);
                probability = 0;
                for k = 0 : numberOfSamples
                    subSetProbability = pmf(k + 1);
                    if k == 0
                        conditionalCDF = double(t <= 0);
                    else
                        switch obj.globalFusionRule
                            case "SLC" % Square Law Combining (LRT under equal SNR)
                                rate = k*uniqueSNR;
                                shift = k*uniqueLocalThreshold;
                                conditionalCDF = chiCCDF(t, rate, shift, k);
                        end
                    end
                    probability = probability + subSetProbability*conditionalCDF;
                end
            end

            function probability = generalizedSummationCCDF(t)
                if numberOfSamples <= 15
                    allSubsets = dec2bin(0 : 2^numberOfSamples - 1) - '0';
                    probability = 0;
                    for subSetID = 1 : size(allSubsets, 1)
                        indices = logical(allSubsets(subSetID, :));
                        subSetProbability = prod(localProbability(indices))*prod(1 - localProbability(~indices));
                        if subSetID == 1
                            conditionalCDF = double(t <= 0);
                        else
                            switch obj.globalFusionRule
                                case "SLC" % Square Law Combining (LRT under equal SNR)
                                    rates = 1./(1 + snr(indices));
                                    shifts = sum(localThreshold(indices));
                                    conditionalCDF = hypoexponentialCCDF(t, rates, shifts);
                                case "WSLC" % Weighted Sqaure Law Combining (LRT under different SNR)
                                    rates = 1./fusionWeights(indices)./(1 + snr(indices));
                                    shifts = fusionWeights(indices)*localThreshold(indices);
                                    conditionalCDF = hypoexponentialCCDF(t, rates, shifts);
                                case "NSLC" % Normalized Square Law Combining
                                    K = nnz(indices);
                                    rates = K./(1 + snr(indices));
                                    shifts = sum(localThreshold(indices))/K;
                                    conditionalCDF = hypoexponentialCCDF(t, rates, shifts);
                            end
                        end
                        probability = probability + subSetProbability*conditionalCDF;
                    end
                else
                    error('#samples is huge for the calculation');
                end
            end

            function probability = orderStatisticsCCDF(t)
                if numberOfSamples == 0
                    probability = double(t <= 0);
                else
                    probability = 0;
                    switch obj.globalFusionRule
                        case "SC" % Selective Combining (Partial Selective Combining K = 1)
                            rate = 1./(1 + snr);
                            shift = uniqueLocalThreshold;
                            argument = max(t, shift);
                            probability = 1 - prod(gammainc(rate*argument, 1, 'lower'), 2);
                        case "PSC"
                            % First term: special k = 0 case
                            rate = 1./(1 + uniqueSNR);
                            shift = uniqueLocalThreshold;
                            argument = max(t - shift, 0);
                            term0 = gammainc(rate*argument, min(numberOfSamples, numberOfRequiredDetections), 'lower');
                            if numberOfRequiredDetections > numberOfSamples
                                probability = 1 - term0;
                            else
                                sum_terms = 0;
                                if numberOfRequiredDetections > 0
                                    for k = 1 : (numberOfSamples - numberOfRequiredDetections)
                                        coeff = (-1)^k*nchoosek(numberOfSamples - numberOfRequiredDetections, k);
                                        term1 = gammainc(rate*argument*(k + 2)/numberOfRequiredDetections, 1, 'lower')/(k + 2);
                                        sum_terms = sum_terms + coeff*term1;
                                    end
                                end
                                if numberOfRequiredDetections > 1
                                    for k = 1 : (numberOfSamples - numberOfRequiredDetections)
                                        coeff = (-1)^k*nchoosek(numberOfSamples - numberOfRequiredDetections, k);
                                        term2 = gammainc(rate*argument*(k + 0)/numberOfRequiredDetections, 1, 'lower')/(k + 0).*exp(-rate*argument);
                                        sum_terms = sum_terms - coeff*term2;
                                    end
                                end
                                probability = 1 - nchoosek(numberOfSamples, numberOfRequiredDetections)*(term0 + numberOfRequiredDetections*sum_terms);
                            end
                    end
                end
            end

            %%%% probability mass functions

            function probability = poissonbinomialpdf(k, p)
                l = (0 : numberOfSamples).'/(numberOfSamples + 1); % [K x 1] vector
                phi = prod(1 + (exp(-1j*2*pi.*l) - 1).*p, 2); % [K x 1] vector
                probability = real(exp(-1j*2*pi*k.*l)'*phi).'/(numberOfSamples + 1);
                probability(isnan(probability) & k == 0) = 1;
                probability(isnan(probability) & k ~= 0) = 0;
            end

            %%%% complementary CDF functions

            function probability = gammaCCDF(t, rate, shift, degree)
                % complementary CDF of shifted gamma distribution
                argument = max(t - shift, 0);
                probability = gammainc(rate*argument, degree, 'upper');
            end

            function probability = chiCCDF(t, rate, shift, degree)
                % complementary CDF of shifted gamma distribution
                argument = max(2*t - shift, 0);
                probability = ncx2cdf(argument, 2*degree, 2*rate, 'upper');
            end

            function probability = hypoexponentialCCDF(t, rates, shifts)
                % complementary CDF of shifted hypoexponential distribution
                % with potentially repeated weights using the matrix exponential method.
                % for unique rates {hypoexponentialCCDF(t, rates) = gammainc(t.*unique(rates), numel(rates), 'upper')}
                if isempty(rates)
                    probability = double(t <= 0);
                else
                    argument = max(t - shifts, 0);
                    mustBeNonnegative(rates);
                    N = length(rates);
                    alpha = [1 zeros(1, N - 1)]; % Initial distribution (starting in state 1)
                    e = ones(N, 1); % Exit vector (absorbing state after last)
                    
                    %%% main calculation
                    Q = diag(-rates); % Construct generator matrix Q
                    Q(N + 1 : N + 1 : end) = rates(1 : end - 1);
                    probability = alpha*expm(Q*argument)*e;
                end
            end

            function probability = gammaCCDFcorrelatedAmplitude(t, rate, shift, degree)
                if any(shift)
                    %%% cannot implemented
                    probability = nan;
                else
                    alpha = 1/(1 - rate);
                    probability = gammainc(t, degree - 1, 'upper') + alpha^(degree - 1)*exp(-t*rate).*gammainc(t/alpha, degree - 1, 'lower'); % Richards
                end

                % argument = max(t - shift, 0);
                % n = 1 : (degree - 1);
                % summand = alpha.^(flip(n)).*gammainc(argument, n, 'upper');
                % probability = alpha^(degree - 1)*gammainc(rate*argument, 1, 'upper') - rate*sum(summand);
            end
        end

        %%% visualization

        function fig = visualize(obj, options)
            arguments
                obj
                options.dataType (1, :) string {mustBeMember(options.dataType, ["analytical", "empirical"])} = ["analytical", "empirical"]
                options.x_axis (1, :) string {mustBeMember(options.x_axis, ["globalPFA", "numberOfSensors", "SNR", "globalThreshold", "localPFA", "localThreshold", "localPD"])} = "numberOfSensors"
                options.y_axis (1, 1) string {mustBeMember(options.y_axis, ["globalPFA", "globalPD", "globalThreshold", ...
                    "numberOfActiveSensorsExpectationUnderNoise", "numberOfActiveSensorsStandardDeviationUnderNoise", ...
                    "numberOfActiveSensorsExpectationUnderSignal", "numberOfActiveSensorsStandardDeviationUnderSignal"])} = "globalPD"
            end
            averageSNR = cell2mat(cellfun(@(c) 10*log10(mean(10.^(.1*c), 1)), obj.SNR_input_dB, 'UniformOutput', false).');
            averageLocalPFA = permute(cell2mat(cellfun(@(c) mean(c, 1), obj.localPFA, 'UniformOutput', false).'), [1 4 2 3]);
            uniqueSNR = unique(averageSNR, 'rows', 'stable');
            uniqueLocalPFA = unique(averageLocalPFA, 'rows', 'stable'); % [1 x M] cell of [M x 1 x 1 x Nlocal] matrices
            if any(strcmp(options.x_axis, "globalThreshold"))
                if any(strcmp(obj.globalFusionRule, ["BC", "CVBC"]))
                    globalThresholddB = obj.globalThreshold;
                else
                    globalThresholddB = 10*log10(obj.globalThreshold);
                end
            end
            dataDimensions = ["globalPFA", "SNR", "numberOfSensors", "localPFA"];
            switch options.y_axis
                case "globalPFA"
                    yFullDataEmpirical = obj.(options.y_axis + "simulation"); % [NpfaGlobal x Nsnr x M x Nlocal]
                    yFullDataAnalytical = obj.globalPFAanalytical;
                    yLabel = 'P_{FA}^{global}';
                    yLabel = 'Global Probability of False Alarm';
                case "globalPD"
                    yFullDataEmpirical = obj.(options.y_axis + "simulation"); % [NpfaGlobal x Nsnr x M x Nlocal]
                    yFullDataAnalytical = obj.globalPD; % [NpfaGlobal x Nsnr x M x Nlocal]
                    yLabel = 'P_{D}^{global}';
                    yLabel = 'Global Probability of Detection';
                case "globalThreshold"
                    if any(strcmp(obj.globalFusionRule, ["BC", "CVBC"]))
                        yFullDataEmpirical = obj.(options.y_axis + "Simulation"); % [NpfaGlobal x Nsnr x M x Nlocal]
                        yFullDataAnalytical = obj.globalThreshold; % [NpfaGlobal x Nsnr x M x Nlocal]
                        yLabel = 'T_{global}';
                        yLabel = 'Global Threshold';
                    else
                        yFullDataEmpirical = 10*log10(obj.(options.y_axis + "Simulation")); % [NpfaGlobal x Nsnr x M x Nlocal]
                        yFullDataAnalytical = 10*log10(obj.globalThreshold); % [NpfaGlobal x Nsnr x M x Nlocal]
                        yLabel = 'T_{global} (dB)';
                        yLabel = 'Global Threshold (dB)';
                    end
                case "numberOfActiveSensorsExpectationUnderNoise"
                    yFullDataEmpirical = obj.(options.y_axis + "Simulation"); % [1 x 1 x M x Nlocal]
                    yFullDataAnalytical = obj.numberOfActiveSensorsExpectationUnderNoise; % [1 x 1 x M x Nlocal]
                    yLabel = 'E[M|H_0]';
                case "numberOfActiveSensorsStandardDeviationUnderNoise"
                    yFullDataEmpirical = obj.(options.y_axis + "Simulation"); % [1 x 1 x M x Nlocal]
                    yFullDataAnalytical = obj.numberOfActiveSensorsStandardDeviationUnderNoise; % [1 x 1 x M x Nlocal]
                    yLabel = '\sigma_M|H_1';
                case "numberOfActiveSensorsExpectationUnderSignal"
                    yFullDataEmpirical = obj.(options.y_axis + "Simulation"); % [1 x Nsnr x M x Nlocal]
                    yFullDataAnalytical = obj.numberOfActiveSensorsExpectationUnderSignal; % [1 x Nsnr x M x Nlocal]
                    yLabel = 'E[M|H_1]';
                case "numberOfActiveSensorsStandardDeviationUnderSignal"
                    yFullDataEmpirical = obj.(options.y_axis + "Simulation"); % [1 x Nsnr x M x Nlocal]
                    yFullDataAnalytical = obj.numberOfActiveSensorsStandardDeviationUnderSignal; % [1 x Nsnr x M x Nlocal]
                    yLabel = '\sigma_M|H_1';
            end
            if isempty(yFullDataAnalytical) && any(strcmp(options.dataType, "analytical"))
                return;
            end
            if isempty(yFullDataEmpirical) && any(strcmp(options.dataType, "empirical"))
                return;
            end
            yTitle = yLabel;
            dataSize = size(yFullDataAnalytical, [1 2 3 4]);

            for plotID = 1 : length(options.x_axis)
                plotFunc = @plot;
                legendString = '';
                legendTitleString = '';
                switch options.x_axis(plotID)
                    case "globalPFA"
                        plotFunc = @semilogx;
                        xData = obj.globalPFA; % [NpfaGlobal x 1]
                        if isscalar(xData)
                            fprintf('choose non-scalar x-axis \n');
                            continue
                        end
                        xLabel = 'P_{FA}^{global}';
                        xLabel = 'Global Probability of False Alarm';
                        xTitle = xLabel;
                        indexingPriority = ["localPFA", "numberOfSensors", "SNR", "globalPFA"];
                        legendPriority = ["SNR", "numberOfSensors", "localPFA"];
                    case "numberOfSensors"
                        xData = obj.numberOfSensors; % [1 x M]
                        if isscalar(xData)
                            fprintf('choose non-scalar x-axis \n');
                            continue
                        end
                        xLabel = '#sensors';
                        xLabel = 'Total Number of Receiving Nodes';
                        xTitle = xLabel;
                        indexingPriority = ["globalPFA", "localPFA", "SNR", "numberOfSensors"];
                        legendPriority = ["SNR", "localPFA", "globalPFA"];
                    case "SNR"
                        xData = uniqueSNR.'; % [1 x Nsnr]
                        if isscalar(xData)
                            fprintf('choose non-scalar x-axis \n');
                            continue
                        elseif isvector(xData)
                            if max(abs(diff(xData))) < 20*eps
                                fprintf('choose non-scalar x-axis \n');
                                continue
                            end
                        else
                            xData = averageSNR.';
                        end
                        xLabel = 'SNR_{in} (dB)';
                        xLabel = 'SNR after Matched Filtering (dB)';
                        xTitle = xLabel;
                        indexingPriority = ["globalPFA", "localPFA", "numberOfSensors", "SNR"];
                        legendPriority = ["numberOfSensors", "localPFA", "globalPFA"];
                    case "globalThreshold"
                        xData = globalThresholddB; % [NpfaGlobal x 1 x M x Nlocal]
                        if isscalar(xData)
                            fprintf('choose non-scalar x-axis \n');
                            continue
                        elseif isvector(xData)
                            if max(abs(diff(xData))) < 20*eps
                                fprintf('choose non-scalar x-axis \n');
                                continue
                            end
                        end
                        if any(strcmp(obj.globalFusionRule, ["BC", "CVBC"]))
                            xLabel = 'T_{global}';
                            xLabel = 'Global Threshold';
                        else
                            xLabel = 'T_{global} (dB)';
                            xLabel = 'Global Threshold (dB)';
                        end
                        xTitle = xLabel;
                        indexingPriority = ["globalPFA", "localPFA", "SNR", "numberOfSensors"];
                        legendPriority = ["SNR", "localPFA", "numberOfSensors", "globalPFA"];
                    case {"localPFA", "localThreshold", "localPD"}
                        switch options.x_axis(plotID)
                            case "localPFA"
                                plotFunc = @semilogx;
                                xData = uniqueLocalPFA; % [1 x Nlocal] or [M x Nlocal]
                                xLabel = 'P_{FA}^{local}';
                                xLabel = 'Local Probability of False Alarm';
                                if isscalar(xData)
                                    fprintf('choose non-scalar x-axis \n');
                                    continue
                                elseif isvector(xData)
                                    if max(abs(diff(xData))) < 20*eps
                                        fprintf('choose non-scalar x-axis \n');
                                        continue
                                    end
                                else
                                    xData = averageLocalPFA; % [M x Nlocal]
                                end
                                legendPriority = ["numberOfSensors", "SNR", "globalPFA"];
                            case "localThreshold"
                                % [1 x M] cell of [M x 1 x 1 x Nlocal] matrices
                                averageThreshold = permute(cell2mat(cellfun(@(c) mean(c, 1), obj.localThreshold, 'UniformOutput', false).'), [1 4 2 3]);
                                uniqueLocalThreshold = unique(averageThreshold, 'rows', 'stable');
                                xData = uniqueLocalThreshold; % [1 x Nlocal] or [M x Nlocal]
                                xLabel = 'T_{local} (dB)';
                                xLabel = 'Local Threshold (dB)';
                                if isscalar(xData)
                                    fprintf('choose non-scalar x-axis \n');
                                    continue
                                elseif isvector(xData)
                                    if max(abs(diff(xData))) < 20*eps
                                        fprintf('choose non-scalar x-axis \n');
                                        continue
                                    end
                                else
                                    xData = averageThreshold; % [M x Nlocal]
                                end
                                legendPriority = ["numberOfSensors", "SNR", "globalPFA", "localPFA"];
                            case "localPD"
                                % [1 x M] cell of [M x Nsnr x 1 x Nlocal] matrices
                                averageLocalPD = reshape(permute(cell2mat(cellfun(@(c) mean(c, 1), obj.localPD, 'UniformOutput', false).'), [1 2 4 3]), length(obj.numberOfSensors), []);
                                uniqueLocalPD = unique(averageLocalPD, 'rows', 'stable'); % [1 x Nsnr*Nlocal] or % [M x Nsnr*Nlocal]
                                xData = uniqueLocalPD;
                                if isscalar(xData)
                                    fprintf('choose non-scalar x-axis \n');
                                    continue
                                elseif isvector(xData)
                                    if max(abs(diff(xData))) < 20*eps
                                        fprintf('choose non-scalar x-axis \n');
                                        continue
                                    end
                                else
                                    xData = reshape(averageLocalPD, length(obj.numberOfSensors), size(uniqueSNR, 2), []); % [M x Nsnr x Nlocal]
                                end
                                xLabel = 'P_{D}^{local}';
                                xLabel = 'Local Probability of Detection';
                                if isscalar(xData)
                                    fprintf('choose non-scalar x-axis \n');
                                    continue
                                end
                                legendPriority = ["numberOfSensors", "SNR", "globalPFA", "localPFA"];
                        end
                        xTitle = xLabel;
                        indexingPriority = ["globalPFA", "numberOfSensors", "SNR", "localPFA"];
                end
                if ~any(strcmp(options.x_axis(plotID), ["globalPFA", "SNR", "numberOfSensors", "localPFA"]))
                    if numel(xData) == length(uniqueLocalPFA(1, :))
                        legendPriority = setdiff(legendPriority, "globalPFA");
                    end
                    if numel(xData) == length(uniqueSNR)
                        legendPriority = setdiff(legendPriority, "SNR");
                    end
                    if numel(xData) == length(uniqueLocalPFA(1, :))
                        legendPriority = setdiff(legendPriority, "localPFA");
                    end
                end
                if numel(xData) == length(obj.numberOfSensors)
                    legendPriority = setdiff(legendPriority, "numberOfSensors");
                end
                if strcmp(options.y_axis, "globalPFA") || (strcmp(options.y_axis, "globalThreshold") && ~any(strcmp(obj.globalFusionRule, ["BC", "CVBC"])))
                    if isequal(plotFunc, @semilogx)
                        plotFunc = @loglog;
                    else
                        plotFunc = @semilogy;
                    end
                end
                indexing = arrayfun(@(i) 1 : dataSize(i), 1 : length(dataSize), 'UniformOutput', false);
                indexingDimension = ismember(dataDimensions, indexingPriority);
                if nnz(dataSize(indexingDimension) ~= 1) > 2
                    indices = ones(1, 4);
                    indices(3) = ceil(dataSize(3)/2); % # sensors
                    possibleIndicedDimensions = dataDimensions(indexingDimension);
                    dataDimension = nnz(dataSize(indexingDimension) ~= 1);
                    indexingPriority = indexingPriority();
                    indexingDimension = false(1, length(dataDimensions));
                    numberOfIndicedDimensions = 0;
                    for priorityID = 1 : length(dataDimensions)
                        if dataSize(strcmp(dataDimensions, indexingPriority(priorityID))) ~= 1
                            indexingDimension(possibleIndicedDimensions == indexingPriority(priorityID)) = true;
                            numberOfIndicedDimensions = numberOfIndicedDimensions + 1;
                        end
                        if numberOfIndicedDimensions == dataDimension - 2
                            break;
                        end
                    end
                    indexing(indexingDimension) = num2cell(indices(indexingDimension));
                    yDataEmpirical = squeeze(yFullDataEmpirical(indexing{:}));
                    yDataAnalytical = squeeze(yFullDataAnalytical(indexing{:}));
                    legendPriority = legendPriority(find(ismember(legendPriority, dataDimensions(cellfun(@(c) ~isscalar(c), indexing))), 1, 'first'));
                    switch options.x_axis(plotID)
                        case "globalThreshold"
                            indexing{2} = 1; % SNR
                            xData = squeeze(xData(indexing{:}));
                        case {"localPFA", "localThreshold"}
                            if ~isvector(xData) && indexing{1} == 1 % #Sensors
                                xData = squeeze(xData(indexing{[1 4]})); % [M x Nlocal] --> [1 x Nlocal]
                            end
                        case "localPD"
                            if ~isvector(xData)
                                xData = squeeze(xData(indexing{[1 2 4]})); % [M x Nsnr x Nlocal] --> [1 x Nlocal]
                            end
                    end
                    if indexingDimension(1) || dataSize(1) == 1
                        globalPFAString = [', P_{FA}^{global} = ', scinot(obj.globalPFA(1))];
                    else
                        globalPFAString = '';
                    end
                    if indexingDimension(2) || dataSize(2) == 1
                        SNRString = [', SNR = ', sprintf('%g', uniqueSNR(1)), ' dB'];
                    else
                        SNRString = '';
                    end
                    if indexingDimension(3) || dataSize(3) == 1
                        numberOfSensorsString = [', #sensors = ', sprintf('%d', obj.numberOfSensors(ceil(end/2)))];
                    else
                        numberOfSensorsString = '';
                    end
                    if indexingDimension(4) || dataSize(4) == 1
                        localPFAString = [', P_{FA}^{local} = ', scinot(uniqueLocalPFA(1))];
                    else
                        localPFAString = '';
                    end
                else
                    yDataDimensions = dataSize ~= 1;
                    yDataEmpirical = squeeze(yFullDataEmpirical);
                    yDataAnalytical = squeeze(yFullDataAnalytical);
                    switch options.x_axis(plotID)
                        case "globalThreshold"
                            xData = squeeze(xData);
                        case {"localPFA", "localThreshold", "localPD"}
                            if ~isvector(xData)
                                % [M x 1 x Nlocal] --> [Nlocal x M]
                                % [1 x Nsnr x Nlocal] --> [Nlocal x Nsnr]
                                % [M x Nsnr x 1] --> [Nsnr x M]
                                xData = squeeze(xData);
                                xData = xData.';
                                yDataEmpirical = yDataEmpirical.';
                                yDataAnalytical = yDataAnalytical.';
                            end
                    end
                    legendPriority = legendPriority(find(ismember(legendPriority, dataDimensions(yDataDimensions)), 1, 'first'));
                    if dataSize(1) == 1
                        globalPFAString = [', P_{FA}^{global} = ', scinot(obj.globalPFA(1))];
                    else
                        globalPFAString = '';
                    end
                    if dataSize(2) == 1
                        SNRString = [', SNR = ', sprintf('%g', uniqueSNR(1)), ' dB'];
                    else
                        SNRString = '';
                    end
                    if dataSize(3) == 1
                        numberOfSensorsString = [', #sensors = ', sprintf('%d', obj.numberOfSensors(ceil(end/2)))];
                    else
                        numberOfSensorsString = '';
                    end
                    if dataSize(4) == 1
                        localPFAString = [', P_{FA}^{local} = ', scinot(uniqueLocalPFA(1))];
                    else
                        localPFAString = '';
                    end
                end
                indexingTitle = [globalPFAString, SNRString, numberOfSensorsString, localPFAString];
                if ~isvector(yDataEmpirical) && ~isvector(yDataAnalytical)
                    switch legendPriority
                        case "numberOfSensors"
                            legendString = num2str(obj.numberOfSensors');
                            legendTitleString = '#sensors';
                            legendTitleString = 'M';
                        case "SNR"
                            legendString = num2str(uniqueSNR(1, :).');
                            legendTitleString = 'SNR_{in} (dB)';
                            legendTitleString = '\eta (dB)';
                        case "globalPFA"
                            legendString = num2str(obj.globalPFA);
                            legendTitleString = 'P_{FA}^{global}';
                        case "localPFA"
                            legendString = scinot(uniqueLocalPFA(ceil(end/2), :).', 0);
                            legendTitleString = 'P_{FA}^{local}';
                    end
                end
                titleString = [yTitle, ' vs ', xTitle, indexingTitle];
                fig = figure;
                if any(strcmp(options.dataType, "empirical"))
                    subTitleString = ['#trials = ', scinot(obj.numberOfTrials)];
                    plotFunc(xData, yDataEmpirical, '-', 'LineWidth', 2);
                    hold on;
                else
                    subTitleString = '';
                end
                if any(strcmp(options.dataType, "analytical"))
                    if any(strcmp(options.dataType, "empirical"))
                        plotFunc(xData, yDataAnalytical, '--k');
                    else
                        plotFunc(xData, yDataAnalytical, '-', 'LineWidth', 2);
                        hold on;
                    end
                end
                if any(strcmp(options.y_axis, ["globalPFA", "globalPD"]))
                    ylim([0, 1]);
                else
                    ylim tight;
                end
                xlim tight;
                grid off; grid minor; grid on;
                xlabel(xLabel); ylabel(yLabel);
                % title(titleString, subTitleString);
                if ~isempty(legendString)
                    leg = legend(legendString, 'Location', 'best');
                    title(leg, legendTitleString);
                end
                hold off; drawnow;
            end

            function str = scinot(value, significantFigure)
                if nargin < 2
                    significantFigure = 2;
                end
                numberOfValues = numel(value);
                significand = 10.^mod(log10(value), 1);
                exponent = floor(log10(value));
                str = strings(1, numberOfValues);
                for valueID = 1 : numberOfValues
                    if sign(exponent(valueID)) == -1
                        str(valueID) = sprintf(['%0.' num2str(significantFigure) 'f•10^{-%i}'], significand(valueID), abs(exponent(valueID)));
                    else
                        str(valueID) = sprintf(['%0.' num2str(significantFigure) 'f•10^%i'], significand(valueID), exponent(valueID));
                    end
                end
            end
        end
    end
end