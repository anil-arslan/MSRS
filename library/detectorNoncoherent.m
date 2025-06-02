classdef detectorNoncoherent < handle
    %detector Summary of this class goes here
    %   Detailed explanation goes here

    %%% Model parameters

    properties (SetAccess = private, GetAccess = public)
        seed (1, 1) double {mustBeNonnegative, mustBeInteger, mustBeInRange(seed, 0, 4294967295)} = 0
        numberOfTrials (1, 1) double {mustBeNonnegative, mustBeInteger} = 1e3
        SNR_input_dB (1, :) cell = {4} % [1 x M] cell of [M x Nsnr] matrices
        numberOfSensors (1, :) double {mustBeNonnegative, mustBeInteger} = 10 % [1 x M]
    end

    %%% Fusion parameters

    properties (SetAccess = private, GetAccess = public)
        numberOfBinaryDetections double {mustBeNonnegative, mustBeInteger} = 1 % [1 x M]
        globalFusionRule (1, 1) string {mustBeMember(globalFusionRule, ["EGC", "MRC", "LDC", "GLDC", "BC", "CVBC", "NEGC", "SC", "PSC"])} = "EGC"
        numberOfBinaryDetectionsRule (1, 1) string {mustBeMember(numberOfBinaryDetectionsRule, ["and", "or", "majority", "fixedGlobalPFA", "userDefined"])} = "fixedGlobalPFA"
        %%% LRT
        % EGC: Equal Gain Combining (LRT under equal SNR)
        % MRC: Maximal Ratio Combining (LRT under different SNR)

        %%% GLRT
        % LDC: Log-Divergence Combining (GLRT w/ joint MLE)
        % GLDC: Generalized Log-Divergence Combining (GLRT w/ independent MLE)

        %%% Binary
        % BC: Binary Combining (LRT under equal SNR)
        % CVBC: Chair-Varshney Binary Combining (LRT under different SNR)

        %%%
        % NEGC: Normalized Equal Gain Combining
        % SC: Selective Combining (Partial Selective Combining K = 1)
        % PSC: Partial Selective Combining (Combine largest K samples)
        fusionWeights (1, :) cell % [1 x M] cell of [M x Nsnr x 1 x Nlocal] matrices
        fusionBias (1, :) cell % [1 x M] cell of [M x Nsnr x 1 x Nlocal] matrices
    end

    properties (SetAccess = private, GetAccess = public)
        globalPFA (:, 1) double {mustBeNonnegative, mustBeInRange(globalPFA, 0, 1)} = 1e-6 % [NpfaGlobal x 1]
        globalPFAsimulation double {mustBeNonnegative, mustBeInRange(globalPFAsimulation, 0, 1)} = [] % [NpfaGlobal x 1 x M x Nlocal]
    end

    properties (Dependent)
        globalThreshold double {mustBeNonnegative} % [NpfaGlobal x Nsnr x M x Nlocal]
        globalPFAanalytical double {mustBeNonnegative} % [NpfaGlobal x Nsnr x M x Nlocal]
        globalPD double {mustBeNonnegative, mustBeInRange(globalPD, 0, 1)} % [NpfaGlobal x Nsnr x M x Nlocal]
    end

    properties (SetAccess = private, GetAccess = public)
        globalPDsimulation double {mustBeNonnegative, mustBeInRange(globalPDsimulation, 0, 1)} = [] % [NpfaGlobal x Nsnr x M x Nlocal]
    end

    %%% Local detection parameters

    properties (SetAccess = private, GetAccess = public)
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

    %%% Data generation

    properties (Dependent, Hidden)
        noise
        signal
    end

    methods
        function obj = detectorNoncoherent(options)
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
            obj.numberOfBinaryDetections = ceil((obj.numberOfSensors + 1)/2);
            numberOfScansPFAglobal = length(obj.globalPFA);
            numberOfScans = length(obj.numberOfSensors);
            obj.SNR_input_dB = cell(1, numberOfScans); % [1 x M] cell
            if ~iscell(options.SNR_input_dB)
                options.SNR_input_dB = {options.SNR_input_dB};
            else
                mustBeVector(options.SNR_input_dB);
            end
            if isscalar(options.SNR_input_dB)
                options.SNR_input_dB = cell2mat(options.SNR_input_dB);
                for scanID = 1 : numberOfScans
                    if size(options.SNR_input_dB, 1) == 1 % same SNR across all sensors
                        obj.SNR_input_dB{scanID} = repmat(options.SNR_input_dB, obj.numberOfSensors(scanID), 1);
                    else
                        obj.SNR_input_dB{scanID} = options.SNR_input_dB(1 : obj.numberOfSensors(scanID), :);
                    end
                end
            else
                mustBeScalarOrEmpty(unique((cellfun(@(c) size(c, 2), options.SNR_input_dB)))); % to ensure same Nsnr
                for scanID = 1 : numberOfScans
                    SNRmatrix = options.SNR_input_dB{scanID};
                    firstSizeSNRmatrix = size(SNRmatrix, 1);
                    assert(firstSizeSNRmatrix == obj.numberOfSensors(scanID) || firstSizeSNRmatrix == 1);
                    if firstSizeSNRmatrix == 1 % same SNR across all sensors
                        SNRmatrix = repmat(SNRmatrix, obj.numberOfSensors(scanID), 1);
                    end
                    obj.SNR_input_dB{scanID} = SNRmatrix;
                end
            end
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
            obj.globalPFAsimulation = zeros(numberOfScansPFAglobal, 1, numberOfScans, numberOfScansLocalPFA);
            obj.globalPDsimulation = zeros(numberOfScansPFAglobal, numberOfScansSNR, numberOfScans, numberOfScansLocalPFA);
            obj.localPFAsimulation = cell(1, numberOfScans);
            obj.localPDsimulation = cell(1, numberOfScans);
            for scanID = 1 : numberOfScans
                obj.localPFAsimulation{scanID} = zeros(numberOfScans, 1, 1, numberOfScansLocalPFA); % [1 x M] cell of [M x 1 x 1 x Nlocal] matrices
                obj.localPDsimulation{scanID} = zeros(numberOfScans, numberOfScansSNR, 1, numberOfScansLocalPFA); % [1 x M] cell of [M x Nsnr x 1 x Nlocal] matrices
            end
        end

        %%% Fusion rule

        function K = get.numberOfBinaryDetections(obj)
            switch obj.numberOfBinaryDetectionsRule
                case "and"
                    K = obj.numberOfSensors;
                case "or"
                    K = ones(1, length(obj.numberOfSensors));
                case "majority"
                    K = ceil((obj.numberOfSensors + 1)/2);
                case "fixedGlobalPFA"
                    K = obj.globalThreshold;
                case "userDefined"
                    K = obj.numberOfBinaryDetections;
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
                case "MRC"
                    for scanID = 1 : numberOfScans
                        snr = 10.^(.1*obj.SNR_input_dB{scanID});
                        w{scanID} = repmat(snr./(1 + snr), [1 1 1 numberOfScansLocalPFA]);
                    end
                case "CVBC"
                    for scanID = 1 : numberOfScans
                        pfaLocal = obj.localPFA{scanID};
                        pdLocal = obj.localPD{scanID}; % requires knowledge of SNR
                        w{scanID} = log((pdLocal.*(1 - pfaLocal))./(pfaLocal.*(1 - pdLocal)));
                        % lets make the sum of weights equal to #sensors so that CVBC and BC has the same maximum test statistics value
                        w{scanID} = w{scanID}./sum(w{scanID}, 1)*obj.numberOfSensors(scanID);
                        w{scanID}(isnan(w{scanID})) = 1;
                    end
                otherwise
                    for scanID = 1 : numberOfScans
                        w{scanID} = ones(obj.numberOfSensors(scanID), numberOfScansSNR, 1, numberOfScansLocalPFA);
                    end
            end
        end

        function bias = get.fusionBias(obj)
            % [1 x M] cell of [M x Nsnr x 1 x Nlocal] matrices
            numberOfScansSNR = unique(cellfun(@(c) size(c, 2), obj.SNR_input_dB));
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            bias = cell(1, numberOfScans);
            switch obj.globalFusionRule
                case "MRC"
                    for scanID = 1 : numberOfScans
                        snr = 10.^(.1*obj.SNR_input_dB{scanID});
                        bias{scanID} = repmat(-log(1 + snr), [1 1 1 numberOfScansLocalPFA]);
                    end
                otherwise
                    for scanID = 1 : numberOfScans
                        bias{scanID} = zeros(obj.numberOfSensors(scanID), numberOfScansSNR, 1, numberOfScansLocalPFA);
                    end
            end
        end

        function Tglobal = get.globalThreshold(obj)
            % [NpfaGlobal x Nsnr x M x Nlocal]
            Tlocal = obj.localThreshold;
            weights = obj.fusionWeights;
            biases = obj.fusionBias;
            numberOfScansSNR = unique(cellfun(@(c) size(c, 2), obj.SNR_input_dB));
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            Tglobal = zeros(length(obj.globalPFA), 1, length(obj.numberOfSensors), numberOfScansLocalPFA);
            for scanID = 1 : numberOfScans
                M = obj.numberOfSensors(scanID);
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
                            b = biases{scanID}(:, snrID, 1, localPFAID);
                            F = obj.totalCDF(M, lambda, pfaLocal, w, b);
                            switch obj.globalFusionRule
                            %%% LRT
                                % EGC: Equal Gain Combining (LRT under equal SNR)
                                % MRC: Maximal Ratio Combining (LRT under different SNR)
                        
                            %%% GLRT
                                % LDC: Log-Divergence Combining (GLRT w/ joint MLE)
                                % GLDC: Generalized Log-Divergence Combining (GLRT w/ independent MLE)
                        
                            %%% Binary
                                % BC: Binary Combining (LRT under equal SNR)
                                % CVBC: Chair-Varshney Binary Combining (LRT under different SNR)
                        
                            %%%
                                % NEGC: Normalized Equal Gain Combining
                                % SC: Selective Combining (Partial Selective Combining K = 1)
                                % PSC: Partial Selective Combining (Combine largest K samples)
                                case {"EGC", "MRC", "LDC", "GLDC", "NEGC", "SC", "PSC"}
                                    thresholdNSC = gammaincinv(pfaGlobal, M, 'upper'); % equal gain non-Selective combiner threshold, conventional
                                    thresholdSearchSpace = [0, sum(w.*lambda) + 2*thresholdNSC] + sum(b);
                                case {"BC", "CVBC"} %%% Binary
                                    thresholdSearchSpace = [0, sum(w)] + sum(b);
                            end
                            switch obj.globalFusionRule
                                case {"EGC", "MRC", "LDC", "NEGC"} % float threshold
                                    Tglobal(PFAID, snrID, scanID, localPFAID) = fzero(@(t) F(t) - pfaGlobal, thresholdSearchSpace);
                                case "BC" % integer threshold
                                    switch obj.numberOfBinaryDetectionsRule
                                        case "fixedGlobalPFA"
                                            pfaGlobalK = zeros(1, M + 1);
                                            for K = 0 : M
                                                pfaGlobalK(K + 1) = F(K);
                                            end
                                            % since global PFA monotonically decreasing with K
                                            Tglobal(PFAID, snrID, scanID, localPFAID) = find(pfaGlobalK < pfaGlobal, 1, 'first') - 0.5;
                                            % - 0.5 prevents boundary problems
                                        otherwise
                                            Tglobal(PFAID, snrID, scanID, localPFAID) = obj.numberOfBinaryDetections(scanID) - 0.5;
                                    end
                                case "CVBC"
                                    for subSetID = 1 : numberOfSubsets
                                        indices = logical(allSubsets(subSetID, :));
                                        pmf(subSetID) = prod(pfaLocal(indices))*prod(1 - pfaLocal(~indices));
                                        weightSet(subSetID) = sum(w(indices));
                                    end
                                    [weightSet, sortIndices] = sort(weightSet);
                                    cdf = 1 - cumsum(pmf(sortIndices));
                                    thresholdIndex = find(cdf < pfaGlobal, 1, 'first');
                                    if thresholdIndex ~= numberOfSubsets
                                        threshold = mean(weightSet([thresholdIndex, thresholdIndex + 1]));
                                    else
                                        threshold = weightSet(thresholdIndex) + 0.5;
                                    end
                                    Tglobal(PFAID, snrID, scanID, localPFAID) = threshold;
                                case "SC" % inverse exists
                                    Tglobal(PFAID, snrID, scanID, localPFAID) = max(unique(lambda), -log(1 - (1 - pfaGlobal).^(1/M)));
                                case "GLDC"
                                    error('not implemented');
                                case "PSC"
                                    error('not implemented');
                            end
                        end
                    end
                end
            end
        end
        
        function pfa = get.globalPFAanalytical(obj)
            % [NpfaGlobal x Nsnr x M x Nlocal]
            Tglobal = obj.globalThreshold;
            Tlocal = obj.localThreshold;
            pFAlocal = obj.localPFA;
            weights = obj.fusionWeights;
            biases = obj.fusionBias;
            numberOfScansSNR = unique(cellfun(@(c) size(c, 2), obj.SNR_input_dB));
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            pfa = zeros(length(obj.globalPFA), numberOfScansSNR, length(obj.numberOfSensors), numberOfScansLocalPFA);
            for scanID = 1 : numberOfScans
                M = obj.numberOfSensors(scanID);
                for PFAID = 1 : length(obj.globalPFA)
                    for localPFAID = 1 : numberOfScansLocalPFA
                        gamma = Tglobal(PFAID, 1, scanID, localPFAID);
                        lambda = Tlocal{scanID}(:, :, :, localPFAID);
                        pfaLocal = pFAlocal{scanID}(:, :, :, localPFAID);
                        for snrID = 1 : numberOfScansSNR
                            w = weights{scanID}(:, snrID, 1, localPFAID);
                            b = biases{scanID}(:, snrID, 1, localPFAID);
                            F = obj.totalCDF(M, lambda, pfaLocal, w, b);
                            pfa(PFAID, snrID, scanID, localPFAID) = F(gamma);
                        end
                        % evaluating @(gamma - 1) is important for binary fusion with integer threshold
                    end
                end
            end
        end

        function ROC = get.globalPD(obj)
            % [NpfaGlobal x Nsnr x M x Nlocal]
            Tglobal = obj.globalThreshold;
            Tlocal = obj.localThreshold;
            pDlocal = obj.localPD;
            weights = obj.fusionWeights;
            biases = obj.fusionBias;
            numberOfScansSNR = unique(cellfun(@(c) size(c, 2), obj.SNR_input_dB));
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            ROC = zeros(length(obj.globalPFA), numberOfScansSNR, numberOfScans, numberOfScansLocalPFA);
            for scanID = 1 : numberOfScans
                M = obj.numberOfSensors(scanID);
                snr = 10.^(.1*obj.SNR_input_dB{scanID});
                for PFAID = 1 : length(obj.globalPFA)
                    for localPFAID = 1 : numberOfScansLocalPFA
                        gamma = Tglobal(PFAID, 1, scanID, localPFAID);
                        lambda = Tlocal{scanID}(:, :, :, localPFAID);
                        for snrID = 1 : numberOfScansSNR
                            pdLocal = pDlocal{scanID}(:, snrID, 1, localPFAID);
                            w = weights{scanID}(:, snrID, 1, localPFAID);
                            b = biases{scanID}(:, snrID, 1, localPFAID);
                            F = obj.totalCDF(M, lambda, pdLocal, w, b, snr(:, snrID));
                            ROC(PFAID, snrID, scanID, localPFAID) = F(gamma);
                            % evaluating @(gamma - 1) is important for binary fusion with integer threshold
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
            theta = cellfun(@(c) 1./(1 + 10.^(.1*c)), obj.SNR_input_dB, 'UniformOutput', false); % [1 x M] cell of [M x Nsnr] matrices
            ROC = cellfun(@(c, t) c.^t, obj.localPFA, theta, 'UniformOutput', false);
        end

        %%% Data generation

        function n = get.noise(obj)
            % [M x 1 x Nmc]
            n = (randn(max(obj.numberOfSensors), 1, obj.numberOfTrials) + 1j*randn(max(obj.numberOfSensors), 1, obj.numberOfTrials))/sqrt(2);
        end

        function s = get.signal(obj)
            % [1 x M] cell of [M x Nsnr] matrices
            % swerling-2 like spatially independent unit power signal
            complexGaussian = (randn(max(obj.numberOfSensors), 1, obj.numberOfTrials) + 1j*randn(max(obj.numberOfSensors), 1, obj.numberOfTrials))/sqrt(2);
            numberOfScans = length(obj.numberOfSensors);
            s = cell(1, numberOfScans);
            for scanID = 1 : numberOfScans
                s{scanID} = 10.^(.05*obj.SNR_input_dB{scanID}).*complexGaussian(1 : obj.numberOfSensors(scanID), 1, :);
            end
        end

        %%% set functions

        function setalgorithm(obj, options)
            arguments
                obj
                options.globalFusionRule (1, 1) string {mustBeMember(options.globalFusionRule, ["EGC", "MRC", "LDC", "GLDC", "BC", "CVBC", "NEGC", "SC", "PSC"])} = obj.globalFusionRule
                %%% LRT
                % EGC: Equal Gain Combining (LRT under equal SNR)
                % MRC: Maximal Ratio Combining (LRT under different SNR)
        
                %%% GLRT
                % LDC: Log-Divergence Combining (GLRT w/ joint MLE)
                % GLDC: Generalized Log-Divergence Combining (GLRT w/ independent MLE)
        
                %%% Binary
                % BC: Binary Combining (LRT under equal SNR)
                % CVBC: Chair-Varshney Binary Combining (LRT under different SNR)
        
                %%%
                % NEGC: Normalized Equal Gain Combining
                % SC: Selective Combining (Partial Selective Combining K = 1)
                % PSC: Partial Selective Combining (Combine largest K samples)
                options.numberOfBinaryDetectionsRule (1, 1) string {mustBeMember(options.numberOfBinaryDetectionsRule, ["and", "or", "majority", "fixedGlobalPFA", "userDefined"])} = obj.numberOfBinaryDetectionsRule
            end
            obj.globalFusionRule = options.globalFusionRule;
            obj.numberOfBinaryDetectionsRule = options.numberOfBinaryDetectionsRule;
        end

        function setuserdefinednumberofdetections(obj, options)
            arguments
                obj
                options.userDefinedNumberOfDetections (1, :) double {mustBeNonnegative, mustBeInteger} % [1 x M]
            end
            if isscalar(options.userDefinedNumberOfDetections)
                obj.numberOfBinaryDetections = repmat(options.userDefinedNumberOfDetections, 1, length(obj.numberOfSensors));
            elseif length(options.userDefinedNumberOfDetections) == length(obj.numberOfSensors)
                obj.numberOfBinaryDetections = options.userDefinedNumberOfDetections;
            else
                error('Length of "#detections for binary Combining algorithm" vector must be 1 or %d', length(obj.numberOfSensors));
            end
        end

        function setmontecarlo(obj, options)
            arguments
                obj
                options.numberOfTrials (1, 1) double {mustBeNonnegative, mustBeInteger} = obj.numberOfTrials
                options.seed (1, 1) double {mustBeNonnegative, mustBeInteger, mustBeInRange(options.seed, 0, 4294967295)} = obj.seed
            end
            obj.numberOfTrials = options.numberOfTrials;
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

        %%% simulation

        function simulate(obj, options)
            arguments
                obj
                options.simulationData (1, :) string {mustBeMember(options.simulationData, ["globalPFA", "globalPD"])} = ["globalPFA", "globalPD"]
                options.printStatus (1, 1) logical {mustBeMember(options.printStatus, [0, 1])} = true;
            end
            rng(obj.seed);
            Tglobal = obj.globalThreshold; % [NpfaGlobal x Nsnr x M x Nlocal]
            Tlocal = obj.localThreshold; % [1 x M] cell of [M x 1 x 1 x Nlocal] matrices
            numberOfSensorsMax = max(obj.numberOfSensors);
            n = obj.noise; % [M x 1 x Nmc] matrix
            if any(strcmp(options.simulationData, "globalPD"))
                s = obj.signal; % [M x 1] cell of [1 x Nsnr x Nmc] matrices
                x = cellfun(@(s) s + n(1 : size(s, 1), 1, :), s, 'UniformOutput', false); % [M x Nsnr x Nmc] matrix
            end
            for scanID = 1 : length(obj.numberOfSensors)
                if options.printStatus
                    fprintf('#sensors = %d/%d\n', obj.numberOfSensors(scanID), numberOfSensorsMax);
                end
                % no signal is present
                if any(strcmp(options.simulationData, "globalPFA"))
                    localTestStatisticsH0 = abs(n(1 : obj.numberOfSensors(scanID), 1, :)).^2; % [M x 1 x Nmc]
                    indicatorFunction = localTestStatisticsH0 > Tlocal{scanID}; % [M x 1 x Nmc x Nlocal]
                    globalTestStatisticsH0 = fusioncenter(localTestStatisticsH0);
                    obj.localPFAsimulation{scanID} = mean(indicatorFunction, 3); % [M x 1 x 1 x Nlocal]
                    obj.globalPFAsimulation(:, :, scanID, :) = mean(globalTestStatisticsH0 >= Tglobal(:, :, scanID, :), 3);
                    % equality is important for binary fusion with integer threshold
                end
                % signal + noise is present
                if any(strcmp(options.simulationData, "globalPD"))
                    localTestStatisticsH1 = abs(x{scanID}).^2; % [M x Nsr x Nmc]
                    indicatorFunction = localTestStatisticsH1 > Tlocal{scanID}; % [M x Nsr x Nmc x Nlocal]
                    globalTestStatisticsH1 = fusioncenter(localTestStatisticsH1);
                    obj.localPDsimulation{scanID} = mean(indicatorFunction, 3); % [M x Nsnr x 1 x Nlocal]
                    obj.globalPDsimulation(:, :, scanID, :) = mean(globalTestStatisticsH1 >= Tglobal(:, :, scanID, :), 3);
                    % equality is important for binary fusion with integer threshold
                end
            end
            function globalTestStatistics = fusioncenter(localTestStatistics)
                localData = indicatorFunction.*localTestStatistics; % [M x Nsnr x Nmc x Nlocal]
                switch obj.globalFusionRule
                %%% LRT
                    case "EGC" % Equal Gain Combining (LRT under equal SNR)
                        globalTestStatistics = sum(localData, 1); % [1 x Nsnr x Nmc x Nlocal]
                    case "MRC" % Maximal Ratio Combining (LRT under different SNR)
                        globalTestStatistics = sum(obj.fusionWeights{scanID}.*localData, 1) + sum(indicatorFunction.*obj.fusionBias{scanID}); % [1 x Nsnr x Nmc x Nlocal]
                 %%% GLRT
                    case "LDC" % Log-Divergence Combining (GLRT w/ joint MLE)
                        K = sum(indicatorFunction, 1);
                        T = sum(localData, 1)./K; % [1 x Nsnr x Nmc x Nlocal]
                        globalTestStatistics = K.*(T - log(T) - 1);
                    case "GLDC" % Generalized Log-Divergence Combining (GLRT w/ independent MLE)
                        globalTestStatistics = sum(localData - log(localData) - 1, 1);
                %%% Binary
                    case "BC" % Binary Combining (LRT under equal SNR)
                        globalTestStatistics = sum(indicatorFunction, 1); % [1 x Nsnr x Nmc x Nlocal]
                    case "CVBC" % Chair-Varshney Binary Combining (LRT under different SNR)
                        globalTestStatistics = sum(obj.fusionWeights{scanID}.*indicatorFunction, 1); % [1 x Nsnr x Nmc x Nlocal]
                %%%
                    case "NEGC" % Normalized Equal Gain Combining
                        globalTestStatistics = sum(localData, 1)./sum(indicatorFunction, 1); % [1 x Nsnr x Nmc x Nlocal]
                    case "SC" % Selective Combining (Partial Selective Combining K = 1)
                        globalTestStatistics = max(localTestStatistics, [], 1); % [1 x Nsnr x Nmc x Nlocal]
                    case "PSC" % Partial Selective Combining (Combine largest K samples)
                        localData = sort(localData, 1, 'descend');
                        globalTestStatistics = sum(localData(1 : min(obj.numberOfBinaryDetections(scanID), obj.numberOfSensors(scanID)), :, :, :), 1); % [1 x Nsnr x Nmc x Nlocal]
                end
            end
        end

        %%% utility functions

        function probability = totalCDF(obj, numberOfSamples, localThreshold, localProbability, fusionWeights, fusionBiases, snr)
            if nargin < 7
                snr = zeros(numberOfSamples, 1);
            end
            localThreshold = localThreshold(:);
            %%% ensure row vectors
            localProbability = localProbability(:).';
            fusionWeights = fusionWeights(:).';
            fusionBiases = fusionBiases(:).';
            snr = snr(:).';

            uniqueLocalThreshold = unique(localThreshold);
            uniqueSNR = unique(snr);
            uniqueRandomVariable = isscalar(uniqueSNR) && isscalar(uniqueLocalThreshold);

            switch obj.globalFusionRule
            %%% LRT
                case "EGC" % Equal Gain Combining (LRT under equal SNR)
                    if uniqueRandomVariable
                        probability = @(t) summationCCDF(t);
                    else
                        probability = @(t) generalizedSummationCCDF(t);
                    end
                case "MRC" % Maximal Ratio Combining (LRT under different SNR)
                    probability = @(t) generalizedSummationCCDF(t);
            %%% GLRT
                case "LDC" % Log-Divergence Combining (GLRT w/ joint MLE)
                    if uniqueRandomVariable
                        probability = @(t) summationCCDF(t);
                    else
                        probability = @(t) generalizedSummationCCDF(t);
                    end
                case "GLDC" % Generalized Log-Divergence Combining (GLRT w/ independent MLE)
                    error('not implemented');
                    % single sensor CDF
                    % CDF of sum is not known, find it derive it
                    % argument1 = max(-lambertw(0, -exp(-(t + 1)))/(1 + uniqueSNR), 0);
                    % argument2 = max(-lambertw(-1, -exp(-(t + 1)))/(1 + uniqueSNR), 0);
                    % probability = gammainc(argument1, 1, 'lower') + gammainc(argument2, 1, 'upper');
            %%% Binary
                case "BC" % Binary Combining (LRT under equal SNR)
                    probability = @(t) summationBinaryCCDF(t);
                case "CVBC" % Chair-Varshney Binary Combining (LRT under different SNR)
                    probability = @(t) generalizedSummationBinaryCCDF(t);
            %%%
                case "NEGC" % Normalized Equal Gain Combining
                    if uniqueRandomVariable
                        probability = @(t) summationCCDF(t);
                    else
                        probability = @(t) generalizedSummationCCDF(t);
                    end
                case "SC" % Selective Combining (Partial Selective Combining K = 1)
                    probability = @(t) orderStatisticsCCDF(t);
                case "PSC" % Partial Selective Combining (Combine largest K samples)
                    error('not implemented');
                    % K = obj.numberOfMaximumDetections(1);
                    % L = min(k, K);
                    % argument = max(t, L*uniqueLocalThreshold)/(1 + uniqueSNR);
                    % if argument < 0
                    %     probability = 1;
                    % else
                    %     if L < K
                    %         probability = gammainc(argument, L, 'upper');
                    %     else % order statistics (Simon & Alouini)
                    %         probability = 1 - sum(arrayfun(@(j) (-1).^j*nchoosek(k - K, j)*gammainc((t - L*uniqueLocalThreshold)/(1 + uniqueSNR), L + j, 'lower'), 0 : (k - K)));
                    %         % probability = probability*K*nchoosek(k, K);
                    %     end
                    % end
            end

            function probability = summationBinaryCCDF(t)
                pmf = poissonbinomialpdf(0 : min(floor(t), numberOfSamples));
                probability = 1 - sum(pmf);
            end

            function probability = generalizedSummationBinaryCCDF(t)
                if numberOfSamples <= 15
                    allSubsets = dec2bin(0 : 2^numberOfSamples - 1) - '0';
                    pmf = zeros(size(allSubsets, 1), 1);
                    weights = zeros(size(allSubsets, 1), 1);
                    for subSetID = 1 : size(allSubsets, 1)
                        indices = logical(allSubsets(subSetID, :));
                        pmf(subSetID) = prod(localProbability(indices))*prod(1 - localProbability(~indices));
                        weights(subSetID) = sum(fusionWeights(indices));
                    end
                    idx = weights <= t;
                    probability = 1 - sum(pmf(idx));
                else
                    error('#samples is huge for the calculation');
                end
            end

            function probability = summationCCDF(t)
                pmf = poissonbinomialpdf(0 : numberOfSamples);
                probability = 0;
                for k = 0 : numberOfSamples
                    subSetProbability = pmf(k + 1);
                    if k == 0
                        conditionalCDF = double(t <= 0);
                    else
                        switch obj.globalFusionRule
                            case "EGC" % Equal Gain Combining (LRT under equal SNR)
                                rate = 1/(1 + uniqueSNR);
                                shift = k*uniqueLocalThreshold;
                                conditionalCDF = gammaCCDF(t, rate, shift, k);
                            case "LDC" % Log-Divergence Combining (GLRT w/ joint MLE)
                                rate = k/(1 + uniqueSNR);
                                shift = uniqueLocalThreshold;
                                conditionalCDF = gammaKLdivCCDF(t, rate, shift, k);
                            case "NEGC" % Normalized Equal Gain Combining
                                rate = k/(1 + uniqueSNR);
                                shift = uniqueLocalThreshold;
                                conditionalCDF = gammaCCDF(t, rate, shift, k);
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
                                case "EGC" % Equal Gain Combining (LRT under equal SNR)
                                    rates = 1./(1 + snr(indices));
                                    shifts = sum(localThreshold(indices));
                                    conditionalCDF = hypoexponentialCCDF(t, rates, shifts);
                                case "MRC" % Maximal Ratio Combining (LRT under different SNR)
                                    rates = 1./fusionWeights(indices)./(1 + snr(indices));
                                    shifts = fusionWeights(indices)*localThreshold(indices) + sum(fusionBiases(indices));
                                    conditionalCDF = hypoexponentialCCDF(t, rates, shifts);
                                case "LDC" % Log-Divergence Combining (GLRT w/ joint MLE)
                                    K = nnz(indices);
                                    rates = K./(1 + snr(indices));
                                    shifts = sum(localThreshold(indices))/K;
                                    conditionalCDF = hypoexponentialKLdivCCDF(t/K, rates, shifts);
                                case "NEGC" % Normalized Equal Gain Combining
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
                    end
                end
            end

            %%%% probability mass functions

            function probability = poissonbinomialpdf(k)
                l = (0 : numberOfSamples).'/(numberOfSamples + 1); % [K x 1] vector
                phi = prod(1 + (exp(-1j*2*pi.*l) - 1).*localProbability, 2); % [K x 1] vector
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

            function probability = gammaKLdivCCDF(t, rate, shift, degree)
                % complementary CDF of Kullback-Leibler divergence of shifted gamma distribution
                argument1 = max(-lambertw(0, -exp(-(t/degree + 1))) - shift, 0);
                argument2 = max(-lambertw(-1, -exp(-(t/degree + 1))) - shift, 0);
                probability = gammainc(rate*argument1, degree, 'lower') + gammainc(rate*argument2, degree, 'upper');
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

            function probability = hypoexponentialKLdivCCDF(t, rates, shifts)
                % complementary CDF of Kullback-Leibler divergence of shifted hypoexponential distribution
                probability = 1 - hypoexponentialCCDF(-lambertw(0, -exp(-(t + 1))), rates, shifts) + hypoexponentialCCDF(-lambertw(-1, -exp(-(t + 1))), rates, shifts);
            end
        end

        %%% visualization

        function visualize(obj, options)
            arguments
                obj
                options.dataType (1, :) string {mustBeMember(options.dataType, ["analytical", "empirical"])} = ["analytical", "empirical"]
                options.x_axis (1, :) string {mustBeMember(options.x_axis, ["globalPFA", "numberOfSensors", "SNR", "globalThreshold", "localPFA", "localThreshold", "localPD"])} = "numberOfSensors"
                options.y_axis (1, 1) string {mustBeMember(options.y_axis, ["globalPFA", "globalPD"])} = "globalPD"
            end
            uniqueSNR = unique(cell2mat(cellfun(@(c) 10*log10(mean(10.^(.1*c), 1)), obj.SNR_input_dB, 'UniformOutput', false).'), 'rows', 'stable');
            uniqueLocalPFA = unique(permute(cell2mat(cellfun(@(c) mean(c, 1), obj.localPFA, 'UniformOutput', false).'), [1 4 2 3]), 'rows', 'stable'); % [1 x M] cell of [M x 1 x 1 x Nlocal] matrices
            for plotID = 1 : length(options.x_axis)
                plotFunc = @plot;
                legendString = '';
                legendTitleString = '';
                switch options.x_axis(plotID)
                    case "globalPFA"
                        plotFunc = @semilogx;
                        xData = obj.globalPFA; % [NpfaGlobal x 1]
                        if isscalar(xData) || max(abs(diff(xData))) < 10*eps
                            fprintf('choose non-scalar x-axis \n');
                            continue
                        end
                        xLabel = 'P_{FA}^{global}';
                        xTitle = xLabel;
                        indexingPriority = ["localPFA", "numberOfSensors", "SNR", "globalPFA"];
                        legendPriority = ["SNR", "numberOfSensors", "localPFA"];
                    case "numberOfSensors"
                        xData = obj.numberOfSensors; % [1 x M]
                        if isscalar(xData) || max(abs(diff(xData))) < 10*eps
                            fprintf('choose non-scalar x-axis \n');
                            continue
                        end
                        xLabel = '#sensors';
                        xTitle = xLabel;
                        indexingPriority = ["globalPFA", "localPFA", "SNR", "numberOfSensors"];
                        legendPriority = ["SNR", "localPFA", "globalPFA"];
                    case "SNR"
                        xData = uniqueSNR.'; % [1 x Nsnr]
                        if isscalar(xData)
                            fprintf('choose non-scalar x-axis \n');
                            continue
                        elseif isvector(xData)
                            if max(abs(diff(xData))) < 10*eps
                                fprintf('choose non-scalar x-axis \n');
                                continue
                            end
                        end
                        xLabel = 'SNR_{in} (dB)';
                        xTitle = xLabel;
                        indexingPriority = ["globalPFA", "localPFA", "numberOfSensors", "SNR"];
                        legendPriority = ["numberOfSensors", "localPFA", "globalPFA"];
                    case "globalThreshold"
                        xData = 10*log10(obj.globalThreshold); % [NpfaGlobal x 1 x M x Nlocal]
                        if isscalar(xData)
                            fprintf('choose non-scalar x-axis \n');
                            continue
                        elseif isvector(xData)
                            if max(abs(diff(xData))) < 10*eps
                                fprintf('choose non-scalar x-axis \n');
                                continue
                            end
                        end
                        xLabel = 'T_{global} (dB)';
                        xTitle = xLabel;
                        indexingPriority = ["globalPFA", "localPFA", "SNR", "numberOfSensors"];
                        legendPriority = ["SNR", "localPFA", "numberOfSensors", "globalPFA"];
                    case {"localPFA", "localThreshold", "localPD"}
                        switch options.x_axis(plotID)
                            case "localPFA"
                                plotFunc = @semilogx;
                                xData = uniqueLocalPFA; % [1 x Nlocal] or [M x Nlocal]
                                xLabel = 'P_{FA}^{local}';
                                if isscalar(xData)
                                    fprintf('choose non-scalar x-axis \n');
                                    continue
                                elseif isvector(xData)
                                    if max(abs(diff(xData))) < 10*eps
                                        fprintf('choose non-scalar x-axis \n');
                                        continue
                                    end
                                else
                                    xData = xData(1, :);
                                end
                                legendPriority = ["SNR", "numberOfSensors", "globalPFA"];
                            case "localThreshold"
                                uniqueLocalThreshold = unique(permute(cell2mat(cellfun(@(c) mean(c, 1), obj.localThreshold, 'UniformOutput', false).'), [1 4 2 3]), 'rows', 'stable'); % [1 x M] cell of [M x 1 x 1 x Nlocal] matrices
                                xData = uniqueLocalThreshold; % [1 x Nlocal] or [M x Nlocal]
                                xLabel = 'T_{local} (dB)';
                                if isscalar(xData)
                                    fprintf('choose non-scalar x-axis \n');
                                    continue
                                elseif isvector(xData)
                                    if max(abs(diff(xData))) < 10*eps
                                        fprintf('choose non-scalar x-axis \n');
                                        continue
                                    end
                                else
                                    xData = xData(1, :);
                                end
                                legendPriority = ["localPFA", "SNR", "numberOfSensors", "globalPFA"];
                            case "localPD"
                                uniqueLocalPD = unique(reshape(permute(cell2mat(cellfun(@(c) mean(c, 1), obj.localPD, 'UniformOutput', false).'), [1 2 4 3]), length(obj.numberOfSensors), []), 'rows', 'stable'); % [1 x M] cell of [M x Nsnr x 1 x Nlocal] matrices
                                xData = uniqueLocalPD;
                                if isscalar(xData)
                                    fprintf('choose non-scalar x-axis \n');
                                    continue
                                elseif isvector(xData)
                                    if max(abs(diff(xData))) < 10*eps
                                        fprintf('choose non-scalar x-axis \n');
                                        continue
                                    end
                                else
                                    xData = reshape(xData(1, :), [], length(uniqueLocalPFA(1, :))); % [Nsnr x Nlocal]
                                end
                                xLabel = 'P_{D}^{local}';
                                if isscalar(xData)
                                    fprintf('choose non-scalar x-axis \n');
                                    continue
                                end
                                legendPriority = ["localPFA", "SNR", "numberOfSensors", "globalPFA"];
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
                dataDimensions = ["globalPFA", "SNR", "numberOfSensors", "localPFA"];
                yDataEmpirical = obj.(options.y_axis + "simulation"); % [NpfaGlobal x Nsnr x M x Nlocal]
                if strcmp(options.y_axis, "globalPFA")
                    yDataAnalytical = obj.globalPFAanalytical;
                    if isequal(plotFunc, @semilogx)
                        plotFunc = @loglog;
                    else
                        plotFunc = @semilogy;
                    end
                    yLabel = 'P_{FA}^{global}';
                elseif strcmp(options.y_axis, "globalPD")
                    yDataAnalytical = obj.globalPD; % [NpfaGlobal x Nsnr x M x Nlocal]
                    yLabel = 'P_{D}^{global}';
                end
                yTitle = yLabel;
                dataSize = size(yDataAnalytical, [1 2 3 4]);
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
                    yDataEmpirical = squeeze(yDataEmpirical(indexing{:}));
                    yDataAnalytical = squeeze(yDataAnalytical(indexing{:}));
                    legendPriority = legendPriority(find(ismember(legendPriority, dataDimensions(cellfun(@(c) ~isscalar(c), indexing))), 1, 'first'));
                    switch options.x_axis(plotID)
                        case "globalThreshold"
                            indexing{2} = 1; % SNR
                            xData = squeeze(xData(indexing{:}));
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
                    yDataEmpirical = squeeze(yDataEmpirical);
                    yDataAnalytical = squeeze(yDataAnalytical);
                    switch options.x_axis(plotID)
                        case "globalThreshold"
                            xData = squeeze(xData);
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
                        case "SNR"
                            legendString = num2str(uniqueSNR(1, :).');
                            legendTitleString = 'SNR_{in} (dB)';
                        case "globalPFA"
                            legendString = num2str(obj.globalPFA);
                            legendTitleString = 'P_{FA}^{global}';
                        case "localPFA"
                            legendString = num2str(uniqueLocalPFA(1, :).');
                            legendTitleString = 'P_{FA}^{local}';
                    end
                end
                titleString = [yTitle, ' vs ', xTitle, indexingTitle];
                figure;
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
                ylim([0, 1]); xlim tight;
                grid off; grid minor; grid on;
                xlabel(xLabel); ylabel(yLabel);
                title(titleString, subTitleString);
                if ~isempty(legendString)
                    leg = legend(legendString, 'Location', 'best');
                    title(leg, legendTitleString);
                end
                hold off; drawnow;
            end

            function str = scinot(value)
                significand = 10^mod(log10(value), 1);
                exponent = floor(log10(value));
                if sign(exponent) == -1
                    str = sprintf('%0.2f•10^{-%i}', significand, abs(exponent));
                else
                    str = sprintf('%0.2f•10^%i', significand, exponent);
                end
            end
        end
    end
end