classdef detectorNoncoherent < handle
    %detector Summary of this class goes here
    %   Detailed explanation goes here

    %%% Model parameters

    properties (SetAccess = private, GetAccess = public)
        seed (1, 1) double {mustBeNonnegative, mustBeInteger, mustBeInRange(seed, 0, 4294967295)} = 0
        numberOfTrials (1, 1) double {mustBeNonnegative, mustBeInteger} = 1e3
        numberOfSensors (1, :) double {mustBeNonnegative, mustBeInteger} = 10 % [1 x M]
        SNR_input_dB (1, :) cell = {4} % [1 x M] cell of [M x Nsnr] matrices
    end

    %%% Fusion parameters

    properties (SetAccess = private, GetAccess = public)
        globalFusionRule (1, 1) string {mustBeMember(globalFusionRule, ["EGC", "NEGC", "MRC"])} = "EGC"
        % EGC: Equal Gain Combining
        % NEGC: Normalized Equal Gain Combining
        % MRC: Maximal Ratio Combining
    end

    properties (SetAccess = private, GetAccess = public)
        globalPFA (:, 1) double {mustBeNonnegative, mustBeInRange(globalPFA, 0, 1)} = 1e-6 % [NpfaGlobal x 1]
        globalPFAsimulation double {mustBeNonnegative, mustBeInRange(globalPFAsimulation, 0, 1)} = [] % [NpfaGlobal x 1 x M x Nlocal]
    end

    properties (Dependent)
        globalThreshold double {mustBeNonnegative} % [NpfaGlobal x 1 x M x Nlocal]
        globalPFAanalytical double {mustBeNonnegative} % [NpfaGlobal x 1 x M x Nlocal]
        globalPD double {mustBeNonnegative, mustBeInRange(globalPD, 0, 1)} % [NpfaGlobal x Nsnr x M x Nlocal]
    end

    properties (SetAccess = private, GetAccess = public)
        globalPDsimulation double {mustBeNonnegative, mustBeInRange(globalPDsimulation, 0, 1)} = [] % [NpfaGlobal x Nsnr x M x Nlocal]
    end

    %%% Local detection parameters

    properties (SetAccess = private, GetAccess = public)
        localDataSent (1, 1) string {mustBeMember(localDataSent, ["power", "binary"])} = "power"
    end

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

        %%% Fusion parameters

        function Tglobal = get.globalThreshold(obj)
            % [NpfaGlobal x 1 x M x Nlocal]
            Tlocal = obj.localThreshold;
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            Tglobal = zeros(length(obj.globalPFA), 1, length(obj.numberOfSensors), numberOfScansLocalPFA);
            for scanID = 1 : numberOfScans
                M = obj.numberOfSensors(scanID);
                for PFAID = 1 : length(obj.globalPFA)
                    pfaGlobal = obj.globalPFA(PFAID);
                    for localPFAID = 1 : numberOfScansLocalPFA
                        pfaLocal = unique(obj.localPFA{scanID}(:, :, :, localPFAID)); % unique local PFA across all sensors
                        lambda = unique(Tlocal{scanID}(:, :, :, localPFAID)); % unique local threshold across all sensors
                        thresholdNSC = gammaincinv(pfaGlobal, M, 'upper'); % equal gain non-selective combiner threshold, conventional
                        thresholdSearchSpace = [0, M*lambda + thresholdNSC];
                        if lambda
                            F = obj.totalCDF(M, lambda, pfaLocal);
                            Tglobal(PFAID, 1, scanID, localPFAID) = fzero(@(t) F(t) - pfaGlobal, thresholdSearchSpace);
                        else
                            switch obj.globalFusionRule
                                case "NEGC"
                                    thresholdNSC = thresholdNSC/M;
                            end
                            Tglobal(PFAID, 1, scanID, localPFAID) = thresholdNSC;
                        end
                    end
                end
            end
        end
        
        function pfa = get.globalPFAanalytical(obj)
            % [NpfaGlobal x Nsnr x M x Nlocal]
            Tglobal = obj.globalThreshold;
            Tlocal = obj.localThreshold;
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            pfa = zeros(length(obj.globalPFA), 1, length(obj.numberOfSensors), numberOfScansLocalPFA);
            for scanID = 1 : numberOfScans
                M = obj.numberOfSensors(scanID);
                for PFAID = 1 : length(obj.globalPFA)
                    for localPFAID = 1 : numberOfScansLocalPFA
                        pfaLocal = unique(obj.localPFA{scanID}(:, :, :, localPFAID)); % unique local PFA across all sensors
                        lambda = unique(Tlocal{scanID}(:, :, :, localPFAID)); % unique local threshold across all sensors
                        F = obj.totalCDF(M, lambda, pfaLocal);
                        pfa(PFAID, 1, scanID, localPFAID) = F(Tglobal(PFAID, 1, scanID, localPFAID));
                    end
                end
            end
        end

        function ROC = get.globalPD(obj)
            % [NpfaGlobal x Nsnr x M x Nlocal]
            Tglobal = obj.globalThreshold;
            Tlocal = obj.localThreshold;
            pDlocal = obj.localPD;
            numberOfScansSNR = unique((cellfun(@(c) size(c, 2), obj.SNR_input_dB)));
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 4), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            ROC = zeros(length(obj.globalPFA), numberOfScansSNR, numberOfScans, numberOfScansLocalPFA);
            for scanID = 1 : numberOfScans
                M = obj.numberOfSensors(scanID);
                uniqueSNR = 10.^(.1*unique(obj.SNR_input_dB{scanID}, 'rows', 'stable'));
                for PFAID = 1 : length(obj.globalPFA)
                    for localPFAID = 1 : numberOfScansLocalPFA
                        gamma = Tglobal(PFAID, 1, scanID, localPFAID);
                        lambda = unique(Tlocal{scanID}(:, :, :, localPFAID)); % unique local threshold across all sensors
                        if size(uniqueSNR, 1) % same SNR across all sensors
                            for snrID = 1 : numberOfScansSNR
                                pdLocal = unique(pDlocal{scanID}(:, snrID, 1, localPFAID)); % unique local PD across all sensors
                                F = obj.totalCDF(M, lambda, pdLocal, uniqueSNR(snrID));
                                ROC(PFAID, snrID, scanID, localPFAID) = F(gamma);
                            end
                        else
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
                options.globalFusionRule (1, 1) string {mustBeMember(options.globalFusionRule, ["EGC", "NEGC", "MRC"])} = obj.globalFusionRule
                % EGC: Equal Gain Combining
                % NEGC: Normalized Equal Gain Combining
                % MRC: Maximal Ratio Combining
                options.localDataSent (1, 1) string {mustBeMember(options.localDataSent, ["power", "binary"])} = obj.localDataSent
            end
            obj.globalFusionRule = options.globalFusionRule;
            obj.localDataSent = options.localDataSent;
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
            Tglobal = obj.globalThreshold; % [NpfaGlobal x 1 x M x Nlocal]
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
                    obj.localPFAsimulation{scanID} = mean(indicatorFunction, 3); % [M x 1 x 1 x Nlocal]
                    switch obj.localDataSent
                        case "power"
                            localData = indicatorFunction.*localTestStatisticsH0; % [M x 1 x Nmc x Nloca]
                        case "binary"
                            localData = indicatorFunction; % [M x 1 x Nmc x Nlocal]
                    end
                    switch obj.globalFusionRule
                        case "EGC"
                            globalTestStatisticsH0 = sum(localData, 1); % [1 x 1 x Nmc x Nlocal]
                        case "NEGC"
                            globalTestStatisticsH0 = sum(localData, 1)./sum(indicatorFunction, 1); % [1 x 1 x Nmc x Nlocal]
                        case "MRC"
                    end
                    obj.globalPFAsimulation(:, :, scanID, :) = mean(globalTestStatisticsH0 > Tglobal(:, :, scanID, :), 3);
                end
                % signal + noise is present
                if any(strcmp(options.simulationData, "globalPD"))
                    localTestStatisticsH1 = abs(x{scanID}).^2; % [M x Nsr x Nmc]
                    indicatorFunction = localTestStatisticsH1 > Tlocal{scanID}; % [M x Nsr x Nmc x Nlocal]
                    obj.localPDsimulation{scanID} = mean(indicatorFunction, 3); % [M x Nsnr x 1 x Nlocal]
                    switch obj.localDataSent
                        case "power"
                            localData = indicatorFunction.*localTestStatisticsH1; % [M x Nsnr x Nmc x Nloca]
                        case "binary"
                            localData = indicatorFunction; % [M x Nsnr x Nmc x Nlocal]
                    end
                    switch obj.globalFusionRule
                        case "EGC"
                            globalTestStatisticsH1 = sum(localData, 1); % [1 x Nsnr x Nmc x Nlocal]
                        case "NEGC"
                            globalTestStatisticsH1 = sum(localData, 1)./sum(indicatorFunction, 1); % [1 x Nsnr x Nmc x Nlocal]
                        case "MRC"
                    end
                    obj.globalPDsimulation(:, :, scanID, :) = mean(globalTestStatisticsH1 > Tglobal(:, :, scanID, :), 3);
                end
            end
        end

        function visualize(obj, options)
            arguments
                obj
                options.dataType (1, :) string {mustBeMember(options.dataType, ["analytical", "empirical"])} = ["analytical", "empirical"]
                options.x_axis (1, :) string {mustBeMember(options.x_axis, ["globalPFA", "numberOfSensors", "SNR", "globalThreshold", "localPFA", "localThreshold", "localPD"])} = "numberOfSensors"
                options.y_axis (1, 1) string {mustBeMember(options.y_axis, ["globalPFA", "globalPD"])} = "globalPD"
            end
            uniqueSNR = unique(cell2mat(cellfun(@(c) mean(c, 1), obj.SNR_input_dB, 'UniformOutput', false).'), 'rows', 'stable');
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
                    % if strcmp(options.y_axis, "globalPD") && ~any(strcmp(options.dataType, "analytical"))
                    %     maxValues = max(yDataEmpirical, [], find(size(yDataEmpirical) ~= length(xData)));
                    %     plotFunc(xData, maxValues, '-ok');
                    % end
                else
                    subTitleString = '';
                end
                if any(strcmp(options.dataType, "analytical"))
                    if any(strcmp(options.dataType, "empirical"))
                        plotFunc(xData, yDataAnalytical, '--k');
                    else
                        plotFunc(xData, yDataAnalytical, '-', 'LineWidth', 2);
                        hold on;
                        if strcmp(options.y_axis, "globalPD")
                            [maxValues, maxIndices] = max(yDataAnalytical, [], find(size(yDataAnalytical) == length(xData)));
                            boundaries = maxIndices == 1 | maxIndices == length(xData);
                            plotFunc(xData(maxIndices(~boundaries)), maxValues(~boundaries), '-ok');
                        end
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

        function probability = totalCDF(obj, numberOfSamples, localThreshold, localPD, snr)
            if nargin < 5
                snr = 0;
            end
            probability = @(t) sum(arrayfun(@(k) binopdf(k, numberOfSamples, localPD)*gammaincshifted(t, k, localThreshold, snr), 0 : numberOfSamples));
            function probability = gammaincshifted(t, k, localThreshold, snr)
                % conditional CDF
                if nargin < 4
                    snr = 0;
                end
                if k == 0
                    probability = double(t <= 0);
                else
                    switch obj.globalFusionRule
                        case "EGC"
                            if t < k*localThreshold
                                probability = 1;
                            else
                                probability = gammainc((t - k*localThreshold)./(1 + snr), k, 'upper');
                            end
                        case "NEGC"
                            if t < localThreshold
                                probability = 1;
                            else
                                probability = gammainc((k*t - k*localThreshold)./(1 + snr), k, 'upper');
                            end
                    end
                end
            end
        end
    end
end