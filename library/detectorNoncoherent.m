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
        globalFusionRule (1, 1) string {mustBeMember(globalFusionRule, ["EGN", "MRC"])} = "EGN"
        % EGN: Equal Gain Combining
        % MRC: Maximal Ratio Combining
    end

    properties (Dependent)
        globalFusionFunction (1, 1) function_handle
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

        function f = get.globalFusionFunction(obj)
            switch obj.globalFusionRule
                case "EGN"
                    f = @(x) sum(x, 1);
                case "MRC"
            end
        end

        function Tglobal = get.globalThreshold(obj)
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
                        lambda = unique(Tlocal{scanID}(:, localPFAID)); % unique local threshold across all sensors
                        if lambda
                            F = obj.totalCDF(M, lambda, pfaLocal);
                            Tglobal(PFAID, 1, scanID, localPFAID) = fzero(@(t) F(t) - pfaGlobal, [0, M*lambda + gammaincinv(pfaGlobal, M, 'upper')]);
                        else
                            Tglobal(PFAID, 1, scanID, localPFAID) = gammaincinv(pfaGlobal, M, 'upper');
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
                        lambda = unique(Tlocal{scanID}(:, localPFAID)); % unique local threshold across all sensors
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
                        lambda = unique(Tlocal{scanID}(:, 1, 1, localPFAID)); % unique local threshold across all sensors
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

        function setmontecarlo(obj, options)
            arguments
                obj
                options.numberOfTrials (1, 1) double {mustBeNonnegative, mustBeInteger} = obj.numberOfTrials
                options.seed (1, 1) double {mustBeNonnegative, mustBeInteger, mustBeInRange(options.seed, 0, 4294967295)} = obj.seed
            end
            obj.numberOfTrials = options.numberOfTrials;
            obj.seed = options.seed;
        end

        function setconstraint(obj, options)
            arguments
                obj
                options.constraint (1, :) string {mustBeMember(options.constraint, ["dataRate", "transmittedPower"])} = ["dataRate", "transmittedPower"]
            end
            if any(strcmp(options.constraint, "transmittedPower"))
                for scanID = 2 : length(obj.numberOfSensors)
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
                            localData = localTestStatisticsH0; % [M x 1 x Nmc]
                        case "binary"
                            localData = indicatorFunction; % [M x 1 x Nmc x Nlocal]
                    end
                    switch obj.globalFusionRule
                        case "EGN"
                            dataSize = size(indicatorFunction, [2 3 4]);
                            globalTestStatisticsH0 = zeros([1, dataSize]);
                            for idx1 = 1 : dataSize(1)
                                for idx2 = 1 : dataSize(2)
                                    for idx3 = 1 : dataSize(3)
                                        globalTestStatisticsH0(:, idx1, idx2, idx3) = obj.globalFusionFunction(localData(indicatorFunction(:, idx1, idx2, idx3), idx1, idx2)); % [1 x 1 x Nmc x Nlocal]
                                    end
                                end
                            end
                        case "MRC"
                    end
                    obj.globalPFAsimulation(:, :, scanID, :) = mean(globalTestStatisticsH0 > Tglobal(:, :, scanID), 3);
                end
                % signal + noise is present
                if any(strcmp(options.simulationData, "globalPD"))
                    localTestStatisticsH1 = abs(x{scanID}).^2; % [M x Nsr x Nmc]
                    indicatorFunction = localTestStatisticsH1 > Tlocal{scanID}; % [M x Nsr x Nmc x Nlocal]
                    obj.localPDsimulation{scanID} = mean(indicatorFunction, 3); % [M x Nsnr x 1 x Nlocal]
                    switch obj.localDataSent
                        case "power"
                            localData = localTestStatisticsH1; % [M x Nsnr x Nmc]
                        case "binary"
                            localData = indicatorFunction; % [M x Nsnr x Nmc x Nlocal]
                    end
                    switch obj.globalFusionRule
                        case "EGN"
                            dataSize = size(indicatorFunction, [2 3 4]);
                            globalTestStatisticsH1 = zeros([1, dataSize]);
                            for idx1 = 1 : dataSize(1)
                                for idx2 = 1 : dataSize(2)
                                    for idx3 = 1 : dataSize(3)
                                        globalTestStatisticsH1(:, idx1, idx2, idx3) = obj.globalFusionFunction(localData(indicatorFunction(:, idx1, idx2, idx3), idx1, idx2)); % [1 x Nsnr x Nmc x Nlocal]
                                    end
                                end
                            end
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
                options.x_axis (1, 1) string {mustBeMember(options.x_axis, ["globalPFA", "numberOfSensors", "SNR", "globalThreshold", "localPFA", "localThreshold", "localPD"])} = "numberOfSensors"
                options.y_axis (1, 1) string {mustBeMember(options.y_axis, ["globalPFA", "globalPD"])} = "globalPD"
            end
            plotFunc = @plot;
            legendString = '';
            legendTitleString = '';
            switch options.x_axis
                case "globalPFA"
                    plotFunc = @semilogx;
                    xData = obj.globalPFA; % [NpfaGlobal x 1]
                    if isscalar(xData) || max(abs(diff(xData))) < eps
                        fprintf('choose non-scalar x-axis \n');
                        return
                    end
                    xLabel = 'P_{FA}^{global}';
                    xTitle = xLabel;
                    indexingPriority = ["numberOfSensors", "SNR"];
                    legendPriority = "SNR";
                    indexingTitle = [', #sensors = ', sprintf('%d', obj.numberOfSensors(ceil(end/2)))];
                case "numberOfSensors"
                    xData = obj.numberOfSensors; % [1 x M]
                    if isscalar(xData) || max(abs(diff(xData))) < eps
                        fprintf('choose non-scalar x-axis \n');
                        return
                    end
                    xLabel = '#sensors';
                    xTitle = xLabel;
                    indexingPriority = ["globalPFA", "localPFA"];
                    legendPriority = "SNR";
                    uniqueLocalPFA = unique(permute(cell2mat(cellfun(@(c) mean(c, 1), obj.localPFA, 'UniformOutput', false).'), [1 4 2 3]), 'rows');
                    indexingTitle = [', P_{FA}^{global} = ', scinot(obj.globalPFA(1)), ', P_{FA}^{local} = ', scinot(uniqueLocalPFA(1))];
                case "SNR"
                    uniqueSNR = unique(cell2mat(cellfun(@(c) mean(c, 1), obj.SNR_input_dB, 'UniformOutput', false).'), 'rows');
                    xData = uniqueSNR; % [1 x Nsnr]
                    if isscalar(xData) || max(abs(diff(xData))) < eps
                        fprintf('choose non-scalar x-axis \n');
                        return
                    end
                    xLabel = 'SNR_{in} (dB)';
                    xTitle = xLabel;
                    indexingPriority = ["globalPFA", "localPFA"];
                    legendPriority = "numberOfSensors";
                    uniqueLocalPFA = unique(permute(cell2mat(cellfun(@(c) mean(c, 1), obj.localPFA, 'UniformOutput', false).'), [1 4 2 3]), 'rows');
                    indexingTitle = [', P_{FA}^{global} = ', scinot(obj.globalPFA(1)), ', P_{FA}^{local} = ', scinot(uniqueLocalPFA(1))];
                case "globalThreshold"
                    xData = 10*log10(obj.globalThreshold); % [NpfaGlobal x 1 x M x Nlocal]
                    if isscalar(xData)
                        fprintf('choose non-scalar x-axis \n');
                        return
                    elseif ~isvector(xData)
                        xData = squeeze(xData(1, 1, :, 1)); % [NpfaGlobal x 1 x M x Nlocal] --> [M x 1] after indexing
                    end
                    if max(abs(diff(xData))) < eps
                        fprintf('choose non-scalar x-axis \n');
                        return
                    end
                    xLabel = 'T_{global} (dB)';
                    xTitle = xLabel;
                    indexingPriority = ["globalPFA", "localPFA"];
                    legendPriority = "SNR";
                    uniqueLocalPFA = unique(permute(cell2mat(cellfun(@(c) mean(c, 1), obj.localPFA, 'UniformOutput', false).'), [1 4 2 3]), 'rows');
                    indexingTitle = [', P_{FA}^{global} = ', scinot(obj.globalPFA(1)), ', P_{FA}^{local} = ', scinot(uniqueLocalPFA(1))];
                case {"localPFA", "localThreshold", "localPD"}
                    switch options.x_axis
                        case "localPFA"
                            plotFunc = @semilogx;
                            uniqueLocalPFA = unique(permute(cell2mat(cellfun(@(c) mean(c, 1), obj.localPFA, 'UniformOutput', false).'), [1 4 2 3]), 'rows', 'stable'); % [1 x M] cell of [M x 1 x 1 x Nlocal] matrices
                            xData = uniqueLocalPFA; % [1 x Nlocal]
                            xLabel = 'P_{FA}^{local}';
                            if isscalar(xData)
                                fprintf('choose non-scalar x-axis \n');
                                return
                            elseif isvector(xData)
                                if max(abs(diff(xData))) < eps
                                    fprintf('choose non-scalar x-axis \n');
                                    return
                                end
                            end
                        case "localThreshold"
                            uniqueLocalThreshold = unique(permute(cell2mat(cellfun(@(c) mean(c, 1), obj.localThreshold, 'UniformOutput', false).'), [1 4 2 3]), 'rows', 'stable'); % [1 x M] cell of [M x 1 x 1 x Nlocal] matrices
                            xData = uniqueLocalThreshold; % [1 x Nlocal]
                            xLabel = 'T_{local} (dB)';
                            if isscalar(xData)
                                fprintf('choose non-scalar x-axis \n');
                                return
                            elseif isvector(xData)
                                if max(abs(diff(xData))) < eps
                                    fprintf('choose non-scalar x-axis \n');
                                    return
                                end
                            end
                        case "localPD"
                            uniqueLocalPD = unique(reshape(permute(cell2mat(cellfun(@(c) mean(c, 1), obj.localPD, 'UniformOutput', false).'), [1 2 4 3]), length(obj.numberOfSensors), []), 'rows', 'stable'); % [1 x M] cell of [M x Nsnr x 1 x Nlocal] matrices
                            xData = reshape(uniqueLocalPD(1, :), [], length(obj.localPFA)); % [Nsnr x Nlocal]
                            xLabel = 'P_{D}^{local}';
                            if isscalar(xData)
                                fprintf('choose non-scalar x-axis \n');
                                return
                            end
                    end
                    xTitle = xLabel;
                    indexingPriority = ["globalPFA", "numberOfSensors", "SNR"];
                    legendPriority = ["SNR", "numberOfSensors", "globalPFA"];
                    indexingTitle = [', P_{FA}^{global} = ', scinot(obj.globalPFA(1)), ', #sensors = ', sprintf('%d', obj.numberOfSensors(ceil(end/2)))];
            end
            switch options.y_axis
                case "globalPFA"
                    if isequal(plotFunc, @semilogx)
                        plotFunc = @loglog;
                    else
                        plotFunc = @semilogy;
                    end
                    % dataDimensions = ["globalPFA", "", "numberOfSensors", "localPFA"];
                    yDataEmpirical = obj.globalPFAsimulation; % [NpfaGlobal x 1 x M x Nlocal]
                    indexingTitle = '';
                    switch options.x_axis
                        case "globalPFA"
                            yDataAnalytical = repmat(obj.globalPFA, 1, length(obj.numberOfSensors)); % [NpfaGlobal x M]
                            legendString = num2str(obj.numberOfSensors.');
                            legendTitleString = '#sensors';
                        case "numberOfSensors"
                            yDataAnalytical = repmat(obj.globalPFA, 1, length(obj.numberOfSensors)); % [NpfaGlobal x M]
                            legendString = num2str(obj.globalPFA);
                            legendTitleString = 'P_{FA}^{global}';
                        case "SNR"
                            yDataAnalytical = repmat(obj.globalPFA, 1, length(xData)); % [NpfaGlobal x Nsnr]
                            yDataEmpirical = repmat(yDataEmpirical, 1, length(xData)); % [NpfaGlobal x Nsnr x M]
                            yDataEmpirical = yDataEmpirical(:, :, ceil(end/2));
                            indexingTitle = [', #sensors = ', sprintf('%d', obj.numberOfSensors(ceil(end/2)))];
                            legendString = num2str(obj.globalPFA);
                            legendTitleString = 'P_{FA}^{global}';
                        case "localPFA"
                            yDataAnalytical = repmat(obj.globalPFA, 1, length(xData)); % [NpfaGlobal x Nlocal]
                            yDataEmpirical = squeeze(yDataEmpirical(1, 1, :, :));
                            indexingTitle = [', P_{FA}^{global} = ', sprintf('%d', obj.globalPFA(1))];
                            legendString = num2str(obj.numberOfSensors.');
                            legendTitleString = '#sensors';
                    end
                    yLabel = 'P_{FA}^{global}';
                    yTitle = yLabel;
                case "globalPD"
                    dataDimensions = ["globalPFA", "SNR", "numberOfSensors", "localPFA"];
                    yDataEmpirical = obj.globalPDsimulation; % [NpfaGlobal x Nsnr x M x Nlocal]
                    yDataAnalytical = obj.globalPD; % [NpfaGlobal x Nsnr x M x Nlocal]
                    yLabel = 'P_{D}^{global}';
                    yTitle = yLabel;
                    dataSize = size(yDataAnalytical, [1 2 3 4]);
                    indexing = arrayfun(@(i) 1 : dataSize(i), 1 : length(dataSize), 'UniformOutput', false);
                    indexingDimension = ismember(dataDimensions, indexingPriority);
                    if nnz(dataSize(indexingDimension) ~= 1) > 2
                        indices = ones(1, 4);
                        indices(3) = ceil(dataSize(3)/2); % # sensors
                        possibleIndicedDimensions = dataDimensions(indexingDimension);
                        dataDimension = length(dataDimensions);
                        indexingDimension = false(1, dataDimension);
                        for priorityID = 1 : dataDimension - 2
                            indexingDimension(possibleIndicedDimensions == indexingPriority(priorityID)) = true;
                        end
                        indexing(indexingDimension) = num2cell(indices(indexingDimension));
                        yDataEmpirical = squeeze(yDataEmpirical(indexing{:}));
                        yDataAnalytical = squeeze(yDataAnalytical(indexing{:}));
                        legendPriority = legendPriority(1);
                    else
                        yDataDimensions = dataSize ~= 1;
                        yDataEmpirical = squeeze(yDataEmpirical);
                        yDataAnalytical = squeeze(yDataAnalytical);
                        legendPriority = legendPriority(find(ismember(legendPriority, dataDimensions(yDataDimensions)), 1, 'first'));
                        switch legendPriority
                            case "numberOfSensors"
                                uniqueSNR = unique(cell2mat(cellfun(@(c) mean(c, 1), obj.SNR_input_dB, 'UniformOutput', false).'), 'rows', 'stable');
                                indexingTitle = [', P_{FA}^{global} = ', scinot(obj.globalPFA(1)), ', SNR = ', sprintf('%g', uniqueSNR(1)), ' dB'];
                            case "globalPFA"
                                uniqueSNR = unique(cell2mat(cellfun(@(c) mean(c, 1), obj.SNR_input_dB, 'UniformOutput', false).'), 'rows', 'stable');
                                indexingTitle = [', SNR = ', sprintf('%g', uniqueSNR(1)), ' dB'];
                        end
                    end
                    if ~isvector(yDataEmpirical) && ~isvector(yDataAnalytical)
                        switch legendPriority
                            case "numberOfSensors"
                                legendString = num2str(obj.numberOfSensors');
                                legendTitleString = '#sensors';
                            case "SNR"
                                uniqueSNR = unique(cell2mat(cellfun(@(c) mean(c, 1), obj.SNR_input_dB, 'UniformOutput', false).'), 'rows', 'stable');
                                if size(uniqueSNR, 1) == 1
                                    legendString = num2str(uniqueSNR.');
                                    legendTitleString = 'SNR_{in} (dB)';
                                end
                            case "globalPFA"
                                legendString = num2str(obj.globalPFA);
                                legendTitleString = 'P_{FA}^{global}';
                        end
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
                    plotFunc(xData, yDataAnalytical, '-');
                end
            end
            grid off; grid minor; grid on;
            xlabel(xLabel); ylabel(yLabel);
            title(titleString, subTitleString);
            if ~isempty(legendString)
                leg = legend(legendString, 'Location', 'best');
                title(leg, legendTitleString);
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

    methods (Static)
        function probability = totalCDF(numberOfSamples, localThreshold, localPD, snr)
            if nargin < 4
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
                    if t < k*localThreshold
                        probability = 1;
                    else
                        probability = gammainc((t - k*localThreshold)./(1 + snr), k, 'upper');
                    end
                end
            end
        end
    end
end