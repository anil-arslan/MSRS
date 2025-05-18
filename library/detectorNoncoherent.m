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
        globalFusionFunction function_handle
    end

    properties (SetAccess = private, GetAccess = public)
        globalPFA (:, 1) double {mustBeNonnegative, mustBeInRange(globalPFA, 0, 1)} = 1e-6 % [NpfaGlobal x 1]
        globalPFAsimulation double {mustBeNonnegative, mustBeInRange(globalPFAsimulation, 0, 1)} = [] % [NpfaGlobal x 1 x M x Nlocal]
    end

    properties (Dependent)
        globalThreshold double {mustBeNonnegative} % [NpfaGlobal x 1 x M x Nlocal]
        globalPD double {mustBeNonnegative, mustBeInRange(globalPD, 0, 1)} % [NpfaGlobal x Nsnr x M x Nlocal]
    end

    properties (SetAccess = private, GetAccess = public)
        globalPDsimulation double {mustBeNonnegative, mustBeInRange(globalPDsimulation, 0, 1)} = [] % [NpfaGlobal x Nsnr x M]
    end

    %%% Local detection parameters

    properties (SetAccess = private, GetAccess = public)
        constraint (1, :) string {mustBeMember(constraint, ["dataRate", "transmittedPower", "none"])} = ["dataRate", "transmittedPower"]
        localPFA (1, :) cell = {0.1} % [1 x M] cell of [M x Nlocal] matrices
    end

    properties (Dependent)
        localThreshold (1, :) cell % [1 x M] cell of [M x Nlocal] matrices
    end

    %%% Data generation

    properties (Dependent, Hidden)
        noise
        signal
    end

    properties (Access = private)
        numberOfSensorsLoopID (1, :) double {mustBeNonnegative, mustBeInteger} = 1
    end

    properties (Dependent, Access = private)
        currentNumberOfSensors (1, :) double {mustBeNonnegative, mustBeInteger}
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
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 2), obj.localPFA)));
            obj.globalPFAsimulation = zeros(numberOfScansPFAglobal, 1, numberOfScans, numberOfScansLocalPFA);
            obj.globalPDsimulation = zeros(numberOfScansPFAglobal, numberOfScansSNR, numberOfScans, numberOfScansLocalPFA);
        end

        %%% Fusion parameters

        function f = get.globalFusionFunction(obj)
            switch obj.globalFusionRule
                case "EGN"
                    f = @(x) sum(abs(x).^2, 1);
                case "MRC"
            end
        end

        function T = get.globalThreshold(obj)
            numberOfScansLocalPFA = unique((cellfun(@(c) size(c, 2), obj.localPFA)));
            numberOfScans = length(obj.numberOfSensors);
            T = zeros(length(obj.globalPFA), 1, length(obj.numberOfSensors), numberOfScansLocalPFA);
            for scanID = 1 : numberOfScans
                for PFAID = 1 : length(obj.globalPFA)
                    for localPFAID = 1 : numberOfScansLocalPFA
                        gamma = obj.globalPFA(PFAID);
                        M = obj.numberOfSensors(scanID);
                        p = unique(obj.localPFA{scanID}(:, localPFAID)); % unique local PFA across all sensors
                        lambda = unique(obj.localThreshold{scanID}(:, localPFAID)); % unique local PFA across all sensors
                        totalCDF = @(t) sum(arrayfun(@(k) obj.binomialcoef(M, k, p)*gammashifted(t, k, lambda), 0 : M));
                        T(PFAID, 1, scanID, localPFAID) = fzero(@(t) totalCDF(t) - (1 - gamma), [0, M*lambda + gammaincinv(gamma, M, 'upper')]);
                    end
                end
            end

            function cdfValue = gammashifted(t, k, lambda, theta)
                if nargin < 4
                    theta = 1;
                end
                % Shifted gamma CDF
                if k == 0
                    cdfValue = double(t >= 0); % Degenerate at T = 0
                    return;
                end
                if t < k*lambda
                    cdfValue = 0;
                else
                    cdfValue = 1 - gammainc((t - k*lambda).*theta, k, 'upper');
                end
            end
        end

        function ROC = get.globalPD(obj)
            % [NpfaGlobal x Nsnr x M x Nlocal]
            numberOfScansSNR = unique((cellfun(@(c) size(c, 2), obj.SNR_input_dB)));
            numberOfScans = length(obj.numberOfSensors);
            ROC = zeros(length(obj.globalPFA), numberOfScansSNR, numberOfScans);
            for scanID = 1 : numberOfScans
                uniqueSNR = 10.^(.1*unique(obj.SNR_input_dB{scanID}, 'rows'));
                for PFAID = 1 : length(obj.globalPFA)
                    if size(uniqueSNR, 1) % same SNR across all sensors
                        ROC(PFAID, :, scanID) = gammainc(obj.globalThreshold(PFAID, 1, scanID)./(1 + uniqueSNR), obj.numberOfSensors(scanID), 'upper');
                    else
                    end
                end
            end



            theta = 1./(1 + SNR);
            
            % Define the local threshold
            localThreshold = -log(PFA_local);
            
            % Local probability of detection
            PD_local = PFA_local.^theta;
            
            % Define the total CDF function F_T(gamma)
            globalPD = sum(arrayfun(@(m) binomial_term(M, m, PD_local)*conditional_cdf_gamma_shifted(globalThreshold, m, localThreshold, theta), 0 : M));



            function cdfValue = gammashifted(t, k, lambda, theta)
                if nargin < 4
                    theta = 1;
                end
                % Shifted gamma CDF
                if k == 0
                    cdfValue = double(t <= 0); % Degenerate at T = 0
                    return;
                end
                if t < k*lambda
                    cdfValue = 0;
                else
                    cdfValue = gammainc((t - k*lambda).*theta, k, 'upper');
                end
            end
        end

        %%% Local detection parameters

        function T = get.localThreshold(obj)
            T = cellfun(@(c) -log(c), obj.localPFA, 'UniformOutput', false);
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

        function M = get.currentNumberOfSensors(obj)
            M = obj.numberOfSensors(obj.numberOfSensorsLoopID);
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

        %%% simulation

        function simulate(obj, options)
            arguments
                obj
                options.simulationData (1, :) string {mustBeMember(options.simulationData, ["globalPFA", "globalPD"])} = ["globalPFA", "globalPD"]
                options.printStatus (1, 1) logical {mustBeMember(options.printStatus, [0, 1])} = true;
            end
            rng(obj.seed);
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
                    switch obj.globalFusionRule
                        case "EGN"
                            testStatisticsH0 = obj.globalFusionFunction(n(1 : obj.numberOfSensors(scanID), 1, :));
                        case "MRC"
                    end
                    obj.globalPFAsimulation(:, :, scanID) = mean(testStatisticsH0 > obj.globalThreshold(:, :, scanID), 3);
                end
                % signal + noise is present
                if any(strcmp(options.simulationData, "globalPD"))
                    switch obj.globalFusionRule
                        case "EGN"
                            testStatisticsH1 = obj.globalFusionFunction(x{scanID});
                        case "MRC"
                    end
                    obj.globalPDsimulation(:, :, scanID) = mean(testStatisticsH1 > obj.globalThreshold(:, :, scanID), 3);
                end
            end
        end

        function visualize(obj, options)
            arguments
                obj
                options.dataType (1, :) string {mustBeMember(options.dataType, ["analytical", "empirical"])} = ["analytical", "empirical"]
                options.x_axis (1, 1) string {mustBeMember(options.x_axis, ["globalPFA", "numberOfSensors", "SNR", "globalThreshold"])} = "numberOfSensors"
                options.y_axis (1, 1) string {mustBeMember(options.y_axis, ["globalPFA", "globalPD"])} = "globalPD"
            end
            plotFunc = @plot;
            legendString = '';
            legendTitleString = '';
            switch options.x_axis
                case "globalPFA"
                    plotFunc = @semilogx;
                    xData = obj.globalPFA; % [NpfaGlobal x 1]
                    if isscalar(xData)
                        fprintf('choose non-scalar x-axis \n');
                        return
                    end
                    xLabel = 'P_{FA}^{global}';
                    xTitle = xLabel;
                    indexingPriority = "numberOfSensors";
                    legendPriority = "SNR";
                    indexingTitle = [', #sensors = ', sprintf('%d', obj.numberOfSensors(ceil(end/2)))];
                case "numberOfSensors"
                    xData = obj.numberOfSensors; % [1 x M]
                    xLabel = '#sensors';
                    xTitle = xLabel;
                    indexingPriority = "globalPFA";
                    legendPriority = "SNR";
                    indexingTitle = [', P_{FA}^{global} = ', scinot(obj.globalPFA(1))];
                case "SNR"
                    uniqueSNR = unique(cell2mat(cellfun(@(c) mean(c, 1), obj.SNR_input_dB, 'UniformOutput', false).'), 'rows');
                    xData = uniqueSNR; % [1 x Nsnr]
                    xLabel = 'SNR_{in} (dB)';
                    xTitle = xLabel;
                    indexingPriority = "globalPFA";
                    legendPriority = "numberOfSensors";
                    indexingTitle = [', P_{FA}^{global} = ', scinot(obj.globalPFA(1))];
                case "globalThreshold"
                    xData = 10*log10(obj.globalThreshold); % [NpfaGlobal x 1 x M]
                    if isscalar(xData)
                        fprintf('choose non-scalar x-axis \n');
                        return
                    elseif ~isvector(xData)
                        xData = squeeze(xData(1, 1, :)); % [NpfaGlobal x 1 x M] --> [M x 1] after indexing
                    end
                    xLabel = 'T_{global} (dB)';
                    xTitle = xLabel;
                    indexingPriority = "globalPFA";
                    legendPriority = "SNR";
                    indexingTitle = [', P_{FA}^{global} = ', scinot(obj.globalPFA(1))];
            end
            switch options.y_axis
                case "globalPFA"
                    if isequal(plotFunc, @semilogx)
                        plotFunc = @loglog;
                    else
                        plotFunc = @semilogy;
                    end
                    yDataEmpirical = obj.globalPFAsimulation; % [NpfaGlobal x 1 x M]
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
                            yDataAnalytical = repmat(obj.globalPFA, 1, length(uniqueSNR)); % [NpfaGlobal x Nsnr]
                            yDataEmpirical = repmat(yDataEmpirical, 1, length(uniqueSNR)); % [NpfaGlobal x Nsnr x M]
                            yDataEmpirical = yDataEmpirical(:, :, ceil(end/2));
                            indexingTitle = [', #sensors = ', sprintf('%d', obj.numberOfSensors(ceil(end/2)))];
                            legendString = num2str(obj.globalPFA);
                            legendTitleString = 'P_{FA}^{global}';
                    end
                    yLabel = 'P_{FA}^{global}';
                    yTitle = yLabel;
                case "globalPD"
                    dataDimensions = ["globalPFA", "SNR", "numberOfSensors"];
                    yDataEmpirical = obj.globalPDsimulation; % [NpfaGlobal x Nsnr x M]
                    yDataAnalytical = obj.globalPD; % [NpfaGlobal x Nsnr x M]
                    yLabel = 'P_{D}^{global}';
                    yTitle = yLabel;

                    dataSize = size(yDataEmpirical);

                    indexing = arrayfun(@(i) 1 : dataSize(i), 1 : length(dataSize), 'UniformOutput', false);
                    indexingDimension = dataDimensions == indexingPriority;
                    indexing{indexingDimension} = 1;
                    yDataEmpirical = squeeze(yDataEmpirical(indexing{:}));
                    yDataAnalytical = squeeze(yDataAnalytical(indexing{:}));
                    switch legendPriority
                        case "numberOfSensors"
                            legendString = num2str(obj.numberOfSensors.');
                            legendTitleString = '#sensors';
                        case "SNR"
                            uniqueSNR = unique(cell2mat(cellfun(@(c) mean(c, 1), obj.SNR_input_dB, 'UniformOutput', false).'), 'rows');
                            if size(uniqueSNR, 1) == 1
                                legendString = num2str(uniqueSNR.');
                                legendTitleString = 'SNR_{in} (dB)';
                            end
                    end
            end
            titleString = [yTitle, ' vs ', xTitle, indexingTitle];
            subTitleString = ['#trials = ', scinot(obj.numberOfTrials)];
            figure;
            if any(strcmp(options.dataType, "empirical"))
                plotFunc(xData, yDataEmpirical, '-', 'LineWidth', 2);
                hold on;
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
        function b = binomialcoef(m, k, p)
            b = nchoosek(m, k)*p^k*(1 - p)^(m - k);
        end
    end
end