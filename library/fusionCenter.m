classdef fusionCenter < handle
    %fusionCenter Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = private, GetAccess = public)
        interfaces (1, :) interface = interface.empty()
        monteCarlo (1, 1) struct = struct( ...
            'seed', 0, ...
            'numberOfTrials', 1, ...
            'numberOfTrialsParallel', 1)
        configuration (1, 1) struct = struct( ...
            'globalPFA', 1e-6, ...
            'localPFA', 1, ...
            'numberOfTrainingCells', 20, ...
            'numberOfGuardCells', 5);
        configurationMonostatic (1, 1) struct = struct( ...
            'removeBoundaryDetectionDF', 1, ...
            'removeBoundaryDetectionMF', 0, ...
            'interpolationDF', 1, ...
            'interpolationMF', 0, ...
            'groupingCloseness', 2, ... 
            'CFAR', 1);
        gridResolution (3, 1) double {mustBeNonnegative} = 100*ones(3, 1)
        globalIntegrationRule (1, 1) string {mustBeMember(globalIntegrationRule, ["SLC", "WSLC", "BC", "CVBC", "LLC", "WLLC", "EGC", "MRC"])} = "SLC"
        %%% LRT: widely seperated receivers | spatially i.i.d circularly symmetric complex gaussian signal model
            % SLC: Square Law Combining (LRT under equal SNR)
            % WSLC: Weighted Sqaure Law Combining (LRT under different SNR)

        %%% Binary
            % BC: Binary Combining (LRT under equal SNR)
            % CVBC: Chair-Varshney Binary Combining (LRT under different SNR)

        %%% Other
            % LLC: Linear Law Combining
            % WLLC: Weighted Linear Law Combining

        %%% LRT: phase synchronous closely spaced receivers | spatially fully correlated circularly symmetric complex gaussian signal model
            % EGC: Equal Gain Combining
            % MRC: Matched Filter
    end

    properties (Dependent)
        network (1, 1) radarNetwork
        gridPoints cell
        gridSize (1, 3) double {mustBeNonnegative}
        gridPointsMesh struct
    end

    properties (Dependent)
        averageSNR_dB double % dB power, % (Ntx x Nrx x Nt x Nmcp matrix)
        averageSNR_lin dobule % linear scale power, % (Ntx x Nrx x Nt x Nmcp matrix)
        noisePowersPerSample_W double {mustBePositive} % linear scale noise power
        dictionaryComplexWeights double % (Ntx x Nrx x Ni matrix)
    end

    properties (SetAccess = private, GetAccess = public)
        hypothesizedVisibleInterface interface
        blindZone logical {mustBeNumericOrLogical, mustBeMember(blindZone, [0, 1])} = []; 
        expectedCellSNRs double % (Ntx x Nrx x Ni matrix)
        dictionaryCellMapping = [] % (Ni x Nrx x Ntx matrix)
    end

    properties (SetAccess = private, GetAccess = public)
        signalsMatchFiltered cell = {} % (1 x Nrx cell of Ns + L - 1 x Nrxch x Ntxch x Nmcp matrix) Ncmp : number of parallel trials
        thresholdCFAR cell = {} % (1 x Nrx cell of Ns + L - 1 x Nrxch x Nmcp matrix)
    end

    properties (SetAccess = private, GetAccess = public)
        detectionReport struct = struct( ...
            'cellIDs', [], ...
            'detectedRawPower', [], ...
            'detectedProcessedPower', [], ...
            'estimatedPositions', [], ...
            'estimationError', [], ...
            'localDetectionNodes', [], ...
            'numberOfUtilizedNodes', [], ...
            'integratedPower', []);
        coverageSimulationReport struct = struct( ...
            'targetCellIDs', [], ...
            'globalPDrealized', [], ...
            'globalPDmodel', [], ...
            'globalPFArealized', [], ...
            'globalPFAmodel', [], ...
            'meanSNRsModel', [], ...
            'meanSNRsRealized', [])
        detectionSimulationReport struct = struct( ...
            'globalPDrealized', [], ...
            'globalPDmodeled', [], ...
            'globalPFArealized', [], ...
            'globalPFAmodel', [], ...
            'meanSNRsRealized', [], ...
            'meanSNRsModel', [], ...
            'numberOfUtilizedNodes', []);
    end

    properties (Access = private)
        seedShuffle (1, 1) logical {mustBeNumericOrLogical, mustBeMember(seedShuffle, [0, 1])} = 0
    end

    properties (Constant)
        speedOfLight (1, 1) double {mustBePositive} = physconst('LightSpeed');
    end

    methods
        function obj = fusionCenter(options)
            arguments
                options.interfaces (1, 1) interface
                options.gridResolution (3, 1) double {mustBeNonnegative} = zeros(3, 1)
            end
            %fusionCenter Construct an instance of this class
            %   Detailed explanation goes here
            obj.interfaces = options.interfaces;
            if ~any(options.gridResolution)
                obj.gridResolution = obj.interfaces.network.defaultGridResolution;
            else
                obj.gridResolution = options.gridResolution;
            end
            obj.interfaces.numberOfTrialsParallel = obj.monteCarlo.numberOfTrialsParallel;
        end

        %%% get methods

        function mc = get.monteCarlo(obj)
            mc.numberOfTrials = obj.monteCarlo.numberOfTrials;
            mc.numberOfTrialsParallel = obj.monteCarlo.numberOfTrialsParallel;
            if obj.seedShuffle
                mc.seed = 'shuffle';
            else
                mc.seed = obj.monteCarlo.seed;
            end
        end

        function net = get.network(obj)
            net = obj.interfaces.network;
        end

        function SNRin = get.averageSNR_lin(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            SNRin = 10.^(.1*obj.averageSNR_dB);
        end

        function SNRin = get.averageSNR_dB(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            SNRin = obj.interfaces.averageSNR_dB + obj.network.processingGain_dB.';
        end

        function cfg = get.configuration(obj)
            cfg = obj.configuration;
            cfg.detector = detector( ...
                "globalPFA", cfg.globalPFA, ...
                "SNR_input_dB", reshape(obj.averageSNR_dB(:, :, :, 1), obj.network.numberOfActiveBistaticPairs, []), ...
                "numberOfSensors", obj.network.numberOfActiveBistaticPairs, ...
                "localPFA", cfg.localPFA);
            switch obj.interfaces.spatialCoherency
                case "deterministic"
                    signalAmplitudeModel = "deterministic";
                    signalPhaseModel = "correlatedUniform";
                case "coherent"
                    signalAmplitudeModel = "correlatedExponential";
                    signalPhaseModel = "correlatedUniform";
                case "correlatedAmplitudeFluctuations"
                    signalAmplitudeModel = "correlatedExponential";
                    signalPhaseModel = "decorrelatedUniform";
                case "noncoherent"
                    signalAmplitudeModel = "decorrelatedExponential";
                    signalPhaseModel = "decorrelatedUniform";
            end
            cfg.detector.setalgorithm( ...
                "binaryDetectionRule", "notSpecified", ...
                "binaryDetectionConstraint", "fixedGlobalPFA", ...
                "binaryDetectionConstraint", "fixedGlobal|LocalPFA", ...
                "signalAmplitudeModel", signalAmplitudeModel, ...
                "signalPhaseModel", signalPhaseModel, ...
                "globalFusionRule", obj.globalIntegrationRule ...
                );
        end

        function n = get.noisePowersPerSample_W(obj)
            n = 10.^(.1*[obj.network.activeReceivingNodes.noisePowerPerSample_dB]); % kTBN
        end

        function w = get.dictionaryComplexWeights(obj)
            % (Ntx x Nrx x Ni matrix)
            switch obj.globalIntegrationRule
                case "EGC"
                    w = obj.expectedCellSNRs./abs(obj.expectedCellSNRs); % phase only
                case "MRC"
                    w = obj.expectedCellSNRs;
                case "WSLC"
                    w = abs(obj.expectedCellSNRs).^2./(1 + abs(obj.expectedCellSNRs).^2);
                case "WLLC"
                    w = abs(obj.expectedCellSNRs);
                otherwise
                    w = ones(obj.network.numberOfActiveTransmittingNodes, obj.network.numberOfActiveReceivingNodes, size(obj.dictionaryCellMapping, 1));
            end
        end

        function posScan = get.gridPoints(obj)
            networkBoundary = obj.network.boundaryListened;
            posScan = cell(3, 1);
            for i = 1 : 3
                boundaryCenter = mean(networkBoundary(i, :));
                posScan{i} = [flip(boundaryCenter : -obj.gridResolution(i) : networkBoundary(i, 1)) boundaryCenter + obj.gridResolution(i) : obj.gridResolution(i) : networkBoundary(i, 2)];
                if isscalar(posScan{i})
                    posScan{i} = mean(networkBoundary(i, :));
                end
            end
        end

        function size = get.gridSize(obj)
            scanPoints = obj.gridPoints;
            size = [numel(scanPoints{1}), numel(scanPoints{2}), numel(scanPoints{3})];
        end

        function gridScan = get.gridPointsMesh(obj)
            scanPoints = obj.gridPoints;
            [gridScan.x, gridScan.y, gridScan.z] = meshgrid(scanPoints{1}, scanPoints{2}, scanPoints{3});
        end

        %%% monostatic

        function configuremonostatic(obj, options)
            arguments
                obj
                options.removeBoundaryDetectionDF (1, 1) double {mustBeNonnegative} = obj.configurationMonostatic.removeBoundaryDetectionDF
                options.removeBoundaryDetectionMF (1, 1) double {mustBeNonnegative} = obj.configurationMonostatic.removeBoundaryDetectionMF
                options.interpolationDF (1, 1) double {mustBeNonnegative} = obj.configurationMonostatic.interpolationDF
                options.interpolationMF (1, 1) double {mustBeNonnegative} = obj.configurationMonostatic.interpolationMF
                options.groupingCloseness (1, 1) double {mustBePositive} = obj.configurationMonostatic.groupingCloseness
                options.CFAR (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.CFAR, [0, 1])} = obj.configurationMonostatic.CFAR
            end
            obj.configurationMonostatic.removeBoundaryDetectionDF = options.removeBoundaryDetectionDF;
            obj.configurationMonostatic.removeBoundaryDetectionMF = options.removeBoundaryDetectionMF;
            obj.configurationMonostatic.interpolationDF = options.interpolationDF;
            obj.configurationMonostatic.interpolationMF = options.interpolationMF;
            obj.configurationMonostatic.groupingCloseness = options.groupingCloseness;
            obj.configurationMonostatic.CFAR = options.CFAR;
        end

        %%% set methods

        function setdictionary(obj)
            % (Ns_i x Ni x Nrx x Ntx matrix)
            gridScan = obj.gridPointsMesh;
            hypothesizedPositions = reshape(permute(cat(4, gridScan.x, gridScan.y, gridScan.z), [4 1 2 3]), [3 prod(obj.gridSize)]);
            targets = target( ...
                "position", hypothesizedPositions, ...
                "meanRCS_dbsm", 0);
            hypothesizedFullInterface = interface( ...
                'network', obj.network, ...
                'targets', targets);
            %%% add back of array and outside of the CPI to blind zone
            averageSNRdB = hypothesizedFullInterface.averageSNR_dB;
            obj.blindZone = squeeze(all(isinf(permute(averageSNRdB, [5 2 3 1 4])), [1 2 4]));
            visibleZone = ~obj.blindZone;
            numberOfVisibleCells = nnz(visibleZone);
            obj.expectedCellSNRs = 10.^(.5*averageSNRdB(:, :, visibleZone)); % (Ntx x Nrx x Ni matrix)
            targets = target( ...
                "position", hypothesizedPositions(:, visibleZone), ...
                "meanRCS_dbsm", 0);
            obj.hypothesizedVisibleInterface = interface( ...
                'network', obj.network, ...
                'targets', targets);
            obj.hypothesizedVisibleInterface.configure('pathLoss', 0);
            signals = obj.hypothesizedVisibleInterface.signalReceivedFromScatterers; % Ns x 1 x Nt x Ntx x M
            % has the information of relative phases as well
            % pos = hypothesizedPositions(:, (all(abs(squeeze(signals{1}(:, :, :, :, 1))) == 0)));
            % figure; plot(pos(1, :), pos(2, :), '.');
            matchFilters = obj.network.matchFilter; % L x Nrx x Ntx matrix
            switch obj.network.networkMode
                case 'multiStatic'
                    % (Ni x Nrx x Ntx matrix)
                    obj.dictionaryCellMapping = zeros(numberOfVisibleCells, obj.network.numberOfActiveReceivingNodes, obj.network.numberOfActiveTransmittingNodes);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                            matchFilter = conj(matchFilters(:, rxID, txID));
                            for cellID = 1 : numberOfVisibleCells
                                [~, obj.dictionaryCellMapping(cellID, rxID, txID)] = max(abs(conv(signals{rxID}(:, 1, cellID, txID), matchFilter)));
                            end
                        end
                    end
                case 'monoStatic'
                    obj.dictionaryCellMapping = cell(1, obj.network.numberOfActiveReceivingNodes);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        matchFilter = conj(matchFilters(:, rxID, obj.network.monoStaticTransmitterIDs(rxID)));
                        obj.dictionaryCellMapping{rxID} = zeros(numberOfVisibleCells, obj.network.activeReceivingNodes(rxID).numberOfTotalChannels);
                        for channelID = 1 : obj.network.activeReceivingNodes(rxID).numberOfTotalChannels
                            for cellID = 1 : numberOfVisibleCells
                                [~, obj.dictionaryCellMapping{rxID}(cellID, channelID)] = max(abs(conv(signals{rxID}(:, 1, cellID, obj.network.monoStaticTransmitterIDs(rxID), channelID), matchFilter)));
                            end
                        end
                    end
            end
        end

        function configure(obj, options)
            arguments
                obj
                options.globalPFA (1, 1) double {mustBeNonnegative} = obj.configuration.globalPFA
                options.localPFA (1, 1) double {mustBeNonnegative} = obj.configuration.localPFA
                options.globalIntegrationRule (1, 1) string {mustBeMember(options.globalIntegrationRule, ["SLC", "WSLC", "BC", "CVBC", "LLC", "WLLC", "EGC", "MRC"])} = "SLC"
                    %%% LRT: widely seperated receivers | spatially i.i.d circularly symmetric complex gaussian signal model
                    % SLC: Square Law Combining (LRT under equal SNR)
                    % WSLC: Weighted Sqaure Law Combining (LRT under different SNR)
    
                    %%% Binary
                    % BC: Binary Combining (LRT under equal SNR)
                    % CVBC: Chair-Varshney Binary Combining (LRT under different SNR)
    
                    %%% Other
                    % LLC: Linear Law Combining
                    % WLLC: Weighted Linear Law Combining
    
                    %%% LRT: phase synchronous closely spaced receivers | spatially fully correlated circularly symmetric complex gaussian signal model
                    % EGC: Equal Gain Combining
                    % MRC: Matched Filter
                options.numberOfTrainingCells (1, 1) {mustBePositive, mustBeInteger} = obj.configuration.numberOfTrainingCells
                options.numberOfGuardCells (1, 1) {mustBeNonnegative, mustBeInteger} = obj.configuration.numberOfGuardCells
                options.numberOfTrials (1, 1) {mustBeNonnegative, mustBeInteger} = obj.monteCarlo.numberOfTrials
                options.numberOfTrialsParallel (1, 1) {mustBeNonnegative, mustBeInteger} = obj.monteCarlo.numberOfTrialsParallel
                options.seed (1, 1) {mustBeNonnegative, mustBeInteger, mustBeLessThan(options.seed, 4294967296)} = 0
                options.seedShuffle (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.seedShuffle, [0, 1])} = obj.seedShuffle
            end
            obj.configuration.globalPFA = options.globalPFA;
            obj.configuration.localPFA = options.localPFA;
            obj.configuration.numberOfTrainingCells = options.numberOfTrainingCells;
            obj.configuration.numberOfGuardCells = options.numberOfGuardCells;
            obj.globalIntegrationRule = options.globalIntegrationRule;
            mc = obj.monteCarlo;
            mc.numberOfTrials = options.numberOfTrials;
            mc.numberOfTrialsParallel = options.numberOfTrialsParallel;
            obj.interfaces.numberOfTrialsParallel = options.numberOfTrialsParallel;
            mc.seed = options.seed;
            obj.monteCarlo = mc;
            obj.seedShuffle = options.seedShuffle;
            rng(obj.monteCarlo.seed);
        end

        %%% utility methods

        function cellID = neighbourssquarewindow(obj, cellID, options)
            arguments
                obj
                cellID (1, 1) {mustBePositive, mustBeInteger}
                options.offset (1, 1) {mustBeNonnegative, mustBeInteger} = 1
            end
            if options.offset
                [xInd, yInd, zInd] = ind2sub(obj.gridSize, cellID);
                xNeigbours = validneighbours(xInd, 1);
                yNeigbours = validneighbours(yInd, 2);
                zNeigbours = validneighbours(zInd, 3);
                [xNeigbours, yNeigbours, zNeigbours] = meshgrid(xNeigbours, yNeigbours, zNeigbours);
                cellID = sub2ind(obj.gridSize, xNeigbours(:), yNeigbours(:), zNeigbours(:));
            end
            function ngbhrs = validneighbours(center, index)
                ngbhrs = [center, center + (1 : options.offset), center - (1 : options.offset)];
                ngbhrs = ngbhrs(ngbhrs > 0 & ngbhrs <= obj.gridSize(index));
            end
        end

        function applyspatialprocessing(obj, options)
            arguments
                obj
                options.saveSpectrum (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.saveSpectrum, [0, 1])} = 0
                options.doesPrint (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.doesPrint, [0, 1])} = 0
            end
            if isempty(obj.dictionaryCellMapping)
                fprintf('dictionary is not set\n');
                return;
            end
            % (Ni x Nmcp matrix)
            mc = obj.monteCarlo;
            config = obj.configuration;
            gridScan = obj.gridPointsMesh;
            visibleZone = ~obj.blindZone;
            L = obj.network.pulseWidthSample;
            Ns = size([obj.network.activeReceivingNodes.samplingInstants], 1);
            Nmf = Ns + L - 1;
            matchFilters = obj.network.matchFilter; % L x Nrx x Ntx matrix
            NTC = config.numberOfTrainingCells; % CA-CFAR
            NGC = config.numberOfGuardCells; % CA-CFAR
            signalBeamformed = obj.interfaces.signalBeamformed;
            report = struct( ...
                'cellIDs', [], ...
                'detectedRawPower', [], ...
                'detectedProcessedPower', [], ...
                'estimatedPositions', [], ...
                'estimationError', [], ...
                'localDetectionNodes', [], ...
                'numberOfUtilizedNodes', [], ...
                'integratedPower', []);
            report = repmat(report, 1, mc.numberOfTrialsParallel);
            %%% local processing START
            cellMapping = obj.dictionaryCellMapping; % (Ni x Nrxch matrix)
            switch obj.network.networkMode
                case 'multiStatic'
                    Nrx = obj.network.numberOfActiveReceivingNodes;
                    Ntxrx = obj.network.numberOfActiveTransmittingNodes;
                    dict = obj.dictionaryComplexWeights; % (Ntx x Nrx x Ni matrix)
                    dict = permute(dict, [3 4 2 1]);
                case 'monoStatic'
                    Nrx = obj.network.numberOfParallelSignalFusions;
                    Ntxrx = 1;
                    dict = 1;
            end
            numOfSpatialChannels = Nrx*Ntxrx;
            observations = cell(1, Nrx);
            obj.signalsMatchFiltered = cell(1, Nrx);
            for rxID = 1 : Nrx
                numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
                observations{rxID} = signalBeamformed{rxID}./sqrt(obj.noisePowersPerSample_W(rxID)); % Ns x Nrxch x Ntxch x Nmcp matrix
                observations{rxID}(isinf(observations{rxID})) = 0;
                obj.signalsMatchFiltered{rxID} = zeros(Nmf(rxID), numberOfTotalChannels, Ntxrx, mc.numberOfTrialsParallel);
                obj.thresholdCFAR{rxID} = zeros(Nmf(rxID), numberOfTotalChannels, Ntxrx, mc.numberOfTrialsParallel);
                for txID = 1 : Ntxrx
                    switch obj.network.networkMode
                        case 'multiStatic'
                            matchFilter = conj(matchFilters(:, rxID, txID));
                        case 'monoStatic'
                            matchFilter = conj(matchFilters(:, rxID, obj.network.monoStaticTransmitterIDs(rxID)));
                            % if isempty(obj.network.activeReceivingNodes(rxID).directionFinder.table)
                            %     obj.network.setdirectionfinders;
                            % end
                            % DF = obj.network.activeReceivingNodes(rxID).directionFinder.table;
                            % azimuthDF = obj.network.activeReceivingNodes(rxID).directionFinder.azimuth;
                            % azimuthResolution = obj.network.activeReceivingNodes(rxID).directionFinder.azimuthResolution;
                            % elevationSteering = obj.network.activeReceivingNodes(rxID).array.steeringElevation;
                            % azimuthSteering = obj.network.activeReceivingNodes(rxID).array.steeringAzimuth;
                            % sizeDF = size(DF, 1);
                            % rotationMat = obj.network.activeReceivingNodes(rxID).array.rotationMatrix;
                            % origin = obj.network.activeReceivingNodes(rxID).position;
                    end
                    for chID = 1 : numberOfTotalChannels
                        for mcID = 1 : mc.numberOfTrialsParallel
                            matchFilteredSignal = conv(observations{rxID}(:, chID, txID, mcID), matchFilter);
                            % localThreshold = config.detector.localThreshold{1}(rxID); % constant local threshold
                            %%% CA-CFAR
                            numberOfSamples = size(matchFilteredSignal, 1);
                            localThreshold = zeros(numberOfSamples, 1);
                            for sampleID = 1 : numberOfSamples % CA-CFAR
                                trainingStart = sampleID + NGC + 1;
                                upperTraining = unique(min(trainingStart : trainingStart + NTC - 1, numberOfSamples));
                                trainingStart = sampleID - NGC - 1;
                                lowerTraining = unique(max(trainingStart - NTC + 1 : trainingStart, 1));
                                trainingCellIDs = [lowerTraining, upperTraining];
                                numberOfTrainingCells = length(trainingCellIDs);
                                alpha = numberOfTrainingCells.*(config.localPFA.^(-1/numberOfTrainingCells) - 1);
                                localThreshold(sampleID) = alpha.*mean(abs(matchFilteredSignal(trainingCellIDs)).^2);
                            end
                            %%% CA-CFAR
                            obj.thresholdCFAR{rxID}(:, chID, txID, mcID) = localThreshold;
                            obj.signalsMatchFiltered{rxID}(:, chID, txID, mcID) = matchFilteredSignal; % matched filter ciktisini kirpmak gerekecek
                            % figure; plot(20*log10(abs(matchFilteredSignal)));
                            % hold on; plot(10*log10(preThresholdCFAR));
                        end
                    end
                end
            end
            observations = cell2mat(obj.signalsMatchFiltered); % [(Ns x 1 x Nrxch x Ntx)]
            binaryIndicator = abs(observations).^2 > cell2mat(obj.thresholdCFAR); % logical Ns x 1



            % detectionSamples = mod(detectionIndex - 1, Ns) + 1;
            % if obj.configurationMonostatic.removeBoundaryDetectionMF
            %     boundaryDetections = detectionSamples ~= 1 & detectionSamples ~= Ns;
            %     detectionIndex = detectionIndex(boundaryDetections);
            %     detectionSamples = detectionSamples(boundaryDetections);
            % end
            % numberOfDetections = length(detectionIndex);
            % if numberOfDetections
            %     detectionPower = 20*log10(absObs(detectionIndex, mcID));
            %     idx = 1;
            %     while idx < numberOfDetections
            %         group = abs(detectionSamples - detectionSamples(idx)) < obj.configurationMonostatic.groupingCloseness;
            %         groupPowers = detectionPower(group);
            %         groupSamples = detectionSamples(group);
            %         groupIndices = detectionIndex(group);
            %         [powerRepresenter, groupRepresenterIndex] = max(groupPowers);
            %         indexRepresenter = groupIndices(groupRepresenterIndex);
            %         sampleRepresenter = groupSamples(groupRepresenterIndex);
            %         detectionPower = [powerRepresenter; detectionPower(~group)];
            %         detectionIndex = [indexRepresenter; detectionIndex(~group)];
            %         detectionSamples = [sampleRepresenter; detectionSamples(~group)];
            %         idx = idx + 1 - sum(group(idx - 1 : -1 : 1));
            %         numberOfDetections = length(detectionIndex);
            %     end
            %     if obj.configurationMonostatic.interpolationMF
            %         detectionSamplesInterpolated = zeros(numberOfDetections, 1);
            %         interpolatables = detectionSamples ~= 1 & detectionSamples ~= Ns;
            %         rightNgbh = absObs(detectionIndex(interpolatables) + 1, mcID);
            %         leftNgbh = absObs(detectionIndex(interpolatables) - 1, mcID);
            %         center = absObs(detectionIndex(interpolatables), mcID);
            %         detectionSamplesInterpolated(interpolatables) = detectionSamples(interpolatables) + (leftNgbh - rightNgbh)./(2*center - (rightNgbh + leftNgbh))/2;
            %     else
            %         detectionSamplesInterpolated = detectionSamples;
            %     end
            %     spectrum = abs(DF*obs(detectionSamples, :, mcID)');
            %     [maxVal, idxDF] = max(spectrum);
            %     if obj.configurationMonostatic.removeBoundaryDetectionDF
            %         boundaryDetections = idxDF == 1 | idxDF == length(azimuthDF);
            %         maxVal = maxVal(~boundaryDetections);
            %         idxDF = idxDF(~boundaryDetections);
            %         detectionPower = detectionPower(~boundaryDetections);
            %         detectionSamplesInterpolated = detectionSamplesInterpolated(~boundaryDetections);
            %         numberOfDetections = length(idxDF);
            %     end
            %     if numberOfDetections
            %         elevation = zeros(1, length(idxDF));
            %         azimuth = azimuthDF(idxDF);
            %         if obj.configurationMonostatic.interpolationDF
            %             interpolatables = idxDF ~= 1 & idxDF ~= length(azimuthDF);
            %             specVec = sizeDF*(0 : numberOfDetections - 1);
            %             rightNgbh = spectrum(idxDF(interpolatables) + 1 + specVec(interpolatables));
            %             leftNgbh = spectrum(idxDF(interpolatables) - 1 + specVec(interpolatables));
            %             corrections = (leftNgbh - rightNgbh)./(2*maxVal(interpolatables) - (rightNgbh + leftNgbh))/2;
            %             azimuth(interpolatables) = azimuth(interpolatables) + azimuthResolution.*corrections;
            %         end
            %         range = obj.speedOfLight.*Ts*(detectionSamplesInterpolated.' - L - 1)/2;
            %         elevation = asind(sind(elevation) + sind(elevationSteering));
            %         azimuth = asind(sind(azimuth) + sind(azimuthSteering));
            %         imaginaryDetections = imag(azimuth) | imag(elevation) | imag(range);
            %         range = range(~imaginaryDetections);
            %         elevation = elevation(~imaginaryDetections);
            %         azimuth = azimuth(~imaginaryDetections);
            %         numberOfDetections = length(range);
            %         if numberOfDetections
            %             [x, y, z] = sph2cart(pi*azimuth/180, pi*elevation/180, range);
            %             position = [x; y; z];
            %             detection(rxID, mcID).power = detectionPower;
            %             detection(rxID, mcID).range = range;
            %             detection(rxID, mcID).elevation = elevation;
            %             detection(rxID, mcID).azimuth = azimuth;
            %             detection(rxID, mcID).numberOfDetections = numberOfDetections;
            %             detection(rxID, mcID).position = rotationMat*position + origin;
            %         end
            %     end
            % end




            observations = observations.*binaryIndicator;
            switch obj.globalIntegrationRule
                case {"SLC", "WSLC", "BC", "CVBC"}
                    observations = abs(observations).^2; % noncoherent integration
            end
            %%% local processing END

            if iscell(cellMapping)
                cellMapping = cellMapping{1};
            end
            Nc = size(cellMapping, 1); % number of cells
            for mcID = 1 : mc.numberOfTrialsParallel
                report(mcID).detectedRawPower = nan(numOfSpatialChannels, 1);
                observationTrial = permute(observations(:, :, :, mcID), [1 4 2 3]); % (Ns x 1 x Nrxch x Ntxch) or (Ns*Nrxch*Ntxch x 1)
                mappedSignal = zeros(Nc, numOfSpatialChannels);
                for spchID = 1 : numOfSpatialChannels
                    mappedSignal(:, spchID) = observationTrial(cellMapping(:, spchID), :, spchID);
                end
                switch obj.globalIntegrationRule
                    case "EGC"
                        report(mcID).integratedPower = abs(sum(mappedSignal, 2)).^2; % coherent integration
                    case "MRC"
                        report(mcID).integratedPower = abs(sum(dict.*mappedSignal, 2)).^2; % coherent integration
                    case "SLC"
                        report(mcID).integratedPower = sum(mappedSignal, 2); % noncoherent integration
                    case "WSLC"
                        report(mcID).integratedPower = sum(dict.*mappedSignal, 2); % noncoherent integration
                    case "BC"
                        report(mcID).integratedPower = sum(mappedSignal ~= 0, 2); % binary integration
                end

                cellIDs = find(report(mcID).integratedPower > config.detector.globalThreshold);
                [~, idx] = max(report(mcID).integratedPower(cellIDs));

                %%% global processing END

                cellIDs = cellIDs(idx);
                numberOfEstimations = length(cellIDs);

                report(mcID).detectedProcessedPower = report(mcID).integratedPower(cellIDs).';
                report(mcID).detectedRawPower = mappedSignal(cellIDs, :).';
                report(mcID).localDetectionNodes = false(Nrx*Ntxrx, numberOfEstimations);
                for estimationID = 1 : numberOfEstimations
                    report(mcID).localDetectionNodes(:, estimationID) = logical(mappedSignal(cellIDs(estimationID), :));
                    report(mcID).numberOfUtilizedNodes(estimationID) = nnz(report(mcID).localDetectionNodes(:, estimationID));
                end
                if options.doesPrint && mcID == 1
                    fprintf('peak power = %g dB, minimum power exceeding threshold = %g dB, Threshold = %g dB\n', 10*log10(max(report(mcID).detectedProcessedPower)), 10*log10(min(report(mcID).detectedProcessedPower)), 10*log10(config.detector.globalThreshold));
                end
                visibleRegion = false(Nc, 1);
                visibleRegion(cellIDs) = true;
                wholeRegion = false(prod(obj.gridSize), 1);
                wholeRegion(visibleZone) = visibleRegion;
                cellIDsWholeRegion = find(wholeRegion);
                report(mcID).cellIDs = cellIDsWholeRegion.';
                report(mcID).estimatedPositions = [gridScan.x(cellIDsWholeRegion) gridScan.y(cellIDsWholeRegion) gridScan.z(cellIDsWholeRegion)].';
                report(mcID).estimationError = zeros(3, numberOfEstimations, obj.interfaces.numberOfTargets);
                if ~isempty(obj.interfaces.targets)
                    for targetID = 1 : obj.interfaces.numberOfTargets
                        report(mcID).estimationError(:, :, targetID) = obj.interfaces.targets.position(:, targetID) - report(mcID).estimatedPositions;
                    end
                end
            end
            
            obj.detectionReport = report;
        end

        %%% monte carlo simulation

        function simulatedetection(obj, options)
            arguments
                obj
                options.randomOnCell (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.randomOnCell, [0, 1])} = 0
                options.numberOfTargets (1, 1) double {mustBePositive} = inf
            end
            mc = obj.monteCarlo;
            config = obj.configuration;
            numberOfTotalTrials = mc.numberOfTrials*mc.numberOfTrialsParallel;
            previousDetectionReport = obj.detectionReport;
            cleanup = onCleanup(@() cleanupFunction(obj, previousDetectionReport));
            gridScan = obj.gridPointsMesh;
            dimensions = find(size(gridScan.x) ~= 1);
            width = zeros(3, 1);
            width(dimensions) = obj.gridResolution(dimensions);
            obj.detectionSimulationReport = struct( ...
                'globalPDrealized', [], ...
                'globalPDmodeled', [], ...
                'globalPFArealized', [], ...
                'globalPFAmodel', [], ...
                'meanSNRsRealized', [], ...
                'meanSNRsModel', [], ...
                'numberOfUtilizedNodes', []);
            obj.detectionSimulationReport.globalPDrealized = zeros(prod(obj.gridSize), 1);
            obj.detectionSimulationReport.meanSNRsModel = 10*log10(sum(mean(obj.averageSNR_lin, 4), [1 2]));
            obj.detectionSimulationReport.globalPDmodeled = config.detector.globalPD;
            obj.detectionSimulationReport.meanSNRsRealized = zeros(prod(obj.gridSize), 1);
            obj.detectionSimulationReport.numberOfUtilizedNodes = zeros(prod(obj.gridSize), 1);
            obj.detectionSimulationReport.localDetectionNodes = zeros(obj.network.numberOfActiveBistaticPairs, 1);
            for mcID = 1 : mc.numberOfTrials
                if options.randomOnCell
                    obj.interfaces.settargetpositions("width", width);
                end
                obj.applyspatialprocessing;
                switch obj.network.networkMode
                    case 'multiStatic'
                        for mcpID = 1 : mc.numberOfTrialsParallel
                            targetCellIDs = obj.detectionReport(mcpID).cellIDs(1 : min(end, ceil(options.numberOfTargets)));
                            obj.detectionSimulationReport.globalPDrealized(targetCellIDs) = obj.detectionSimulationReport.globalPDrealized(targetCellIDs) + 1/numberOfTotalTrials;
                            obj.detectionSimulationReport.meanSNRsRealized(targetCellIDs) = obj.detectionSimulationReport.meanSNRsRealized(targetCellIDs) + obj.detectionReport(1, mcpID).detectedProcessedPower(1 : min(end, ceil(options.numberOfTargets)));
                            obj.detectionSimulationReport.numberOfUtilizedNodes(targetCellIDs) = obj.detectionSimulationReport.numberOfUtilizedNodes(targetCellIDs) + obj.detectionReport(1, mcpID).numberOfUtilizedNodes(1 : min(end, ceil(options.numberOfTargets)));
                        end
                        obj.detectionSimulationReport.localDetectionNodes = sum([obj.detectionReport.localDetectionNodes], 2);
                    case 'monoStatic'
                        for mcpID = 1 : mc.numberOfTrialsParallel
                            targetCellIDs = cell(1, obj.network.numberOfParallelSignalFusions);
                            SNRmonostatic = cell(1, obj.network.numberOfParallelSignalFusions);
                            for fusionID = 1 : obj.network.numberOfParallelSignalFusions
                                targetCellIDs{fusionID} = obj.detectionReport(fusionID, mcpID).cellIDs;
                                SNRmonostatic{fusionID} = obj.detectionReport(fusionID, mcpID).detectedProcessedPower(1 : end - 1);
                            end
                            targetCellIDs = cell2mat(targetCellIDs);
                            SNRmonostatic = cell2mat(SNRmonostatic);
                            [uniqueTargetCellIDs, uniqueIndices] = unique(targetCellIDs);
                            for detectionID = 1 : length(uniqueTargetCellIDs)
                                sameIndices = targetCellIDs == uniqueTargetCellIDs(detectionID);
                                SNRmonostatic(sameIndices) = max(SNRmonostatic(sameIndices));
                            end
                            targetCellIDs = targetCellIDs(uniqueIndices);
                            SNRmonostatic = SNRmonostatic(uniqueIndices);
                            obj.detectionSimulationReport.globalPDrealized(targetCellIDs) = obj.detectionSimulationReport.globalPDrealized(targetCellIDs) + 1/numberOfTotalTrials;
                            obj.detectionSimulationReport.meanSNRsRealized(targetCellIDs) = obj.detectionSimulationReport.meanSNRsRealized(targetCellIDs) + SNRmonostatic;
                        end
                end
                if ~mod(mcID, 100)
                    fprintf('trial = %d\n', mcID);
                end
            end
            obj.detectionSimulationReport.meanSNRsRealized = 10*log10(obj.detectionSimulationReport.meanSNRsRealized./obj.detectionSimulationReport.globalPDrealized/numberOfTotalTrials);
            obj.detectionSimulationReport.meanSNRsRealized(~logical(obj.detectionSimulationReport.globalPDrealized)) = -inf;
            obj.detectionSimulationReport.numberOfUtilizedNodes = obj.detectionSimulationReport.numberOfUtilizedNodes./obj.detectionSimulationReport.globalPDrealized/numberOfTotalTrials;
            function cleanupFunction(obj, previousDetectionReport)
                obj.detectionReport = previousDetectionReport;
                obj.interfaces.settargetpositions("width", 0);
            end
        end

        function fig = visualizedetectionsimulation(obj, options)
            arguments
                obj
                options.saveFigures (1, 1) {mustBeNumericOrLogical, mustBeMember(options.saveFigures, [0, 1])} = false
                options.monoStaticNetworkRXID (1, 1) {mustBePositive, mustBeInteger} = 1
            end
            if isempty(obj.detectionSimulationReport.globalPDrealized)
                fprintf('detection simulation had not executed\n');
                return;
            end
            gridScan = obj.gridPointsMesh;
            dimensions = {"x", "y", "z"};
            dims = size(gridScan.x) ~= 1;
            dimensions = dimensions(dims);
            xLabel = dimensions{1} + " (km)";
            yLabel = dimensions{2} + " (km)";
            x1 = obj.gridPoints{1}/1e3;
            x2 = obj.gridPoints{2}/1e3;
            visibleZone = reshape(~obj.blindZone, obj.gridSize([2 1 3]));

            globalPDmodeled = obj.detectionSimulationReport.globalPDmodeled(:, options.monoStaticNetworkRXID);
            meanSNRsModel = obj.detectionSimulationReport.meanSNRsModel(:, options.monoStaticNetworkRXID);
            globalPDRealized = reshape(obj.detectionSimulationReport.globalPDrealized(:, options.monoStaticNetworkRXID), obj.gridSize([2 1 3]));
            meanSNRsRealized = reshape(obj.detectionSimulationReport.meanSNRsRealized(:, options.monoStaticNetworkRXID), obj.gridSize([2 1 3]));
            numberOfUtilizedNodes = reshape(obj.detectionSimulationReport.numberOfUtilizedNodes(:, options.monoStaticNetworkRXID), obj.gridSize([2 1 3]));

            if ~isempty(obj.interfaces.targets)
                x = obj.interfaces.targets.position(dims, :)/1e3;
            end
            posRX = [obj.network.activeReceivingNodes.position]/1e3;
            posTX = [obj.network.activeTransmittingNodes.position]/1e3;
            posRX = posRX(dims, :); posTX = posTX(dims, :);

            fig(1) = figure;
            img = imagesc(x1, x2, globalPDRealized);
            ylim([-12 20]); xlim([-20 30]);
            colorbar; colormap('gray'); clim([0 1]);
            ax = gca; set(ax, 'Ydir', 'Normal');
            % set(img, 'AlphaData', visibleZone);
            delete(datatip(img, 2, 2));
            grid on; grid minor;
            xlabel(xLabel); ylabel(yLabel); zlabel('p_D');
            title(sprintf('Probability of detection (total P_d = %g)', sum(globalPDRealized, 'all')), sprintf('modeled P_d = %g', globalPDmodeled)); hold off;
            img.DataTipTemplate.DataTipRows(1).Label = "x";
            img.DataTipTemplate.DataTipRows(1).Value = gridScan.x;
            img.DataTipTemplate.DataTipRows(2).Label = "y";
            img.DataTipTemplate.DataTipRows(2).Value = gridScan.y;
            img.DataTipTemplate.DataTipRows(3).Label = "globalPD";
            img.DataTipTemplate.DataTipRows(3).Value = globalPDRealized;
            hold on;
            % if ~isempty(obj.interfaces.targets)
            %     plot(x(1, :), x(2, :), '+g', 'LineWidth', 2);
            % end
            plot(posRX(1, :), posRX(2, :), 'xb', 'LineWidth', 2);
            plot(posTX(1, :), posTX(2, :), '+r', 'LineWidth', 2);
            hold off; drawnow;

            fig(2) = figure;
            img = imagesc(x1, x2, meanSNRsRealized);
            colorbar; colormap('default');
            ax = gca; set(ax, 'Ydir', 'Normal');
            set(img, 'AlphaData', visibleZone & ~isinf(meanSNRsRealized));
            delete(datatip(img, 2, 2));
            grid off; grid on; grid minor;
            xlabel(xLabel); ylabel(yLabel);
            title('realized SNR with straddle loss averaged over trials', sprintf('modeled SNR = %g', meanSNRsModel)); hold off;
            img.DataTipTemplate.DataTipRows(1).Label = "x";
            img.DataTipTemplate.DataTipRows(1).Value = gridScan.x;
            img.DataTipTemplate.DataTipRows(2).Label = "y";
            img.DataTipTemplate.DataTipRows(2).Value = gridScan.y;
            img.DataTipTemplate.DataTipRows(3).Label = "mean SNR";
            img.DataTipTemplate.DataTipRows(3).Value = meanSNRsRealized;
            hold on;
            if ~isempty(obj.interfaces.targets)
                plot(x(1, :), x(2, :), '+k', 'LineWidth', 2);
            end
            plot(posRX(1, :), posRX(2, :), 'xb', 'LineWidth', 2);
            plot(posTX(1, :), posTX(2, :), '+r', 'LineWidth', 2);
            hold off; drawnow;

            fig(3) = figure;
            img = imagesc(x1, x2, numberOfUtilizedNodes);
            ylim([-12 20]); xlim([-20 30]);
            colorbar; colormap('default');
            ax = gca; set(ax, 'Ydir', 'Normal');
            set(img, 'AlphaData', visibleZone & ~isnan(numberOfUtilizedNodes));
            delete(datatip(img, 2, 2));
            grid off; grid on; grid minor;
            xlabel(xLabel); ylabel(yLabel);
            title('number of utilized nodes averaged over trials'); hold off;
            img.DataTipTemplate.DataTipRows(1).Label = "x";
            img.DataTipTemplate.DataTipRows(1).Value = gridScan.x;
            img.DataTipTemplate.DataTipRows(2).Label = "y";
            img.DataTipTemplate.DataTipRows(2).Value = gridScan.y;
            img.DataTipTemplate.DataTipRows(3).Label = "mean #nodes";
            img.DataTipTemplate.DataTipRows(3).Value = numberOfUtilizedNodes;
            hold on;
            if ~isempty(obj.interfaces.targets)
                plot(x(1, :), x(2, :), '+k', 'LineWidth', 2);
            end
            plot(posRX(1, :), posRX(2, :), 'xb', 'LineWidth', 2);
            plot(posTX(1, :), posTX(2, :), '+r', 'LineWidth', 2);
            text(posRX(1, :), posRX(2, :), num2str(obj.detectionSimulationReport.localDetectionNodes), ...
                "FontWeight", "normal", ...
                "VerticalAlignment", "bottom", ...
                "HorizontalAlignment", "left");
            hold off; drawnow;
        end

        function simulatecoverage(obj, options)
            arguments
                obj
                options.meanRCS_dbsm (1, 1) double = 0
                options.onCellCenters (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.onCellCenters, [0, 1])} = 0;
                options.neighbourOffset (1, 1) double {mustBeNonnegative} = 1
            end
            mc = obj.monteCarlo;
            config = obj.configuration;
            numberOfTotalTrials = mc.numberOfTrials*mc.numberOfTrialsParallel;
            originalTargets = obj.interfaces.targets;
            previousDetectionReport = obj.detectionReport;
            cleanup = onCleanup(@() cleanupFunction(obj, originalTargets, previousDetectionReport));
            gridScan = obj.gridPointsMesh;
            dimensions = find(size(gridScan.x) ~= 1);
            visibleZone = ~obj.blindZone;
            cellPositions = cat(3, gridScan.x, gridScan.y, gridScan.z);
            cellPositions = reshape(cellPositions(:, :, dimensions), [], numel(dimensions)).';
            cellPositions = cellPositions(:, visibleZone);
            width = zeros(3, 1);
            width(dimensions) = obj.gridResolution(dimensions);
            obj.coverageSimulationReport = struct( ...
                'targetCellIDs', [], ...
                'globalPDrealized', [], ...
                'globalPDmodel', [], ...
                'globalPFArealized', [], ...
                'globalPFAmodel', [], ...
                'meanSNRsModel', [], ...
                'meanSNRsRealized', []);
            obj.coverageSimulationReport.globalPDrealized = zeros(prod(obj.gridSize), 1);
            obj.coverageSimulationReport.globalPDmodel = zeros(prod(obj.gridSize), 1);
            obj.coverageSimulationReport.meanSNRsModel = -inf(prod(obj.gridSize), 1);
            obj.coverageSimulationReport.meanSNRsRealized = zeros(prod(obj.gridSize), 1);
            visibleCellIDs = 1 : size(cellPositions, 2); %% ayar
            numberOfCells = length(visibleCellIDs);
            targetCellIDs = find(visibleZone);
            targetCellIDs = targetCellIDs(visibleCellIDs);
            obj.coverageSimulationReport.targetCellIDs = targetCellIDs;
            targetPositions = zeros(3, numberOfCells);
            targetPositions(dimensions, :) = cellPositions(dimensions, visibleCellIDs);
            currentDetector = config.detector;
            for targetID = 1 : numberOfCells
                targetCellID = targetCellIDs(targetID);
                targetNeighbourCellIDs = obj.neighbourssquarewindow(targetCellID, "offset", options.neighbourOffset);
                obj.interfaces.settargets(target( ...
                    "position", targetPositions(:, targetID), ...
                    "meanRCS_dbsm", options.meanRCS_dbsm));
                obj.coverageSimulationReport.meanSNRsModel(targetCellID) = 10*log10(sum(mean(obj.averageSNR_lin, 4), [1 2]));
                currentDetector.setSNR(10*log10(mean(obj.averageSNR_lin, 4)).');
                obj.coverageSimulationReport.globalPDmodel(targetCellID) = currentDetector.globalPD;
                for mcID = 1 : mc.numberOfTrials
                    if ~options.onCellCenters
                        obj.interfaces.settargetpositions("width", width);
                    end
                    obj.applyspatialprocessing;
                    switch obj.network.networkMode
                        case 'multiStatic'
                            for mcpID = 1 : mc.numberOfTrialsParallel
                                targetIndex = ismember(obj.detectionReport(mcpID).cellIDs, targetNeighbourCellIDs);
                                if any(targetIndex)
                                    obj.coverageSimulationReport.globalPDrealized(targetCellID) = obj.coverageSimulationReport.globalPDrealized(targetCellID) + 1/numberOfTotalTrials;
                                    obj.coverageSimulationReport.meanSNRsRealized(targetCellID) = obj.coverageSimulationReport.meanSNRsRealized(targetCellID) + mean(obj.detectionReport(1, mcpID).detectedProcessedPower(find(targetIndex)));
                                end
                            end
                        case 'monoStatic'
                            for mcpID = 1 : mc.numberOfTrialsParallel
                                SNRmax = -inf;
                                for fusionID = 1 : obj.network.numberOfParallelSignalFusions
                                    targetIndex = ismember(obj.detectionReport(fusionID, mcpID).cellIDs, targetNeighbourCellIDs);
                                    if any(targetIndex)
                                        SNRmonostatic = mean(obj.detectionReport(fusionID, mcpID).detectedProcessedPower(find(targetIndex)));
                                        if SNRmonostatic > SNRmax
                                            SNRmax = SNRmonostatic;
                                        end
                                    end
                                end
                                if ~isinf(SNRmax)
                                    obj.coverageSimulationReport.globalPDrealized(targetCellID) = obj.coverageSimulationReport.globalPDrealized(targetCellID) + 1/numberOfTotalTrials;
                                    obj.coverageSimulationReport.meanSNRsRealized(targetCellID) = obj.coverageSimulationReport.meanSNRsRealized(targetCellID) + SNRmax;
                                end
                            end
                    end
                    if ~mod(mcID, 100)
                        fprintf('trial = %d/%d\n', mcID, mc.numberOfTrials);
                    end
                end
                obj.coverageSimulationReport.meanSNRsRealized(targetCellID) = 10*log10(obj.coverageSimulationReport.meanSNRsRealized(targetCellID)./obj.coverageSimulationReport.globalPDrealized(targetCellID)/numberOfTotalTrials);
                if ~mod(targetID, 100)
                    fprintf('target = %d/%d\n', targetID, numberOfCells);
                end
            end
            obj.coverageSimulationReport.meanSNRsRealized(~logical(obj.coverageSimulationReport.globalPDrealized)) = -inf;
            function cleanupFunction(obj, originalTargets, previousDetectionReport)
                obj.detectionReport = previousDetectionReport;
                obj.interfaces.targets = originalTargets;
                obj.interfaces.settargetpositions("width", 0);
            end
        end

        function analyticalcoverage(obj, options)
            arguments
                obj
                options.meanRCS_dbsm (1, 1) double = 0
            end
            originalTargets = obj.interfaces.targets;
            cleanup = onCleanup(@() cleanupFunction(obj, originalTargets));
            gridScan = obj.gridPointsMesh;
            visibleZone = ~obj.blindZone;
            obj.coverageSimulationReport = struct( ...
                'targetCellIDs', [], ...
                'globalPDrealized', [], ...
                'globalPDmodel', [], ...
                'globalPFArealized', [], ...
                'globalPFAmodel', [], ...
                'meanSNRsModel', [], ...
                'meanSNRsRealized', []);
            cellPositions = reshape(permute(cat(4, gridScan.x, gridScan.y, gridScan.z), [4 1 2 3]), [3 prod(obj.gridSize)]);
            cellPositions = cellPositions(:, visibleZone);
            obj.coverageSimulationReport.globalPDmodel = zeros(prod(obj.gridSize), 1);
            obj.coverageSimulationReport.meanSNRsModel = -inf(prod(obj.gridSize), 1);
            obj.coverageSimulationReport.targetCellIDs = find(visibleZone);
            obj.interfaces.settargets(target( ...
                "position", cellPositions, ...
                "meanRCS_dbsm", options.meanRCS_dbsm));
            obj.coverageSimulationReport.meanSNRsModel(obj.coverageSimulationReport.targetCellIDs) = 10*log10(sum(mean(obj.averageSNR_lin, 4), [1 2]));
            obj.coverageSimulationReport.globalPDmodel(obj.coverageSimulationReport.targetCellIDs) = obj.configuration.detector.globalPD;
            function cleanupFunction(obj, originalTargets)
                obj.interfaces.targets = originalTargets;
                obj.interfaces.settargetpositions("width", 0);
            end
        end

        function visualizecoveragesimulation(obj, options)
            arguments
                obj
                options.saveFigures (1, 1) {mustBeNumericOrLogical, mustBeMember(options.saveFigures, [0, 1])} = false
                options.monoStaticNetworkRXID (1, 1) {mustBePositive, mustBeInteger} = 1
                options.contourLevelDetection (1, 1) {mustBeNonnegative, mustBeInRange(options.contourLevelDetection, 0, 1)} = 0.85
                options.saveFigure (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.saveFigure, [0, 1])} = false
                options.header (1, :) char = datetime
            end
            if isempty(obj.coverageSimulationReport.targetCellIDs)
                fprintf('coverage simulation had not executed\n');
                return;
            end
            gridScan = obj.gridPointsMesh;
            dimensions = {"x", "y", "z"};
            dimensions = dimensions(size(gridScan.x) ~= 1);
            xLabel = dimensions{1} + " (km)";
            yLabel = dimensions{2} + " (km)";
            x1 = obj.gridPoints{1}/1e3;
            x2 = obj.gridPoints{2}/1e3;
            visibleZone = reshape(~obj.blindZone, obj.gridSize([2 1 3]));

            globalPDmodel = reshape(obj.coverageSimulationReport.globalPDmodel(:, options.monoStaticNetworkRXID), obj.gridSize([2 1 3]));
            meanSNRsModel = reshape(obj.coverageSimulationReport.meanSNRsModel(:, options.monoStaticNetworkRXID), obj.gridSize([2 1 3]));
            if ~isempty(obj.coverageSimulationReport.globalPDrealized)
                globalPDrealized = reshape(obj.coverageSimulationReport.globalPDrealized(:, options.monoStaticNetworkRXID), obj.gridSize([2 1 3]));
                meanSNRsRealized = reshape(obj.coverageSimulationReport.meanSNRsRealized(:, options.monoStaticNetworkRXID), obj.gridSize([2 1 3]));
            else
                globalPDrealized = globalPDmodel;
                meanSNRsRealized = meanSNRsModel;
            end

            fig1 = figure;
            img = imagesc(x1, x2, globalPDmodel); ylim([-12 20]); xlim([-20 30]);
            colorbar; colormap('gray'); clim([0 1]);
            ax = gca; set(ax, 'Ydir', 'Normal');
            % set(img, 'AlphaData', visibleZone);
            delete(datatip(img, 2, 2));
            % grid on; grid minor;
            grid off;
            xlabel(xLabel); ylabel(yLabel); zlabel('p_D');
            hold off;
            img.DataTipTemplate.DataTipRows(1).Label = "x";
            img.DataTipTemplate.DataTipRows(1).Value = gridScan.x;
            img.DataTipTemplate.DataTipRows(2).Label = "y";
            img.DataTipTemplate.DataTipRows(2).Value = gridScan.y;
            img.DataTipTemplate.DataTipRows(3).Label = "globalPD";
            img.DataTipTemplate.DataTipRows(3).Value = globalPDmodel;
            hold off;
            if options.saveFigure
                figureName = [char(options.header) '_ideal_PD'];
                savefig(fig1, ['C:\GitRepo\MSRS\figuresSim\' figureName '.fig']);
                saveas(fig1, ['C:\GitRepo\MSRS\figuresSim\' figureName '.eps'], 'epsc');
            else
                title('Modeled Probability of Detection'); drawnow;
            end

            fig2 = figure;
            img = imagesc(x1, x2, globalPDrealized); ylim([-12 20]); xlim([-20 30]);
            colorbar; colormap('gray'); clim([0 1]);
            ax = gca; set(ax, 'Ydir', 'Normal');
            % set(img, 'AlphaData', visibleZone);
            delete(datatip(img, 2, 2));
            grid on; grid minor;
            xlabel(xLabel); ylabel(yLabel); zlabel('p_D');
            hold off;
            img.DataTipTemplate.DataTipRows(1).Label = "x";
            img.DataTipTemplate.DataTipRows(1).Value = gridScan.x;
            img.DataTipTemplate.DataTipRows(2).Label = "y";
            img.DataTipTemplate.DataTipRows(2).Value = gridScan.y;
            img.DataTipTemplate.DataTipRows(3).Label = "globalPD";
            img.DataTipTemplate.DataTipRows(3).Value = globalPDrealized;
            hold off;
            if options.saveFigure
                figureName = [char(options.header) '_realized_PD'];
                savefig(fig2, ['C:\GitRepo\MSRS\figuresSim\' figureName '.fig']);
                saveas(fig2, ['C:\GitRepo\MSRS\figuresSim\' figureName '.eps'], 'epsc');
            else
                title('Realized Probability of Detection'); drawnow;
            end

            upperLim = max(max(meanSNRsModel(meanSNRsModel < 300), [], 'all'), max(meanSNRsRealized(meanSNRsRealized < 300), [], 'all'));
            lowerLim = min(min(meanSNRsModel(~isinf(meanSNRsModel)), [], 'all'), min(meanSNRsRealized(~isinf(meanSNRsRealized)), [], 'all'));

            fig3 = figure;
            img = imagesc(x1, x2, meanSNRsModel); ylim([-12 20]); xlim([-20 30]); zlim([-20 inf]);
            colorbar; colormap('default'); clim([lowerLim, upperLim]); clim([-20, upperLim]);
            ax = gca; set(ax, 'Ydir', 'Normal');
            set(img, 'AlphaData', visibleZone);
            delete(datatip(img, 2, 2));
            % grid on; grid minor;
            grid off;
            xlabel(xLabel); ylabel(yLabel);
            hold off;
            img.DataTipTemplate.DataTipRows(1).Label = "x";
            img.DataTipTemplate.DataTipRows(1).Value = gridScan.x;
            img.DataTipTemplate.DataTipRows(2).Label = "y";
            img.DataTipTemplate.DataTipRows(2).Value = gridScan.y;
            img.DataTipTemplate.DataTipRows(3).Label = "mean SNR";
            img.DataTipTemplate.DataTipRows(3).Value = meanSNRsModel;
            hold off;
            if options.saveFigure
                figureName = [char(options.header) '_ideal_SNR'];
                savefig(fig3, ['C:\GitRepo\MSRS\figuresSim\' figureName '.fig']);
                saveas(fig3, ['C:\GitRepo\MSRS\figuresSim\' figureName '.eps'], 'epsc');
            else
                title('Modeled SNR averaged over trials'); drawnow;
            end

            fig4 = figure;
            img = imagesc(x1, x2, meanSNRsRealized); ylim([-12 20]); xlim([-20 30]);
            colorbar; colormap('default'); clim([lowerLim, upperLim]);
            ax = gca; set(ax, 'Ydir', 'Normal');
            set(img, 'AlphaData', visibleZone);
            delete(datatip(img, 2, 2));
            grid off; grid on; grid minor;
            xlabel(xLabel); ylabel(yLabel);
            hold off;
            img.DataTipTemplate.DataTipRows(1).Label = "x";
            img.DataTipTemplate.DataTipRows(1).Value = gridScan.x;
            img.DataTipTemplate.DataTipRows(2).Label = "y";
            img.DataTipTemplate.DataTipRows(2).Value = gridScan.y;
            img.DataTipTemplate.DataTipRows(3).Label = "mean SNR";
            img.DataTipTemplate.DataTipRows(3).Value = meanSNRsRealized;
            hold off;
            if options.saveFigure
                figureName = [char(options.header) '_realized_SNR'];
                savefig(fig4, ['C:\GitRepo\MSRS\figuresSim\' figureName '.fig']);
                saveas(fig4, ['C:\GitRepo\MSRS\figuresSim\' figureName '.eps'], 'epsc');
            else
                title('Realized SNR with straddle loss averaged over trials'); drawnow;
            end

            fig5 = figure;
            switch obj.network.networkMode
                case 'monoStatic'
                    edgeColor = 'k';
                case 'multiStatic'
                    edgeColor = 'r';
            end
            contour(x1, x2, globalPDrealized, [-1 options.contourLevelDetection], 'LineWidth', 2, 'ShowText', 'on', 'EdgeColor', edgeColor); ylim([-12 20]); xlim([-20 30]);
            grid on; grid minor;
            xlabel(xLabel); ylabel(yLabel); zlabel('p_D');
            hold off;
            if options.saveFigure
                figureName = [char(options.header) '_contour_PD'];
                savefig(fig5, ['C:\GitRepo\MSRS\figuresSim\' figureName '.fig']);
                saveas(fig5, ['C:\GitRepo\MSRS\figuresSim\' figureName '.eps'], 'epsc');
            else
                title('Probability of detection'); drawnow;
            end
        end
        
        %%% visualization methods

        function visualizefilteredsignals(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.receivingNodeIDs (1, :) double {mustBeInteger, mustBeNonnegative} = 1 : obj.network.numberOfActiveReceivingNodes
                options.trialID (1, 1) double {mustBePositive, mustBeInteger} = 1
            end
            mustBeInRange(options.receivingNodeIDs, 1, obj.network.numberOfReceivingNodes);
            if isempty(obj.signalsMatchFiltered)
                fprintf('"signalsMatchFiltered" is empty\n');
                return;
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            for rxID = options.receivingNodeIDs
                s = obj.signalsMatchFiltered{rxID}(:, :, :, options.trialID);
                Ts = obj.network.activeReceivingNodes(rxID).samplingPeriod;
                switch obj.network.networkMode
                    case "monoStatic"
                        txIDs = obj.network.monoStaticTransmitterIDs(rxID);
                    case "multiStatic"
                        txIDs = 1 : obj.network.numberOfActiveTransmittingNodes;
                end
                for txID = txIDs
                    switch obj.network.activeTransmittingNodes(txID).transmissionType
                        case "continuous"
                            error('not implemented');
                        case "pulsed"
                            L = obj.network.pulseWidthSample(rxID, txID) - 1;
                            N = obj.network.activeReceivingNodes(rxID).numberOfSamplesPerCPI;
                    end
                    t = (-L + 1 : N)*Ts*1e6;
                    if ~isscalar(options.receivingNodeIDs)
                        plot(t, 20*log10(abs(s(:, ceil(end/2), txID))));
                    else
                        plot(t, 20*log10(abs(s(:, :, txID))));
                    end
                    hold on;
                    T = obj.thresholdCFAR{rxID}(:, :, options.trialID);
                    if ~isscalar(options.receivingNodeIDs)
                        plot(t, 10*log10(T(:, ceil(end/2), txID)), 'LineStyle', '--', 'LineWidth', 2);
                    else
                        plot(t, 10*log10(T), 'LineStyle', '--', 'LineWidth', 2);
                    end
                end
            end
            grid off; grid on; grid minor;
            xlabel('time (\mus)'); ylabel('power (dB)');
            title('filtered signal');
            if ~isscalar(options.receivingNodeIDs)
                leg = legend(num2str(options.receivingNodeIDs.'), 'Location', 'best');
                title(leg, 'RX ID');
            end
            drawnow; hold off;
        end

        function visualizehypothesizedvariables(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.transmittingNodeID (1, 1) double {mustBeInteger, mustBePositive} = 1
                options.receivingNodeID (1, 1) double {mustBeInteger, mustBePositive} = 1
                options.variable (1, 1) string {mustBeMember(options.variable, ["time", "averageSNR", "coveredZone", "indices", "weights", "phase"])} = "time"
                options.dimension (1, 1) string {mustBeMember(options.dimension, ["x-y", "y-z", "z-x", "x-y-z"])} = "x-y"
            end
            mustBeInRange(options.transmittingNodeID, 1, obj.network.numberOfActiveTransmittingNodes);
            mustBeInRange(options.receivingNodeID, 1, obj.network.numberOfActiveReceivingNodes);
            if isempty(obj.dictionaryCellMapping)
                fprintf('hypothesized variables are not set\n');
                return;
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            gridScan = obj.gridPointsMesh;
            switch options.variable
                case "time"
                    var = nan(obj.gridSize);
                    var(~obj.blindZone) = squeeze(obj.hypothesizedVisibleInterface.timeDelay(options.transmittingNodeID, options.receivingNodeID, :))*1e6;
                    labelVariable = 'time delay (\mus)';
                    titleVariable = 'hypothesized time delays';
                case "averageSNR"
                    var = nan(obj.gridSize);
                    var(~obj.blindZone) = obj.hypothesizedVisibleInterface.averageSNR_dB(options.transmittingNodeID, options.receivingNodeID, :);
                    labelVariable = 'amplitude (dB)';
                    titleVariable = 'hypothesized average SNRs for 1 m^2 target';
                case "coveredZone"
                    var = ~obj.blindZone;
                    labelVariable = 'isVisible';
                    titleVariable = 'visibility of cells';
                case "indices"
                    var = nan(prod(obj.gridSize), 1);
                    switch obj.network.networkMode
                        case 'multiStatic'
                            var(~obj.blindZone) = obj.dictionaryCellMapping(:, options.transmittingNodeID, options.receivingNodeID);
                        case 'monoStatic'
                            var(~obj.blindZone) = obj.dictionaryCellMapping{options.transmittingNodeID}(:, options.receivingNodeID);
                    end
                    labelVariable = 'index';
                    titleVariable = 'integration indices';
                case "weights"
                    var = nan(prod(obj.gridSize), 1);
                    var(~obj.blindZone) = 20*log10(abs(squeeze(obj.dictionaryComplexWeights(options.transmittingNodeID, options.receivingNodeID, :))));
                    labelVariable = 'amplitude (dB)';
                    titleVariable = 'amplitude of integration weights';
                case "phase"
                    %%% not implemented
                    % var = 180*squeeze(obj.hypothesizedPhases(options.transmittingNodeID, options.receivingNodeID, :))/pi;
                    labelVariable = 'phase (°)';
                    titleVariable = 'hypothesized phases';
            end
            switch options.dimension
                case "x-y"
                    var = reshape(var, [obj.gridSize]);
                    x = gridScan.x(:, :, 1); x = x(:)/1e3;
                    y = gridScan.y(:, :, 1); y = y(:)/1e3;
                    varPlot = var(:, :, 1);
                    if strcmp(options.variable, "coveredZone")
                        visibleCells = varPlot == 1;
                        plot3(x(visibleCells), y(visibleCells), varPlot(visibleCells), '.b'); hold on;
                        plot3(x(~visibleCells), y(~visibleCells), varPlot(~visibleCells), '.'); view(0, 90);
                    else
                        plot3(x, y, varPlot(:), '.'); hold on;
                    end
                    if size(var, 3) ~= 1
                        for zID = unique(ceil(linspace(2, size(var, 3), 7)))
                            varPlot = var(:, :, zID);
                            plot3(x, y, varPlot(:), '.');
                        end
                    end
                    grid off; grid on; grid minor;
                    xlabel('x (km)'); ylabel('y (km');
                    zlabel(labelVariable);
                    title(titleVariable);
                case "y-z"
                    var = reshape(var, [obj.gridSize]).';
                    y = gridScan.y(:, 1, :); y = y(:)/1e3;
                    z = gridScan.z(:, 1, :); z = z(:)/1e3;
                    varPlot = var(:, 1, :);
                    plot3(y, z, varPlot(:), '.'); hold on;
                    if size(var, 2) ~= 1
                        for xID = unique(ceil(linspace(2, size(var, 2), 7)))
                            varPlot = var(:, xID, :);
                            plot3(y, z, varPlot(:), '.');
                        end
                    end
                    grid off; grid on; grid minor;
                    xlabel('y (km)'); ylabel('z (km)');
                    zlabel(labelVariable);
                    title(titleVariable);
                case "z-x"
                    var = reshape(var, [obj.gridSize]).';
                    z = gridScan.z(1, :, :); z = z(:)/1e3;
                    x = gridScan.x(1, :, :); x = x(:)/1e3;
                    varPlot = var(1, :, :);
                    plot3(z, x, varPlot(:), '.'); hold on;
                    if size(var, 1) ~= 1
                        for yID = unique(ceil(linspace(2, size(var, 1), 7)))
                            varPlot = var(yID, :, :);
                            plot3(z, x, varPlot(:), '.');
                        end
                    end
                    grid off; grid on; grid minor;
                    xlabel('y (km)'); ylabel('z (km)');
                    zlabel(labelVariable);
                    title(titleVariable);
                case "x-y-z"
                    x = gridScan.x(:)/1e3;
                    y = gridScan.y(:)/1e3;
                    z = gridScan.z(:)/1e3;
                    scatter3(x, y, z, 0.25, var);
                    hcb = colorbar;
                    title(hcb, labelVariable);
                    grid off; grid on; grid minor;
                    xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
                    title('hypothesized locations');
            end
        end

        function fig = visualizeintegratedsignals(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.plotMode (1, 1) string {mustBeMember(options.plotMode, ["image", "mesh"])} = "image"
                options.trialID (1, 1) double {mustBePositive, mustBeInteger} = 1
                options.monoStaticNetworkRXID (1, 1) {mustBePositive, mustBeInteger} = 1
            end
            switch obj.network.networkMode
                case "multiStatic"
                    options.monoStaticNetworkRXID = 1;
            end
            config = obj.configuration;
            z = nan(prod(obj.gridSize), 1);
            z(~obj.blindZone) = 10*log10(obj.detectionReport(options.monoStaticNetworkRXID, options.trialID).integratedPower);
            switch obj.network.networkMode
                case "multiStatic"
                    pos = obj.detectionReport(options.trialID).estimatedPositions;
                case "monoStatic"
                    pos = obj.detectionReport(options.monoStaticNetworkRXID, options.trialID).estimatedPositions;
            end
            if isempty(options.figureID)
                fig = figure;
            else
                fig = figure(options.figureID);
            end
            gridScan = obj.gridPointsMesh;
            dims = find(size(gridScan.x) ~= 1);
            dimensions = {"x", "y", "z"};
            dimensions = dimensions(dims);
            switch length(dims)
                case 1
                    xLabel = dimensions{1} + " (km)";
                    plot(gridScan.(dimensions{1})/1e3, z);
                    grid off; grid on; grid minor;
                    xlabel(xLabel); ylabel('power (dB)');
                    title('Response of cells');
                case 2
                    Y = reshape(z, obj.gridSize([2 1 3]));
                    maxY = max(Y, [], 'all');
                    xLabel = dimensions{1} + " (km)";
                    yLabel = dimensions{2} + " (km)";
                    x1 = gridScan.(dimensions{1})/1e3;
                    x2 = gridScan.(dimensions{2})/1e3;
                    posRX = [obj.network.activeReceivingNodes.position]/1e3;
                    posTX = [obj.network.activeTransmittingNodes.position]/1e3;
                    posRX = posRX(dims, :);
                    posTX = posTX(dims, :);
                    switch options.plotMode
                        case "mesh"
                            hold on;
                            plot3(posRX(1, :), posRX(2, :), repmat(maxY, [1, size(posRX, 2)]), 'xm', 'LineWidth', 2);
                            plot3(posTX(1, :), posTX(2, :), repmat(maxY, [1, size(posTX, 2)]), '+r', 'LineWidth', 2);
                            estimations = pos(dims, :)/1e3;
                            plot3(estimations(1, :), estimations(2, :), repmat(maxY, [1, size(estimations, 2)]), 'oy', 'LineWidth', 2);
                            if ~isempty(obj.interfaces.targets)
                                x = obj.interfaces.targets.position(dims, :)/1e3;
                                plot3(x(1, :), x(2, :), repmat(maxY, [1, size(x, 2)]), '+k', 'LineWidth', 2);
                            end
                            plot3(x2(:), x1(:), 10*log10(config.detector.globalThreshold(options.monoStaticNetworkRXID)*ones(1, numel(x2))), 'k');
                            m = mesh(x2, x1, Y);
                            set(m, 'AlphaData', ~isinf(Y) & ~isnan(Y));
                            m.FaceColor = 'flat';
                        case "image"
                            hold on;
                            img = imagesc(obj.gridPoints{1}/1e3, obj.gridPoints{2}/1e3, Y);
                            % set(img, 'AlphaData', ~isinf(Y) & ~isnan(Y));
                            delete(datatip(img, 2, 2));
                            img.DataTipTemplate.DataTipRows(1).Label = "x";
                            img.DataTipTemplate.DataTipRows(1).Value = x1;
                            img.DataTipTemplate.DataTipRows(2).Label = "y";
                            img.DataTipTemplate.DataTipRows(2).Value = x2;
                            img.DataTipTemplate.DataTipRows(3).Label = "power";
                            img.DataTipTemplate.DataTipRows(3).Value = Y;
                            plot(posRX(1, :), posRX(2, :), 'xm', 'LineWidth', 1, 'MarkerSize', 10);
                            plot(posTX(1, :), posTX(2, :), '+r', 'LineWidth', 1, 'MarkerSize', 10);
                            estimations = pos(dims, :)/1e3;
                            plot(estimations(1, :), estimations(2, :), 'oy', 'LineWidth', 2);
                            if ~isempty(obj.interfaces.targets)
                                x = obj.interfaces.targets.position(dims, :)/1e3;
                                plot(x(1, :), x(2, :), '+k', 'LineWidth', 2);
                            end
                    end
                    colorbar; clim([0 inf]); colormap('winter'); view(0, 90);
                    xlim(obj.network.boundaryListened(1, :)/1e3); ylim(obj.network.boundaryListened(2, :)/1e3);
                    xlabel(xLabel); ylabel(yLabel); zlabel('power (dB)');
                    % switch obj.network.networkMode
                    %     case "multiStatic"
                    %         SNR = max(mean(obj.averageSNR_dB, [1 2]), [], [3 4]);
                    %         try
                    %             globalPD = max(config.detector.globalPD, [], 'all');
                    %         catch
                    %             globalPD = nan;
                    %         end
                    %     case "monoStatic"
                    %         SNR = max(obj.averageSNR_dB(:, options.monoStaticNetworkRXID, :, :), [], [3 4]);
                    %         globalPD = max(config.detector.globalPD(:, options.monoStaticNetworkRXID, :, :), [], [3 4]);
                    % end
                    % title(sprintf('average input SNR = %.3g dB, P_{FA}^{global} = %.3e, P_{D}^{global} = %.3g', SNR, config.globalPFA, globalPD)); hold off;
                case 3
            end
            % if isempty(pos)
            %     if ~isempty(obj.interfaces.targets)
            %         legend('RX', 'TX', 'targets', 'Location', 'best');
            %     else
            %         legend('RX', 'TX', 'Location', 'best');
            %     end
            % else
            %     if ~isempty(obj.interfaces.targets)
            %         legend('RX', 'TX', 'estimations', 'targets', 'Location', 'best');
            %     else
            %         legend('RX', 'TX', 'estimations', 'Location', 'best');
            %     end
            % end
            hold off; drawnow;
        end

        function visualizeestimation(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.trialID (1, 1) double {mustBePositive, mustBeInteger} = 1
                options.monoStaticNetworkRXID (1, 1) {mustBePositive, mustBeInteger} = 1
            end
            config = obj.configuration;
            switch obj.network.networkMode
                case "multiStatic"
                    pos = obj.detectionReport(options.trialID).estimatedPositions;
                case "monoStatic"
                    pos = obj.detectionReport(options.monoStaticNetworkRXID, options.trialID).estimatedPositions;
            end
            if isempty(pos)
                fprintf('estimations are empty\n');
                return;
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            gridScan = obj.gridPointsMesh;
            dims = size(gridScan.x) ~= 1;
            dimensions = {"x", "y", "z"};
            dimensions = dimensions(dims);
            xLabel = dimensions{1} + " (km)";
            yLabel = dimensions{2} + " (km)";
            [~, maxPowerDetection] = max(obj.detectionReport(options.monoStaticNetworkRXID, options.trialID).detectedProcessedPower);
            Y = 10*log10(obj.detectionReport(options.monoStaticNetworkRXID, options.trialID).detectedRawPower(:, maxPowerDetection));
            pos = pos(:, maxPowerDetection);
            posRX = [obj.network.activeReceivingNodes.position]/1e3;
            posTX = [obj.network.activeTransmittingNodes.position]/1e3;
            plot3(posRX(1, :), posRX(2, :), posRX(3, :), 'xb', 'LineWidth', 2, 'MarkerSize', 10);
            hold on; plot3(posTX(1, :), posTX(2, :), posTX(3, :), '+r', 'LineWidth', 2, 'MarkerSize', 10);
            if ~isempty(obj.interfaces.targets)
                x = obj.interfaces.targets.position(1, :)/1e3;
                y = obj.interfaces.targets.position(2, :)/1e3;
                z = obj.interfaces.targets.position(3, :)/1e3;
                plot3(x, y, z, '*k', 'LineWidth', 1, 'MarkerSize', 10);
            end
            xEstimation = pos(1)/1e3;
            yEstimation = pos(2)/1e3;
            zEstimation = pos(3)/1e3;
            plot3(xEstimation, yEstimation, zEstimation, 'oc', 'LineWidth', 2, 'MarkerSize', 10);
            text(posRX(1, :), posRX(2, :), posRX(3, :), num2str((1 : obj.network.numberOfActiveReceivingNodes).'), "FontSize", 14, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            posRX = repelem(posRX, 1, obj.interfaces.network.numberOfActiveTransmittingNodes);
            posTX = repmat(posTX, 1, obj.interfaces.network.numberOfActiveReceivingNodes);
            line([posRX(1, :); posTX(1, :)], [posRX(2, :); posTX(2, :)], [posRX(3, :); posTX(3, :)], 'lineStyle', '--', 'Color', 'k');
            text((posRX(1, :) + posTX(1, :))/2, (posRX(2, :) + posTX(2, :))/2, num2str(Y), ...
                'FontSize', 12, ...
                'Color', 'blue', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom');
            xlim(obj.network.boundaryListened(1, :)/1e3); ylim(obj.network.boundaryListened(2, :)/1e3);
            grid off; grid on; grid minor;
            xlabel(xLabel); ylabel(yLabel); zlabel('z (km)')
            switch obj.network.networkMode
                case "multiStatic"
                    SNR = max(mean(obj.averageSNR_dB, [1 2]), [], [3 4]);
                    try
                        globalPD = max(config.detector.globalPD, [], 'all');
                    catch
                        globalPD = nan;
                    end
                case "monoStatic"
                    SNR = max(obj.averageSNR_dB(:, options.monoStaticNetworkRXID, :, :), [], [3 4]);
                    globalPD = max(config.detector.globalPD(:, options.monoStaticNetworkRXID, :, :), [], [3 4]);
            end
            title(sprintf('average input SNR = %.3g dB, PFA_{global} = %.3e, PD_{global} = %.3g', SNR, config.globalPFA, globalPD)); hold off;
            legend('RX', 'TX', 'targets', 'estimations', 'Location', 'best');
            view(0, 90); hold off; drawnow;
        end
    end
end