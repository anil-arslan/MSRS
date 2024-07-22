classdef spu < handle
    %spu Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = private, GetAccess = public)
        interface (1, :) interface = interface.empty()
        monteCarlo (1, 1) struct = struct( ...
            'seed', 0, ...
            'numberOfTrials', 1, ...
            'numberOfTrialsParallel', 1)
        gridResolution (3, 1) double {mustBeNonnegative} = 100*ones(3, 1)
        configuration (1, 1) struct = struct('PFA', 1e-6)
        processingAlgorithm (1, 1) {mustBeInteger, mustBeInRange(processingAlgorithm, 1, 10)} = 1
        % 1: coherent detector (deterministic signal)
        % 2: coherent detector (stochastic signal, dependent phase, dependent amplitude)
        % 3: square-law envelope detector (non-fluctuating weak stochastic signal, independent phase, dependent amplitude)
        % 4: linear-law envelope detector (non-fluctuating strong stochastic signal, independent phase, dependent amplitude)
        % 5:
        % 6:
    end

    properties (Dependent)
        network (1, 1) radarNetwork
        gridPoints cell
        gridSize (1, 3) double {mustBeNonnegative}
        gridPointsMesh struct
        hypothesizedTimeDelays double
    end

    properties (SetAccess = private, GetAccess = public)
        integrationIndices double {mustBeInteger, mustBePositive} % (Ntx x Nrx x Ni vector)
    end

    properties (Dependent)
        inputSNR_dB double % dB power
        outputSNR_dB double % dB power
        outputSNR_lin double % linear scale power
        threshold function_handle
        ROC function_handle
        expectedAmplitudes double {mustBePositive} % linear scale magnitude normalized to first receiver (reference)
        noisePowersPerSample_W double {mustBePositive} % linear scale power
        integrationWeights double {mustBePositive}
        signalsMatchFilteredTrials double % (Ns + L - 1 x Nrx x Nmcp matrix) Ncmp : number of parallel trials
    end

    properties (SetAccess = private, GetAccess = public)
        signalsMatchFiltered double = 0 % (Ns + L - 1 x Nrx x Nmcp matrix) Ncmp : number of parallel trials
    end

    properties (Dependent)
        signalsIntegrated double % (Ni x Ncmp matrix)
        estimatedTimeDelaysFromMatchFiltration double % (Ntx x Nrx x Nmcp matrix) Ncmp : number of parallel trials
        timeOfArrivalErrorFromMatchFiltration double % (Ntx x Nrx x Nmcp matrix) Ncmp : number of parallel trials
        hypothesisTestingResults cell
        estimatedPositionsFromIntegration cell
        positionError double
        positionErrorTotal double
    end

    properties (Access = private)
        seedShuffle (1, 1) logical {mustBeNumericOrLogical, mustBeMember(seedShuffle, [0, 1])}
    end

    methods
        function obj = spu(options)
            arguments
                options.interface (1, 1) interface
                options.gridResolution (3, 1) double {mustBeNonnegative} = zeros(3, 1)
            end
            %spu Construct an instance of this class
            %   Detailed explanation goes here
            obj.interface = options.interface;
            if ~any(options.gridResolution)
                obj.gridResolution = obj.interface.network.defaultGridResolution;
            else
                obj.gridResolution = options.gridResolution;
            end
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
            net = obj.interface.network;
        end

        function SNRin = get.inputSNR_dB(obj)
            SNRin = obj.interface.inputSNR_dB;
        end

        function snr = get.outputSNR_lin(obj)
            snr = sum(10.^(.1*(obj.network.processingGain_dB.' + obj.inputSNR_dB)), 2);
        end

        function snr = get.outputSNR_dB(obj)
            snr = 10*log10(obj.outputSNR_lin);
        end

        function PD = get.ROC(obj)
            if obj.network.multiTransmitter
                PD = @(PFA, SNR) nan;
            else
                switch obj.network.networkCoherency
                    case "coherent"
                        switch obj.processingAlgorithm
                            case 1
                                PD = @(PFA, SNR) .5*erfc(erfcinv(2*PFA) - sqrt(SNR));
                            case 2
                            case 3
                            case 4
                                PD = @(PFA, SNR) marcumq(sqrt(2*SNR), sqrt(-2*log(PFA)));
                            otherwise
                                error('Processing algorithm %d cannot be used for the current configuration', obj.processingAlgorithm);
                        end
                    case "short-term coherent"
                        PD = @(PFA, SNR) nan;
                    case "incoherent"
                        PD = @(PFA, SNR) nan;
                end
            end
        end

        function T = get.threshold(obj)
            if obj.network.multiTransmitter
                T = @(PFA, SNR) nan;
            else
                switch obj.network.networkCoherency
                    case "coherent"
                        switch obj.processingAlgorithm
                            case 1
                                T = @(PFA, SNR) sqrt(SNR)*erfcinv(2*PFA);
                            case 2
                            case 3
                            case 4
                                T = @(PFA, SNR) sqrt(-SNR*log(PFA));
                            otherwise
                                error('Processing algorithm %d cannot be used for the current configuration', obj.processingAlgorithm);
                        end
                    case "short-term coherent"
                        T = @(PFA, SNR) nan;
                    case "incoherent"
                        T = @(PFA, SNR) nan;
                end
            end
        end

        function cfg = get.configuration(obj)
            cfg = obj.configuration;
            if obj.outputSNR_lin
                cfg.PD = obj.ROC(cfg.PFA, obj.outputSNR_lin);
                cfg.threshold = obj.threshold(cfg.PFA, obj.outputSNR_lin);
            else
                cfg.PD = obj.ROC(cfg.PFA, 1./obj.noisePowersPerSample_W);
                cfg.threshold = obj.threshold(cfg.PFA, 1./obj.noisePowersPerSample_W);
            end
            cfg.threshold_dB = 20*log10(abs(cfg.threshold));
        end

        function n = get.noisePowersPerSample_W(obj)
            n = 10.^(.1*[obj.network.activeReceivingNodes.noisePowerPerSample_dB]); % kTBN
        end

        function w = get.expectedAmplitudes(obj)
            w = 10.^(.05*obj.interface.receivedPowerFromScatterers_dBW);
        end

        function w = get.integrationWeights(obj)
            if obj.network.multiTransmitter
                w = nan;
                warning('not implemented');
            else
                switch obj.network.networkCoherency
                    case "coherent"
                        switch obj.processingAlgorithm
                            case 1 % Maximal ratio combining
                                w = obj.expectedAmplitudes./obj.noisePowersPerSample_W;
                                w(w == 0) = 1./obj.noisePowersPerSample_W(w == 0);
                                w(w == inf) = obj.expectedAmplitudes(w == inf);
                            case 4
                                w = obj.expectedAmplitudes./obj.noisePowersPerSample_W;
                                w(w == 0) = 1./obj.noisePowersPerSample_W(w == 0);
                                w(w == inf) = obj.expectedAmplitudes(w == inf);
                            otherwise
                                error('Processing algorithm %d cannot be used for the current configuration', obj.processingAlgorithm);
                        end
                    case "short-term coherent"
                        w = nan;
                        warning('not implemented');
                    case "incoherent"
                        w = nan;
                        warning('not implemented');
                end
            end
        end

        function posScan = get.gridPoints(obj)
            networkBoundary = obj.network.boundaryListened;
            posScan = cell(3, 1);
            for i = 1 : 3
                posScan{i} = networkBoundary(i, 1) : obj.gridResolution(i) : networkBoundary(i, 2);
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

        function tau = get.hypothesizedTimeDelays(obj)
            gridScan = obj.gridPointsMesh;
            hypothesizedPositions = permute(cat(4, gridScan.x, gridScan.y, gridScan.z), [4 1 2 3]);
            hypothesizedPositions = reshape(hypothesizedPositions, [3, 1, 1, prod(obj.gridSize)]);
            rangeReceivers = sqrt(sum((hypothesizedPositions - permute([obj.network.activeReceivingNodes.position], [1 3 2])).^2));
            rangeTransmitters = sqrt(sum((hypothesizedPositions - [obj.network.activeTransmittingNodes.position]).^2));
            tau = permute(rangeReceivers + rangeTransmitters, [2 3 4 1])/obj.network.speedOfLight;
        end

        function y = get.signalsMatchFilteredTrials(obj)
            % (Ns + L - 1 x Nrx x Ncmp matrix) Ncmp : number of parallel trials
            t = [obj.network.activeReceivingNodes.samplingInstants]; % Ns x Nrx matrix
            demodulator = exp(1j*2*pi*shiftdim([obj.network.activeTransmittingNodes.carrierFrequency], -1).*t);
            mf = obj.network.matchFilter; % L x Nrx x Ntx matrix
            y = zeros(size(t, 1) + size(mf, 1) - 1, obj.network.numberOfActiveReceivingNodes);
            for mcID = 1 : obj.monteCarlo.numberOfTrialsParallel
                s = obj.interface.signalSuperposed; % Ns x Nrx matrix
                for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                    y(:, rxID, mcID) = conv(s(: ,rxID).*demodulator(:, rxID), conj((mf(:, rxID, 1))));
                end
            end
        end

        function Y = get.signalsIntegrated(obj)
            % (Ni x Nmcp matrix)
            if obj.network.multiTransmitter
                Y = nan;
                warning('not implemented');
            else
                if isscalar(obj.signalsMatchFiltered) || isempty(obj.signalsMatchFiltered)
                    obj.setmatchfilteredsignals;
                end
                switch obj.network.networkCoherency
                    case "coherent"
                        switch obj.processingAlgorithm
                            case 1
                                yWeighted = obj.integrationWeights.*obj.signalsMatchFiltered;
                                Y = zeros(prod(obj.gridSize), obj.monteCarlo.numberOfTrialsParallel);
                                for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                                    Y = Y + yWeighted(obj.integrationIndices(1, rxID, :), rxID, :);
                                end
                            case 4 % linear-law envelope detector for strong signals
                                yWeighted = obj.integrationWeights.*abs(obj.signalsMatchFiltered);
                                Y = zeros(prod(obj.gridSize), obj.monteCarlo.numberOfTrialsParallel);
                                for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                                    Y = Y + yWeighted(obj.integrationIndices(1, rxID, :), rxID, :);
                                end
                            otherwise
                                error('Processing algorithm %d cannot be used for the current configuration', obj.processingAlgorithm);
                        end
                    case "short-term coherent"
                        Y = nan;
                        warning('not implemented');
                    case "incoherent"
                        Y = nan;
                        warning('not implemented');
                end
            end
        end

        function tau = get.estimatedTimeDelaysFromMatchFiltration(obj)
            if obj.network.multiTransmitter
                tau = nan;
            else
                if isscalar(obj.signalsMatchFiltered) || isempty(obj.signalsMatchFiltered)
                    obj.setmatchfilteredsignals;
                end
                [~, ind] = max(abs(obj.signalsMatchFiltered));
                tau = zeros(obj.network.numberOfActiveTransmittingNodes, obj.network.numberOfActiveReceivingNodes, obj.monteCarlo.numberOfTrialsParallel);
                for mcID = 1 : obj.monteCarlo.numberOfTrialsParallel
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        Ts = obj.network.activeReceivingNodes(rxID).samplingPeriod;
                        for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                            switch obj.network.activeTransmittingNodes(txID).transmissionType
                                case "continuous"
                                    error('not implemented');
                                case "pulsed"
                                    L = obj.network.pulseWidthSample(rxID, txID) - 2;
                                    N = obj.network.activeReceivingNodes(rxID).numberOfSamplesPerCPI;
                                    timeDelays = (-L : N)*Ts;
                                    tau(txID, rxID, mcID) = timeDelays(ind(rxID, mcID));
                            end
                        end
                    end
                end
            end
        end

        function err = get.timeOfArrivalErrorFromMatchFiltration(obj)
            err = obj.interface.timeDelay - obj.estimatedTimeDelaysFromMatchFiltration;
        end

        function idx = get.hypothesisTestingResults(obj)
            if obj.network.multiTransmitter
                idx = nan;
                warning('not implemented');
            else
                switch obj.network.networkCoherency
                    case "coherent"
                        switch obj.processingAlgorithm
                            case {1, 4}
                                z = real(obj.signalsIntegrated);
                                idx = cell(1, obj.monteCarlo.numberOfTrialsParallel);
                                for mcID = 1 : obj.monteCarlo.numberOfTrialsParallel
                                    idx{mcID} = find(z(:, mcID) > obj.configuration.threshold);
                                end
                            otherwise
                                error('Processing algorithm %d cannot be used for the current configuration', obj.processingAlgorithm);
                        end
                    case "short-term coherent"
                        idx = nan;
                        warning('not implemented');
                    case "incoherent"
                        idx = nan;
                        warning('not implemented');
                end
            end
        end

        function pos = get.estimatedPositionsFromIntegration(obj)
            gridScan = obj.gridPointsMesh;
            pos = cell(1, obj.monteCarlo.numberOfTrialsParallel);
            idx = obj.hypothesisTestingResults;
            for mcID = 1 : obj.monteCarlo.numberOfTrialsParallel
                pos{mcID} = [gridScan.x(idx{mcID}) gridScan.y(idx{mcID}) gridScan.z(idx{mcID})].';
            end
        end

        function err = get.positionError(obj)
            err = cell(1, obj.monteCarlo.numberOfTrialsParallel);
            for mcID = 1 : obj.monteCarlo.numberOfTrialsParallel
                estimatedPositions = obj.estimatedPositionsFromIntegration{mcID};
                numberOfEstimations = size(estimatedPositions, 2);
                err{mcID} = zeros(3, numberOfEstimations, obj.interface.numberOfTargets);
                for targetID = 1 : obj.interface.numberOfTargets
                    err{mcID}(:, :, targetID) = obj.interface.targets.position(:, targetID) - estimatedPositions;
                end
            end
        end

        function err = get.positionErrorTotal(obj)
            posErr = obj.positionError;
            err = cell(1, obj.monteCarlo.numberOfTrialsParallel);
            for mcID = 1 : obj.monteCarlo.numberOfTrialsParallel
                err{mcID} = sqrt(sum(posErr{mcID}.^2));
            end
        end

        %%% set methods

        function setintegrationindices(obj)
            obj.integrationIndices = ones(obj.network.numberOfActiveTransmittingNodes, obj.network.numberOfActiveReceivingNodes, prod(obj.gridSize));
            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                Ts = obj.network.activeReceivingNodes(rxID).samplingPeriod;
                for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                    switch obj.network.activeTransmittingNodes(txID).transmissionType
                        case "continuous"
                            error('not implemented');
                        case "pulsed"
                            L = obj.network.pulseWidthSample(rxID, txID) - 2;
                            N = obj.network.activeReceivingNodes(rxID).numberOfSamplesPerCPI;
                    end
                    timeDelays = permute((-L : N)*Ts, [1 3 4 2]);
                    [~, obj.integrationIndices(txID, rxID, :)] = min(abs(obj.hypothesizedTimeDelays(txID, rxID, :) - timeDelays), [], 4);
                end
            end
        end

        function setmatchfilteredsignals(obj)
            obj.signalsMatchFiltered = obj.signalsMatchFilteredTrials;
        end

        function configure(obj, options)
            arguments
                obj
                options.PFA (1, 1) double {mustBeNonnegative} = obj.configuration.PFA
                options.processingAlgorithm (1, 1) {mustBeInteger, mustBeInRange(options.processingAlgorithm, 1, 10)} = obj.processingAlgorithm
                options.numberOfTrials (1, 1) {mustBeNonnegative, mustBeInteger} = obj.monteCarlo.numberOfTrials
                options.numberOfTrialsParallel (1, 1) {mustBeNonnegative, mustBeInteger} = obj.monteCarlo.numberOfTrialsParallel
                options.seed (1, 1) {mustBeNonnegative, mustBeInteger, mustBeLessThan(options.seed, 4294967296)} = 0
                options.seedShuffle (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.seedShuffle, [0, 1])} = obj.seedShuffle
            end
            obj.configuration.PFA = options.PFA;
            obj.processingAlgorithm = options.processingAlgorithm;
            obj.monteCarlo.numberOfTrials = options.numberOfTrials;
            obj.monteCarlo.numberOfTrialsParallel = options.numberOfTrialsParallel;
            obj.monteCarlo.seed = options.seed;
            obj.seedShuffle = options.seedShuffle;
            rng(obj.monteCarlo.seed);
        end

        %%% monte carlo simulation

        function [PD, PFA] = simulatedetection(obj)
            arguments
                obj
            end
            gridScan = obj.gridPointsMesh;
            dims = size(gridScan.x) ~= 1;
            dimensions = {"y", "x", "z"};
            dimensions = dimensions(dims);
            numberOfDetections = zeros(obj.gridSize);
            numberOfTotalTrials = obj.monteCarlo.numberOfTrials*obj.monteCarlo.numberOfTrialsParallel;
            for mcID = 1 : obj.monteCarlo.numberOfTrials
                obj.setmatchfilteredsignals;
                idx = obj.hypothesisTestingResults;
                for mcpID = 1 : obj.monteCarlo.numberOfTrialsParallel
                    numberOfDetections(idx{mcpID}) = numberOfDetections(idx{mcpID}) + 1;
                end
            end
            probabilities = numberOfDetections./numberOfTotalTrials;
            targetIDx = find(sum(abs(obj.interface.targets.position - permute([gridScan.x(:) gridScan.y(:) gridScan.z(:)], [2 3 1])).^2) < sum(abs(obj.gridResolution(dims)).^2));
            noiseIDx = setdiff(1 : prod(obj.gridSize), targetIDx);
            PD = zeros(obj.gridSize([2 1 3])); PFA = zeros(obj.gridSize([2 1 3]));
            PD(targetIDx) = probabilities(targetIDx);
            PFA(noiseIDx) = probabilities(noiseIDx);
            xLabel = dimensions{1} + " (km)";
            yLabel = dimensions{2} + " (km)";
            x1 = gridScan.(dimensions{1})/1e3;
            x2 = gridScan.(dimensions{2})/1e3;
            figure; m = mesh(x2, x1, PD); colorbar;
            m.FaceColor = 'flat'; colorbar; hold on;
            x = obj.interface.targets.position(dims, :)/1e3;
            plot3(x(1, :), x(2, :), ones(1, size(x, 2)), '+k', 'LineWidth', 2);
            grid on; grid minor; view(0, 90);
            xlabel(xLabel); ylabel(yLabel); zlabel('P_D');
            title('Probability of detection'); hold off;
            figure; m = mesh(x2, x1, PFA); colorbar;
            m.FaceColor = 'flat'; colorbar;
            grid on; grid minor; view(0, 90);
            xlabel(xLabel); ylabel(yLabel); zlabel('P_{FA}');
            title('Probability of false alarm'); hold off;
        end
        %%% visualization methods

        function visualizefilteredsignals(obj, options)
            arguments
                obj
                options.receivingNodeIDs (1, :) double {mustBeInteger, mustBeNonnegative} = 1 : obj.network.numberOfActiveReceivingNodes
            end
            mustBeInRange(options.receivingNodeIDs, 1, obj.network.numberOfReceivingNodes);
            mcID = 1;
            figure; hold on;
            if isscalar(obj.signalsMatchFiltered) || isempty(obj.signalsMatchFiltered)
                obj.setmatchfilteredsignals;
            end
            s = obj.signalsMatchFiltered(:, options.receivingNodeIDs, mcID); % first trial
            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                Ts = obj.network.activeReceivingNodes(rxID).samplingPeriod;
                for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                    switch obj.network.activeTransmittingNodes(txID).transmissionType
                        case "continuous"
                            error('not implemented');
                        case "pulsed"
                            L = obj.network.pulseWidthSample(rxID, txID) - 2;
                            N = obj.network.activeReceivingNodes(rxID).numberOfSamplesPerCPI;
                    end
                    t = (-L : N)*Ts;
                    plot(t, 20*log10(abs(s(:, rxID))));
                end
            end
            grid on; grid minor;
            xlabel('time (s)'); ylabel('power (dB)');
            title('filtered signal');
            leg = legend(num2str(options.receivingNodeIDs.'), 'Location', 'best');
            title(leg, 'RX ID');
        end

        function visualizehypothesizedtimedelays(obj, options)
            arguments
                obj
                options.transmittingNodeID (1, 1) double {mustBeInteger, mustBePositive} = 1
                options.receivingNodeID (1, 1) double {mustBeInteger, mustBePositive} = 1
                options.dimension (1, 1) string {mustBeMember(options.dimension, ["x-y", "y-z", "z-x", "x-y-z"])} = "x-y"
            end
            mustBeInRange(options.transmittingNodeID, 1, obj.network.numberOfActiveTransmittingNodes);
            mustBeInRange(options.receivingNodeID, 1, obj.network.numberOfActiveReceivingNodes);
            gridScan = obj.gridPointsMesh;
            tau = obj.hypothesizedTimeDelays(options.transmittingNodeID, options.receivingNodeID, :)*1e6;
            figure;
            switch options.dimension
                case "x-y"
                    tau = reshape(squeeze(tau), [obj.gridSize]);
                    x = gridScan.x(:, :, 1); x = x(:)/1e3;
                    y = gridScan.y(:, :, 1); y = y(:)/1e3;
                    t = tau(:, :, 1);
                    plot3(x, y, t(:), '.'); hold on;
                    for zID = unique(ceil(linspace(2, size(tau, 3), 7)))
                        t = tau(:, :, zID);
                        plot3(x, y, t(:), '.');
                    end
                    grid on; grid minor;
                    xlabel('x (km)'); ylabel('y (km');
                    zlabel('time delay (\mus)');
                    title('hypothesized time delays');
                case "y-z"
                    tau = reshape(squeeze(tau), [obj.gridSize]);
                    y = gridScan.y(:, 1, :); y = y(:)/1e3;
                    z = gridScan.z(:, 1, :); z = z(:)/1e3;
                    t = tau(:, 1, :);
                    plot3(y, z, t(:), '.'); hold on;
                    for xID = unique(ceil(linspace(2, size(tau, 2), 7)))
                        t = tau(:, xID, :);
                        plot3(y, z, t(:), '.');
                    end
                    grid on; grid minor;
                    xlabel('y (km)'); ylabel('z (km)');
                    zlabel('time delay (\mus)');
                    title('hypothesized time delays');
                case "z-x"
                    tau = reshape(squeeze(tau), [obj.gridSize]);
                    z = gridScan.z(1, :, :); z = z(:)/1e3;
                    x = gridScan.x(1, :, :); x = x(:)/1e3;
                    t = tau(1, :, :);
                    plot3(z, x, t(:), '.'); hold on;
                    for yID = unique(ceil(linspace(2, size(tau, 1), 7)))
                        t = tau(yID, :, :);
                        plot3(z, x, t(:), '.');
                    end
                    grid on; grid minor;
                    xlabel('y (km)'); ylabel('z (km)');
                    zlabel('time delay (\mus)');
                    title('hypothesized time delays');
                case "x-y-z"
                    x = gridScan.x(:);
                    y = gridScan.y(:);
                    z = gridScan.z(:);
                    scatter3(x, y, z, 0.25, squeeze(tau));
                    hcb = colorbar;
                    title(hcb, 'time delay (\mus)');
                    grid on; grid minor;
                    xlabel('x (meter)'); ylabel('y (meter)'); zlabel('z (meter)');
                    title('hypothesized locations');
            end
        end

        function visualizeintegrationindices(obj)
            figure; plot(squeeze(obj.integrationIndices).', '.');
            xlabel('grid point'); ylabel('index');
            xlim tight; ylim tight;
            grid on; grid minor;
        end

        function visualizeintegratedsignals(obj, options)
            arguments
                obj
                options.scale (1, 1) string {mustBeMember(options.scale, ["linear", "power", "dB", "SNR"])} = "dB"
                options.plot (1, 1) string {mustBeMember(options.plot, ["magnitude", "real", "imaginary"])} = "real"
            end
            if isnan(obj.integrationWeights)
                return;
            end
            figure; hold on;
            mcID = 1;
            pos = obj.estimatedPositionsFromIntegration(:, mcID); % first trial
            gridScan = obj.gridPointsMesh;
            dims = size(gridScan.x) ~= 1;
            dimensions = {"y", "x", "z"};
            dimensions = dimensions(dims);
            switch sum(dims)
                case 1
                    xLabel = dimensions{1} + " (km)";
                    x1 = gridScan.(dimensions{1})/1e3;
                    plot(x1, 20*log10(abs(obj.signalsIntegrated(:, mcID)))); % first trial
                    grid on; grid minor;
                    xlabel(xLabel); ylabel('power (dB)');
                    title('integrated signal');
                case 2
                    switch options.plot
                        case "magnitude"
                            zLabel = 'magnitude';
                            Y = reshape(abs(obj.signalsIntegrated(:, mcID)), obj.gridSize([2, 1, 3])); % first trial
                        case "real"
                            zLabel = 'real part';
                            Y = reshape(real(obj.signalsIntegrated(:, mcID)), obj.gridSize([2, 1, 3])); % first trial
                        case "imaginary"
                            zLabel = 'imaginary part';
                            Y = reshape(imag(obj.signalsIntegrated(:, mcID)), obj.gridSize([2, 1, 3])); % first trial
                    end
                    switch options.scale
                        case "linear"
                            zLabel = [zLabel, ', linear'];
                        case "power"
                            zLabel = [zLabel, ', power'];
                            Y = Y.^2;
                        case "dB"
                            zLabel = [zLabel, ', power (dB)'];
                            Y = 20*log10(abs(Y));
                        case "SNR"
                            zLabel = [zLabel, ', SNR (dB)'];
                            Y = 20*log10(abs(Y)) - unique([obj.interface.network.activeReceivingNodes.noisePower]);
                    end
                    maxY = max(Y, [], 'all');
                    xLabel = dimensions{1} + " (km)";
                    yLabel = dimensions{2} + " (km)";
                    x1 = gridScan.(dimensions{1})/1e3;
                    x2 = gridScan.(dimensions{2})/1e3;
                    posRX = [obj.network.activeReceivingNodes.position]/1e3;
                    posTX = [obj.network.activeTransmittingNodes.position]/1e3;
                    posRX = posRX(dims, :); posTX = posTX(dims, :);
                    plot3(posRX(1, :), posRX(2, :), repmat(maxY, [1, size(posRX, 2)]), 'vb', 'LineWidth', 2);
                    plot3(posTX(1, :), posTX(2, :), repmat(maxY, [1, size(posTX, 2)]), 'vr', 'LineWidth', 2);
                    for targetID = 1 : length(pos)
                        estimations = pos{targetID}(dims, :)/1e3;
                        plot3(estimations(1, :), estimations(2, :), repmat(maxY, [1, size(estimations, 2)]), 'om', 'LineWidth', 2);
                    end
                    x = obj.interface.targets.position(dims, :)/1e3;
                    plot3(x(1, :), x(2, :), repmat(maxY, [1, size(x, 2)]), '+k', 'LineWidth', 2);
                    plot3(x2, x1, obj.configuration.threshold_dB*ones(size(x2)));
                    m = mesh(x2, x1, Y);
                    m.FaceColor = 'flat'; colorbar;
                    grid on; grid minor; view(0, 90);
                    xlabel(xLabel); ylabel(yLabel); zlabel(zLabel);
                    title('integrated signal'); hold off;
                case 3
            end
            legend('RX', 'TX', 'estimations', 'targets', 'Location', 'best');
        end

        function visualizeestimation(obj)
            mcID = 1;
            posRX = [obj.network.activeReceivingNodes.position]/1e3;
            posTX = [obj.network.activeTransmittingNodes.position]/1e3;
            figure; plot3(posRX(1, :), posRX(2, :), posRX(3, :), 'vb', 'LineWidth', 2, 'MarkerSize', 10);
            hold on; plot3(posTX(1, :), posTX(2, :), posTX(3, :), 'vr', 'LineWidth', 2, 'MarkerSize', 10);
            x = obj.interface.targets.position(1, :)/1e3;
            y = obj.interface.targets.position(2, :)/1e3;
            z = obj.interface.targets.position(3, :)/1e3;
            plot3(x, y, z, '*k', 'LineWidth', 1, 'MarkerSize', 10);
            color = ['m', 'g', 'c'];
            pos = obj.estimatedPositionsFromIntegration(:, mcID);
            for targetID = 1 : obj.configuration.numberOfTargets
                xEstimation = pos{targetID}(1, :)/1e3;
                yEstimation = pos{targetID}(2, :)/1e3;
                zEstimation = pos{targetID}(3, :)/1e3;
                if size(xEstimation, 2) < 9
                    plot3(xEstimation, yEstimation, zEstimation, ['+', color(min(targetID, 3))], 'LineWidth', 2, 'MarkerSize', 10);
                else
                    plot3(xEstimation, yEstimation, zEstimation, ['*', color(min(targetID, 3))], 'LineWidth', 1, 'MarkerSize', 10);
                end
            end
            text(posRX(1, :), posRX(2, :), posRX(3, :), num2str([obj.network.activeReceivingNodes.id].'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            text(posTX(1, :), posTX(2, :), posTX(3, :), num2str([obj.network.activeTransmittingNodes.id].'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            if size(posRX, 2) < 11
                text(x, y, z, num2str((1 : obj.interface.numberOfTargets).'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            end
            posRX = repelem(posRX, 1, obj.interface.network.numberOfActiveTransmittingNodes);
            posTX = repmat(posTX, 1, obj.interface.network.numberOfActiveReceivingNodes);
            line([posRX(1, :); posTX(1, :)], [posRX(2, :); posTX(2, :)], [posRX(3, :); posTX(3, :)], 'lineStyle', '--', 'Color', 'k');
            grid on; grid minor; title('scenario');
            xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
            legend('RX', 'TX', 'targets', 'estimations', 'Location', 'best');
        end

        function visualizereceiveroperatingcharacteristics(obj, options)
            arguments
                obj
                options.snr_dB (1, :) double = [0 3 10 13] % dB
            end
            nop = 1000;
            PFA = logspace(-10, 0, nop);
            PD = zeros(nop, length(options.snr_dB));
            for k = 1 : length(options.snr_dB)
                PD(:, k) = obj.ROC(PFA, 10.^(.1*options.snr_dB(k)));
            end
            figure; semilogx(PFA, PD);
            grid on; grid minor;
            xlabel('PFA'); ylabel('PD');
            leg = legend(num2str(options.snr_dB.'), 'Location', 'best');
            title(leg, 'SNR');
            title("Receiver operating characteristic (ROC) curve of \color{red}" + obj.network.networkCoherency + "\color{black} network with algorithm \color{blue}L_{" + num2str(obj.processingAlgorithm) + "}");
        end
    end
end