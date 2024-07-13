classdef spu < handle
    %spu Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = private, GetAccess = public)
        interface (1, :) interface = interface.empty()
        gridResolution (3, 1) double {mustBeNonnegative} = 100*ones(3, 1)
        configuration (1, 1) struct = struct('numberOfTargets', 1, 'PFA', 1e-6)
        processingAlgorithm (1, 1) {mustBeInteger, mustBeInRange(processingAlgorithm, 1, 10)} = 1
    end

    properties (Dependent)
        network (1, 1) radarNetwork
        gridPoints cell
        gridSize (1, 3) double {mustBeNonnegative}
        gridPointsMesh struct
        hypothesizedTimeDelays double
    end

    properties (SetAccess = private, GetAccess = public)
        integrationIndices double {mustBeInteger, mustBePositive}
    end

    properties (Dependent)
        outputSNR double % dB power
        expectedAmplitudes double {mustBePositive} % linear scale magnitude normalized to first receiver (reference)
        noisePowersPerSample double {mustBePositive} % linear scale power
        integrationWeights double {mustBePositive}
        signalsMatchFiltered double % (Ns + L - 1 x Nrx matrix)
        signalsIntegrated double

        estimatedTimeDelaysFromMatchFiltration double
        timeOfArrivalErrorFromMatchFiltration double
        estimatedPositionsFromIntegration cell
        positionError double
        positionErrorTotal double
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

        function net = get.network(obj)
            net = obj.interface.network;
        end

        function cfg = get.configuration(obj)
            cfg = obj.configuration;
            % cfg.PD = cfg.PFA.^(1/(1 + 10.^(.1*obj.outputSNR)));
            % cfg.threshold = 10*log10(log(1./cfg.PFA.^2)) + 10*log10(sum(obj.noisePowers./(2*obj.expectedAmplitudeRatios.^2))); 
            if obj.network.multiTransmitter
                cfg.threshold = nan;
            else
                switch obj.network.networkCoherency
                    case "coherent"
                        switch obj.processingAlgorithm
                            case 1
                                cfg.threshold = 2*sqrt(obj.outputSNR)*erfcinv(2*cfg.PFA);
                            otherwise
                                error('Processing algorithm %d cannot be used for the current configuration', obj.processingAlgorithm);
                        end
                    case "short-term coherent"
                        cfg.threshold = nan;
                    case "incoherent"
                        cfg.threshold = nan;
                end
            end
            cfg.thresholddB = 20*log10(abs(cfg.threshold));
            cfg.PD = erfc((sqrt(2)*cfg.threshold/sqrt(obj.outputSNR) - sqrt(obj.outputSNR))/2)/2;
        end

        function n = get.noisePowersPerSample(obj)
            n = 10.^(.1*[obj.network.activeReceivingNodes.noisePowerPerSample]); % kTBN
        end

        function w = get.expectedAmplitudes(obj)
            w = 10.^(.05*obj.interface.receivedPowerFromScatterers);
        end

        function w = get.integrationWeights(obj)
            if obj.network.multiTransmitter
                w = nan;
                warning('not implemented');
            else
                switch obj.network.networkCoherency
                    case "coherent"
                        switch obj.processingAlgorithm
                            case 1
                                if obj.interface.configuration.noise
                                    w = obj.expectedAmplitudes./obj.noisePowersPerSample;
                                else
                                    w = obj.expectedAmplitudes;
                                end
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

        function snr = get.outputSNR(obj)
            snr = sum(10.^(.1*obj.network.processingGain.').*obj.expectedAmplitudes.^2./obj.noisePowersPerSample, 2);
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

        function y = get.signalsMatchFiltered(obj)
            % (Ns + L - 1 x Nrx matrix)
            t = [obj.network.activeReceivingNodes.samplingInstants]; % Ns x Nrx matrix
            demodulator = exp(1j*2*pi*shiftdim([obj.network.activeTransmittingNodes.carrierFrequency], -1).*t);
            mf = obj.network.matchFilter; % L x Nrx x Ntx matrix
            s = obj.interface.signalSuperposed; % Ns x Nrx matrix
            y = zeros(size(s, 1) + size(mf, 1) - 1, obj.network.numberOfActiveReceivingNodes);
            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                y(:, rxID) = conv(s(: ,rxID).*demodulator(:, rxID), conj((mf(:, rxID, 1))));
            end
        end

        function Y = get.signalsIntegrated(obj)
            % (Ns + L - 1 x 2*(Ns + L - 1) - 1 matrix)
            if obj.network.multiTransmitter
                Y = nan;
                warning('not implemented');
            else
                switch obj.network.networkCoherency
                    case "coherent"
                        switch obj.processingAlgorithm
                            case 1
                                yWeighted = obj.integrationWeights.*obj.signalsMatchFiltered;
                                Y = zeros(prod(obj.gridSize), 1);
                                for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                                    Y = Y + yWeighted(obj.integrationIndices(1, rxID, :), rxID);
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
                [~, ind] = max(abs(obj.signalsMatchFiltered));
                tau = zeros(obj.network.numberOfActiveTransmittingNodes, obj.network.numberOfActiveReceivingNodes);
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
                                tau(txID, rxID) = timeDelays(ind(rxID));
                        end
                    end
                end
            end
        end

        function err = get.timeOfArrivalErrorFromMatchFiltration(obj)
            err = obj.interface.timeDelay - obj.estimatedTimeDelaysFromMatchFiltration;
        end

        function pos = get.estimatedPositionsFromIntegration(obj)
            if obj.network.multiTransmitter
                pos = nan;
                warning('not implemented');
            else
                switch obj.network.networkCoherency
                    case "coherent"
                        switch obj.processingAlgorithm
                            case 1
                                gridScan = obj.gridPointsMesh;
                                z = abs(obj.signalsIntegrated); % real?
                                pos = cell(1, obj.configuration.numberOfTargets);
                                for targetID = 1 : obj.configuration.numberOfTargets
                                    idx = find(z > obj.configuration.threshold);
                                    pos{targetID} = [gridScan.x(idx) gridScan.y(idx) gridScan.z(idx)].';
                                end
                            otherwise
                                error('Processing algorithm %d cannot be used for the current configuration', obj.processingAlgorithm);
                        end
                    case "short-term coherent"
                        pos = nan;
                        warning('not implemented');
                    case "incoherent"
                        pos = nan;
                        warning('not implemented');
                end
            end
        end

        function err = get.positionError(obj)
            err = zeros(3, obj.interface.numberOfTargets, obj.configuration.numberOfTargets);
            for estimationID = 1 : obj.configuration.numberOfTargets
                for targetID = 1 : obj.interface.numberOfTargets
                    err(:, targetID, estimationID) = obj.interface.targets.position(:, targetID) - obj.estimatedPositionsFromIntegration{estimationID}(:, 1 : obj.configuration.numberOfTargets);
                end
            end
        end

        function err = get.positionErrorTotal(obj)
            err = sqrt(sum(obj.positionError.^2));
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

        function setconfiguration(obj, options)
            arguments
                obj
                options.numberOfTargets (1, 1) double {mustBeInteger, mustBePositive} = obj.configuration.numberOfTargets
                options.PFA (1, 1) double {mustBeNonnegative} = obj.configuration.PFA
            end
            obj.configuration.numberOfTargets = options.numberOfTargets;
            obj.configuration.PFA = options.PFA;
        end

        %%% visualization methods

        function visualizefilteredsignals(obj, options)
            arguments
                obj
                options.receivingNodeIDs (1, :) double {mustBeInteger, mustBeNonnegative} = 1 : obj.network.numberOfActiveReceivingNodes
            end
            mustBeInRange(options.receivingNodeIDs, 1, obj.network.numberOfReceivingNodes);
            figure; hold on;
            s = obj.signalsMatchFiltered(:, options.receivingNodeIDs);
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

        % function visualizeintegratedsignals2(obj)
        %     Y = obj.signalIntegrated;
        %     Ts = obj.network.receivingNodes(1).samplingPeriod;
        %     switch obj.network.transmittingNodes(1).transmissionType
        %         case "continuous"
        %             N = 2*unique([obj.network.receivingNodes.numberOfSamples]) - 1;
        %         case "pulsed"
        %             N = unique([obj.network.receivingNodes.numberOfSamples]) + unique([obj.network.pulseWidthSample]) - 1;
        %     end
        %     t1 = (1 : N)*Ts;
        %     t2 = (-N + 1 : N - 1)*Ts;
        %     figure; mesh(t2, t1, 20*log10(abs(Y)));
        %     grid on; grid minor; colorbar;
        %     xlabel('time difference (s)'); ylabel('time delay RX1 (s)'); zlabel('power (dB)');
        %     title('integrated signal');
        % end

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
            pos = obj.estimatedPositionsFromIntegration;
            gridScan = obj.gridPointsMesh;
            dims = size(gridScan.x) ~= 1;
            dimensions = {"y", "x", "z"};
            dimensions = dimensions(dims);
            switch sum(dims)
                case 1
                    xLabel = dimensions{1} + " (km)";
                    x1 = gridScan.(dimensions{1})/1e3;
                    plot(x1, 20*log10(abs(obj.signalsIntegrated)));
                    grid on; grid minor;
                    xlabel(xLabel); ylabel('power (dB)');
                    title('integrated signal');
                case 2
                    switch options.plot
                        case "magnitude"
                            zLabel = 'magnitude';
                            Y = reshape(abs(obj.signalsIntegrated), obj.gridSize([2, 1, 3]));
                        case "real"
                            zLabel = 'real part';
                            Y = reshape(real(obj.signalsIntegrated), obj.gridSize([2, 1, 3]));
                        case "imaginary"
                            zLabel = 'imaginary part';
                            Y = reshape(imag(obj.signalsIntegrated), obj.gridSize([2, 1, 3]));
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
                    x = obj.interface.targets.position(dims, :)/1e3;
                    plot3(x(1, :), x(2, :), repmat(maxY, [1, size(x, 2)]), '+k', 'LineWidth', 2);
                    for targetID = 1 : obj.configuration.numberOfTargets
                        estimations = pos{targetID}(dims, :)/1e3;
                        plot3(estimations(1, :), estimations(2, :), repmat(maxY, [1, size(estimations, 2)]), 'om', 'LineWidth', 2);
                    end
                    m = mesh(x2, x1, Y);
                    m.FaceColor = 'flat'; colorbar;
                    grid on; grid minor; view(0, 90);
                    xlabel(xLabel); ylabel(yLabel); zlabel(zLabel);
                    title('integrated signal'); hold off;
                case 3
            end
            legend('RX', 'TX', 'targets', 'estimations', 'Location', 'best');
        end

        function visualizeestimation(obj)
            posRX = [obj.network.activeReceivingNodes.position]/1e3;
            posTX = [obj.network.activeTransmittingNodes.position]/1e3;
            figure; plot3(posRX(1, :), posRX(2, :), posRX(3, :), 'vb', 'LineWidth', 2, 'MarkerSize', 10);
            hold on; plot3(posTX(1, :), posTX(2, :), posTX(3, :), 'vr', 'LineWidth', 2, 'MarkerSize', 10);
            x = obj.interface.targets.position(1, :)/1e3;
            y = obj.interface.targets.position(2, :)/1e3;
            z = obj.interface.targets.position(3, :)/1e3;
            if obj.interface.numberOfTargets < 11
                plot3(x, y, z, '+k', 'LineWidth', 2, 'MarkerSize', 10);
            else
                plot3(x, y, z, '*k', 'LineWidth', 1, 'MarkerSize', 10);
            end
            color = ['m', 'g', 'c'];
            pos = obj.estimatedPositionsFromIntegration;
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
            if obj.interface.numberOfTargets < 11
                text(x, y, z, num2str((1 : obj.interface.numberOfTargets).'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            end
            posRX = repelem(posRX, 1, obj.interface.network.numberOfActiveTransmittingNodes);
            posTX = repmat(posTX, 1, obj.interface.network.numberOfActiveReceivingNodes);
            line([posRX(1, :); posTX(1, :)], [posRX(2, :); posTX(2, :)], [posRX(3, :); posTX(3, :)], 'lineStyle', '--', 'Color', 'k');
            grid on; grid minor; title('scenario');
            xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
            legend('RX', 'TX', 'targets', 'estimations', 'Location', 'best');
        end
    end
end