classdef radarNetwork < handle
    %radarNetwork Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = private, GetAccess = public)
        receivingNodes (1, :) receivingNode = receivingNode.empty() % (1 x Nrx vector)
        transmittingNodes (1, :) transmittingNode = transmittingNode.empty() % (1 x Ntx vector)
        receivingNodeActivity (1, :) logical = false % (1 x Nrx vector)
        transmittingNodeActivity (1, :) logical = false % (1 x Ntx vector)
        networkCoherency (1, 1) string {mustBeMember(networkCoherency, ["coherent", "short-term coherent", "incoherent"])} = "coherent"
    end

    properties (Dependent)
        numberOfReceivingNodes (1, 1) double {mustBeInteger, mustBeNonnegative}
        numberOfTransmittingNodes (1, 1) double {mustBeInteger, mustBeNonnegative}
        numberOfBistaticPairs (1, 1) double {mustBeInteger, mustBeNonnegative}
        numberOfActiveReceivingNodes (1, 1) double {mustBeInteger, mustBeNonnegative}
        numberOfActiveTransmittingNodes (1, 1) double {mustBeInteger, mustBeNonnegative}
        numberOfActiveBistaticPairs (1, 1) double {mustBeInteger, mustBeNonnegative}
        multiReceiver (1, 1) logical
        multiTransmitter (1, 1) logical
        passiveNetwork (1, 1) logical
        activeReceivingNodes (1, :) receivingNode % (1 x Nrx vector)
        activeTransmittingNodes (1, :) transmittingNode % (1 x Ntx vector)
        positionReceivingNode (3, :, :) double % (3 x Ntx x Nrx matrix)
        positionTransmittingNode (3, :, :) double % (3 x Ntx x Nrx matrix)
        distanceBaseline double {mustBePositive} % (Ntx x Nrx matrix)
        directPathDelay double {mustBePositive} % (Ntx x Nrx matrix)
        center (3, 1) double
        boundary (3, 2) double
        boundaryListened (3, 2) double
        defaultGridResolution (3, 1) double {mustBeNonnegative}
        pulseWidthSample double % (Nrx x Ntx matrix)
        dutyCycles double % (Nrx x Nrx matrix)
        matchFilter double % (L x Nrx x Ntx matrix)
        processingGain_dB double % (Nrx x Nrx matrix)
        noisePSDmatrix double % W/Hz (Nrx x Nrx matrix)
    end

    properties (Constant)
        speedOfLight (1, 1) double {mustBePositive} = physconst('LightSpeed');
    end

    methods
        function obj = radarNetwork(options)
            arguments
                options.receivingNodes (1, :) receivingNode = receivingNode.empty()
                options.transmittingNodes (1, :) transmittingNode = transmittingNode.empty()
                options.networkCoherency (1, 1) string {mustBeMember(options.networkCoherency, ["coherent", "short-term coherent", "incoherent"])} = "coherent"
            end
            %radarNetwork Construct an instance of this class
            %   Detailed explanation goes here
            obj.receivingNodes = options.receivingNodes;
            obj.transmittingNodes = options.transmittingNodes;
            obj.transmittingNodeActivity = true(1, obj.numberOfTransmittingNodes);
            obj.receivingNodeActivity = true(1, obj.numberOfReceivingNodes);
            obj.networkCoherency = options.networkCoherency;
        end

        %%% get methods

        function N = get.numberOfReceivingNodes(obj)
            N = numel(obj.receivingNodes);
        end

        function N = get.numberOfTransmittingNodes(obj)
            N = numel(obj.transmittingNodes);
        end

        function N = get.numberOfBistaticPairs(obj)
            N = obj.numberOfReceivingNodes*obj.numberOfTransmittingNodes;
        end

        function N = get.numberOfActiveReceivingNodes(obj)
            N = nnz(obj.receivingNodeActivity);
        end

        function N = get.numberOfActiveTransmittingNodes(obj)
            N = nnz(obj.transmittingNodeActivity);
        end

        function N = get.numberOfActiveBistaticPairs(obj)
            N = obj.numberOfActiveReceivingNodes*obj.numberOfActiveTransmittingNodes;
        end

        function tf = get.multiReceiver(obj)
            tf = logical(obj.numberOfActiveReceivingNodes - 1);
        end

        function tf = get.multiTransmitter(obj)
            tf = logical(obj.numberOfActiveTransmittingNodes - 1);
        end

        function tf = get.passiveNetwork(obj)
            tf = ~logical(obj.numberOfActiveTransmittingNodes);
        end

        function rx = get.activeReceivingNodes(obj)
            rx = obj.receivingNodes(obj.receivingNodeActivity);
        end

        function tx = get.activeTransmittingNodes(obj)
            tx = obj.transmittingNodes(obj.transmittingNodeActivity);
        end

        function pos = get.positionReceivingNode(obj)
            pos = repelem([obj.activeReceivingNodes.position], 1, obj.numberOfActiveTransmittingNodes);
            pos = reshape(pos, 3, obj.numberOfActiveTransmittingNodes, obj.numberOfActiveReceivingNodes);
        end

        function pos = get.positionTransmittingNode(obj)
            pos = repmat([obj.activeTransmittingNodes.position], 1, obj.numberOfActiveReceivingNodes);
            pos = reshape(pos, 3, obj.numberOfActiveTransmittingNodes, obj.numberOfActiveReceivingNodes);
        end

        function D = get.distanceBaseline(obj)
            D = sqrt(sum((obj.positionReceivingNode - obj.positionTransmittingNode).^2));
            D = reshape(D, obj.numberOfActiveTransmittingNodes, obj.numberOfActiveReceivingNodes);
        end

        function tau = get.directPathDelay(obj)
            tau = obj.distanceBaseline./obj.speedOfLight;
        end

        function origin = get.center(obj)
            origin = mean([obj.activeTransmittingNodes.position, obj.activeReceivingNodes.position], 2);
        end

        function lims = get.boundary(obj)
            lims = zeros(3, 2);
            nodePositions = [obj.activeTransmittingNodes.position, obj.activeReceivingNodes.position];
            lims(:, 1) = min(nodePositions, [], 2);
            lims(:, 2) = max(nodePositions, [], 2);
        end

        function lims = get.boundaryListened(obj)
            lims = obj.boundary;
            offset = max([obj.activeReceivingNodes.listenedRadius]);
            lims(:, 1) = lims(:, 1) - offset;
            lims(:, 2) = lims(:, 2) + offset;
        end

        function res = get.defaultGridResolution(obj)
            res = diff(obj.boundaryListened, [], 2)/100;
        end

        function N = get.pulseWidthSample(obj)
            % (Nrx x Ntx matrix)
            N = [obj.activeTransmittingNodes.pulseWidth].*[obj.activeReceivingNodes.samplingFrequency].';
        end

        function D = get.dutyCycles(obj)
            D = obj.pulseWidthSample./[obj.activeReceivingNodes.numberOfSamplesPerCPI].';
        end

        function mf = get.matchFilter(obj)
            % (L x Nrx x Ntx matrix)
            Ts = [obj.activeReceivingNodes.samplingPeriod]; % 1 x Nrx matrix
            N = obj.pulseWidthSample; % Nrx x Ntx matrix
            Nunique = unique(N);
            switch length(Nunique)
                case 1
                    switch obj.activeTransmittingNodes(1).transmissionType
                        case "continuous" % infinite pulse width
                            t = [obj.receivingNodes.samplingInstants];
                            mf = ones(size(t, 1), obj.numberOfActiveReceivingNodes, obj.numberOfActiveTransmittingNodes);
                        case "pulsed"
                            t = zeros(Nunique, obj.numberOfActiveReceivingNodes);
                            mf = zeros(Nunique, obj.numberOfActiveReceivingNodes, obj.numberOfActiveTransmittingNodes);
                            for txID = 1 : obj.numberOfActiveTransmittingNodes
                                for rxID = 1 : obj.numberOfActiveReceivingNodes
                                    t(:, rxID) = ((1 - N(rxID, txID) : 0)*Ts(rxID));
                                end
                                mf(:, :, txID) = obj.activeTransmittingNodes(txID).waveform(t);
                            end
                    end
                    mf = mf./sqrt(sum(abs(mf).^2, 1));
                otherwise
                    error('not implemented');
            end
        end

        function G = get.processingGain_dB(obj)
            G = permute(20*log10(sum(abs(obj.matchFilter))), [2 3 1]);
        end

        function B = get.noisePSDmatrix(obj)
            B = diag(10.^(.1*[obj.activeReceivingNodes.noisePowerPerSample_dB])./[obj.activeReceivingNodes.samplingFrequency]);
        end

        %%% set methods

        function activitycontrol(obj, options)
            arguments
                obj
                options.transmittingNodeActivity (1, :) logical = obj.transmittingNodeActivity
                options.receivingNodeActivity (1, :) logical = obj.receivingNodeActivity
            end
            obj.transmittingNodeActivity = options.transmittingNodeActivity;
            obj.receivingNodeActivity = options.receivingNodeActivity;
        end

        %%% visualizeton methods

        function visualizewaveformsampled(obj, options)
            arguments
                obj
                options.domain (1, 1) string {mustBeMember(options.domain, ["ambiguity", "frequency", "time"])} = "ambiguity"
                options.axisAmbiguity (1, 1) string {mustBeMember(options.axisAmbiguity, ["3D", "zeroDoppler", "zeroRange"])} = "zeroDoppler"
                options.plot (1, 1) string {mustBeMember(options.plot, ["magnitude", "phase", "real", "imaginary"])} = "magnitude"
            end
            minMag = -30; % dB
            N = unique(obj.pulseWidthSample);
            nfft = 2.^(nextpow2(N) + 1);
            bins = -nfft/2 : nfft/2 - 1;
            delays = -N + 1 : N - 1;
            if ~strcmpi(options.domain, "ambiguity")
                options.axisAmbiguity = "";
            end
            mfAll = obj.matchFilter;
            figure; hold on;
            if strcmpi(options.axisAmbiguity, "3D")
                txIDs = obj.numberOfActiveTransmittingNodes;
            else
                txIDs = 1 : obj.numberOfActiveTransmittingNodes;
            end
            for txID = txIDs
                for rxID = 1 : obj.numberOfActiveReceivingNodes
                    if strcmpi(options.axisAmbiguity, "3D") && rxID ~= 1
                        break;
                    end
                    mf = mfAll(:, rxID, txID);
                    switch options.domain
                        case "ambiguity"
                            options.plot = "magnitude";
                            switch options.axisAmbiguity
                                case "zeroDoppler"
                                    af = xcorr(mf);
                                    plot(delays, 10*log10(abs(af)));
                                case "zeroRange"
                                    af = fftshift(fft(abs(mf).^2, nfft));
                                    plot(bins, 10*log10(abs(af)));
                                case "3D"
                                    af = zeros(2*N - 1, nfft);
                                    for k = 1 : nfft(1)
                                        af(:, k) = xcorr(mf.*exp(1j*2*pi*bins(k)*(0 : N - 1).'/nfft), mf);
                                    end
                                    af = 10*log10(abs(af));
                                    af(af < minMag) = minMag;
                                    m = mesh(bins, delays, af); colorbar;
                            end
                        case "time"
                            switch options.plot
                                case "magnitude"
                                    plot(0 : N - 1, abs(mf));
                                case "phase"
                                    plot(0 : N - 1, 180*unwrap(angle(mf))/pi);
                                case "real"
                                    plot(0 : N - 1, real(mf));
                                case "imaginary"
                                    plot(0 : N - 1, imag(mf));
                            end
                        case "frequency"
                            switch options.plot
                                case "magnitude"
                                    plot(bins, 20*log10(fftshift(abs(fft(mf, nfft)))));
                                case "phase"
                                    plot(bins, 180*unwrap(angle(fftshift(fft(mf, nfft))))/pi);
                                case "real"
                                    plot(bins, 20*log10(fftshift(abs(fft(mf, nfft)))));
                                case "imaginary"
                                    plot(bins, 20*log10(fftshift(abs(fft(mf, nfft)))));
                            end
                    end
                end
            end
            switch options.domain
                case "ambiguity"
                    title('ambiguity functions');
                case "frequency"
                    title('waveforms in time');
                    xlabel('frequency bins');
                case "time"
                    title('frequency responses of waveform');
                    xlabel('time samples');
            end
            grid on; grid minor;
            if ~strcmpi(options.axisAmbiguity, "3D")
                if strcmpi(options.axisAmbiguity, "zeroDoppler")
                    xlabel('delay samples');
                elseif strcmpi(options.domain, "ambiguity")
                    xlabel('doppler bins');
                end
                switch options.plot
                    case "magnitude"
                        ylabel('power (dB)');
                    case "phase"
                        ylabel('phase (°)');
                    case "real"
                        ylabel('real part');
                    case "imaginary"
                        ylabel('imaginary part');
                end
                [idTx, idRx] = meshgrid(1 : obj.numberOfActiveTransmittingNodes, 1 : obj.numberOfActiveReceivingNodes);
                leg = legend(num2str([idTx(:) idRx(:)]), 'Location', 'best');
                title(leg, 'TX | RX ID');
            else
                m.FaceColor = 'flat';
                xlabel('doppler bins'); ylabel('delay samples');
                switch options.plot
                    case "magnitude"
                        zlabel('power (dB)');
                    case "phase"
                        zlabel('phase (°)');
                    case "real"
                        zlabel('real part');
                    case "imaginary"
                        zlabel('imaginary part');
                end
                zlim([minMag 0]);
            end
            hold off; drawnow;
        end

        function visualizenetwork(obj)
            posRX = [obj.receivingNodes.position]/1e3;
            posTX = [obj.transmittingNodes.position]/1e3;
            figure; plot3(posRX(1, :), posRX(2, :), posRX(3, :), 'vb', 'LineWidth', 3);
            hold on; plot3(posTX(1, :), posTX(2, :), posTX(3, :), 'vr', 'LineWidth', 3);
            text(posRX(1, :), posRX(2, :), posRX(3, :), num2str([obj.receivingNodes.id].'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            text(posTX(1, :), posTX(2, :), posTX(3, :), num2str([obj.transmittingNodes.id].'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            posRX = repelem(posRX, 1, obj.numberOfTransmittingNodes);
            posTX = repmat(posTX, 1, obj.numberOfReceivingNodes);
            line([posRX(1, :); posTX(1, :)], [posRX(2, :); posTX(2, :)], [posRX(3, :); posTX(3, :)], 'lineStyle', '--', 'Color', 'k');
            grid on; grid minor; title('radar network configuration');
            xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
            legend('RX', 'TX', 'Location', 'best');
        end
    end
end