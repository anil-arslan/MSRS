classdef spu < handle
    %spu Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = private, GetAccess = public)
        interfaces (1, :) interface = interface.empty()
        monteCarlo (1, 1) struct = struct( ...
            'seed', 0, ...
            'numberOfTrials', 1, ...
            'numberOfTrialsParallel', 1)
        gridResolution (3, 1) double {mustBeNonnegative} = 100*ones(3, 1)
        configuration (1, 1) struct = struct( ...
            'PFA', 1e-6, ...
            'numberOfTargets', nan, ...
            'referenceActiveReceiverID', 1, ...
            'removeBoundaryDetectionDF', 1, ...
            'removeBoundaryDetectionMF', 0, ...
            'interpolationDF', 1, ...
            'interpolationMF', 0, ...
            'groupingCloseness', 2, ... 
            'CFAR', 1, ... % For coherent monostatic L2
            'numberOfTrainingCells', 20, ...
            'numberOfGuardCells', 5);
        processingAlgorithm (1, 1) {mustBeInteger, mustBeInRange(processingAlgorithm, 1, 6)} = 1
        % 1: Coherent detector - Perfect knowledge of range
        %   Maximal ratio combining
        %   Deterministic signal
        % 2: Coherent detector - Fully correlated spatial observations
        %   Spatially fully correlated phase and amplitude
        %   For weak stochastic signal, it is optimal for any distribution
        %   Otherwise, it is optimal under AWGN
        %   For generalized LRT, it is optimal adaptive processing algorithm
        % 3: Noncoherent detector - Square-law envelope detector
        %   Weak stochastic signal
        %   Spatially independent phase, arbitrarily correlated amplitude
        %   For weak stochastic signal, it is optimal for any distribution
        %   Otherwise, it is an approximation of optimal processing algorithm under AWGN
        % 4: Noncoherent detector - Linear-law envelope detector
        %   Strong stochastic signal
        %   Spatially independent phase, fully correlated amplitude
        %   For generalized LRT, it is an approximation of optimal adaptive processing algorithm
        % 5: Noncoherent detector - Square-law envelope detector
        %   Neither strong nor weak stochastic signal
        %   Spatially independent phase and amplitude
        % 6: Noncoherent detector - Square-law envelope detector
        %   Strong stochastic signal
        %   Spatially independent phase and amplitude
        %   For generalized LRT, it is optimal adaptive processing algorithm
        detectionAlgorithm (1, 1) string {mustBeMember(detectionAlgorithm, ["thresholding", "peak", "CoSaMP", "OMP", "NCOMP"])} = "thresholding"
    end

    properties (Dependent)
        network (1, 1) radarNetwork
        gridPoints cell
        gridSize (1, 3) double {mustBeNonnegative}
        gridPointsMesh struct
        hypothesizedTimeDelays double % (Ntx x Nrx x Ni matrix)
        hypothesizedPhases double % (Ntx x Nrx x Ni matrix)
        hypothesizedPhasesReferenced double % (Ntx x Nrx x Ni matrix)
        hypothesizedAmplitudes double
        hypothesizedAmplitudesReferenced double
    end

    properties (SetAccess = private, GetAccess = public)
        blindZone logical {mustBeNumericOrLogical, mustBeMember(blindZone, [0, 1])} = []; 
        integrationIndices double {mustBeInteger, mustBePositive} % (Ntx x Nrx x Ni vector)
        integrationWeights double % (Ntx x Nrx x Nt matrix)
    end

    properties (Dependent)
        inputSNR_dB double % dB power, (Ntx x Nrx x Nt matrix)
        inputSNR_lin dobule % linear scale power, (Ntx x Nrx x Nt matrix)
        outputSNR_dB cell % dB power
        outputSNR_lin cell % linear scale power
        threshold function_handle
        ROC function_handle
        expectedAmplitudes double {mustBePositive} % (Ntx x Nrx x Nt matrix) linear scale signal magnitude
        expectedAmplitudesReferenced double {mustBePositive} % (Ntx x Nrx x Nt matrix) linear scale signal magnitude
        noisePowersPerSample_W double {mustBePositive} % linear scale noise power
    end

    properties (Dependent, Hidden)
        signalsMatchFilteredTrials cell % (1 x Nrx cell of Ns + L - 1 x Nch x Nmcp matrix) Ncmp : number of parallel trials
        signalsCompressedTrials double % (Ns_i x Ncmp matrix)
    end

    properties (SetAccess = private, GetAccess = public)
        signalsMatchFiltered cell = {0} % (1 x Nrx cell of Ns + L - 1 x Nch x Nmcp matrix) Ncmp : number of parallel trials
        dictionaryCompression double = [] % (Ns_i x Ni matrix)
        compressionReport struct = struct('complexAmplitudes', [], 'cellIDs', [], 'residualHistory', []);
    end

    properties (Dependent)
        signalsIntegrated double % (Ni x 1 x Ncmp matrix)
        signalsIntegratedDirectly double % (Ns_i x Ncmp matrix)

        thresholdCFAR cell % (1 x Nrx cell of Ns + L - 1 x Nch x Nmcp matrix)
        detectionFromMatchFiltration struct % (Nrx x Nmcp struct of 1 x Nd vectors) Nd : number of detections

        hypothesisTestingResults cell % (1 x Nrx cell of Nch x Nmcp cell of Nd x 1 vector) Nd : number of detections
        estimatedPositions cell % (1 x Nrx cell of Nch x Nmcp cell of 3 x Nd matrix)
        detectionFromIntegration struct % (Nrx x Nmcp struct of 1 x Nd vectors) Nd : number of detections

        positionError cell % (1 x Nrx cell of Nch x Nmcp cell of 3 x Nd x Nt matrix)
        positionErrorTotal cell % (1 x Nrx cell of Nch x Nmcp cell of Nd x Nt matrix)

        covarianceMatrixInitial (3, 3) double
    end

    properties (Access = private)
        coverageSimulationReport struct = struct('targetCellIDs', [], 'PD', [], 'PFA', [])
        saveResiduals (1, 1) logical {mustBeNumericOrLogical, mustBeMember(saveResiduals, [0, 1])} = 0
        seedShuffle (1, 1) logical {mustBeNumericOrLogical, mustBeMember(seedShuffle, [0, 1])} = 0
        isWait (1, 1) logical {mustBeNumericOrLogical, mustBeMember(isWait, [0, 1])} = 0;
        isStep (1, 1) logical {mustBeNumericOrLogical, mustBeMember(isStep, [0, 1])} = 0;
        isRestart (1, 1) logical {mustBeNumericOrLogical, mustBeMember(isRestart, [0, 1])} = 0;
        isSimulation (1, 1) logical {mustBeNumericOrLogical, mustBeMember(isSimulation, [0, 1])} = 1;
    end

    properties (Constant)
        speedOfLight (1, 1) double {mustBePositive} = physconst('LightSpeed');
    end

    methods
        function obj = spu(options)
            arguments
                options.interfaces (1, 1) interface
                options.gridResolution (3, 1) double {mustBeNonnegative} = zeros(3, 1)
            end
            %spu Construct an instance of this class
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

        function SNRin = get.inputSNR_lin(obj)
            % (Ntx x Nrx x Nt matrix)
            SNRin = 10.^(.1*obj.inputSNR_dB);
        end

        function SNRin = get.inputSNR_dB(obj)
            % (Ntx x Nrx x Nt matrix)
            SNRin = obj.interfaces.inputSNR_dB;
        end

        function SNRout = get.outputSNR_lin(obj)
            switch obj.network.networkMode
                case "multiStatic"
                    SNRout = sum(10.^(.1*(obj.network.processingGain_dB.' + obj.inputSNR_dB)), [1, 2]);
                    switch obj.processingAlgorithm
                        case 6
                            SNRout = SNRout./obj.network.numberOfActiveBistaticPairs;
                    end
                case "monoStatic"
                    SNRout = cell(1, obj.network.numberOfActiveReceivingNodes);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        SNRout{rxID} = 10.^(.1*(obj.network.processingGain_dB(rxID, obj.network.monoStaticTransmitterIDs(rxID)) + obj.inputSNR_dB(obj.network.monoStaticTransmitterIDs(rxID), rxID, :)));
                    end
            end
        end

        function SNRout = get.outputSNR_dB(obj)
            switch obj.network.networkMode
                case "multiStatic"
                    SNRout = 10*log10(obj.outputSNR_lin);
                case "monoStatic"
                    SNRout = cell(1, obj.network.numberOfActiveReceivingNodes);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        SNRout{rxID} = 10*log10(obj.outputSNR_lin{rxID});
                    end
            end
        end

        function PD = get.ROC(obj)
            switch obj.network.networkMode
                case "multiStatic"
                    switch obj.processingAlgorithm
                        case 1 % Coherent detector - Perfect knowledge of range
                            %   Maximal ratio combining
                            %   Deterministic signal
        
                            % EXACTLY CORRECT FORMULA
                            PD = @(PFA, SNR, N) erfc(erfcinv(2*PFA) - sqrt(SNR))/2;
        
                        case 2 % Coherent detector - Fully correlated spatial observations
                            %   Spatially fully correlated phase and amplitude
                            %   For weak stochastic signal, it is optimal for any distribution
                            %   Otherwise, it is optimal under AWGN
                            %   For generalized LRT, it is optimal adaptive processing algorithm
        
                            % EXACTLY CORRECT FORMULA
                                % If the received signal has the same spatial
                                % probability distributions
                            PD = @(PFA, SNR, N) marcumq(sqrt(2*SNR), sqrt(-2*log(PFA)));
                                % If the received signal has different but fully
                                % correlated spatial probability distributions
                            % PD = @(PFA, SNR) PFA.^(1./(1 + SNR));
        
                        case 3 % Noncoherent detector - Square-law envelope detector
                            %   Weak stochastic signal
                            %   Spatially independent phase, arbitrarily correlated amplitude
                            %   For weak stochastic signal, it is optimal for any distribution
                            %   Otherwise, it is an approximation of optimal processing algorithm under AWGN
        
                            PD = @(PFA, SNR, N) nan;
        
                        case 4 % Noncoherent detector - Linear-law envelope detector
                            %   Strong stochastic signal
                            %   Spatially independent phase, fully correlated amplitude
                            %   For generalized LRT, it is an approximation of optimal adaptive processing algorithm
        
                            PD = @(PFA, SNR, N) nan;
        
                        case 5 % Noncoherent detector - Square-law envelope detector
                            %   Neither strong nor weak stochastic signal
                            %   Spatially independent phase and amplitude
        
                            PD = @(PFA, SNR, N) nan;
        
                        case 6 % Noncoherent detector - Square-law envelope detector
                            %   Strong stochastic signal
                            %   Spatially independent phase and amplitude
                            %   For generalized LRT, it is optimal adaptive processing algorithm
        
                            switch obj.network.numberOfActiveBistaticPairs
                                case 1
                                    PD = @(PFA, SNR, N) exp(-obj.threshold(PFA, SNR, N)./(1 + SNR));
                                otherwise
                                    PD = @(PFA, SNR, N) exp(-obj.threshold(PFA, SNR, N)./(1 + SNR)).*sum((obj.threshold(PFA, SNR, N)./(1 + SNR)).^((0 : N - 1).')./factorial((0 : N - 1).'), 1);
                            end
                    end
                    if ~obj.network.passiveNetwork % Single transmitter networks
                    else % Passive networks
                    end
                case "monoStatic"
                    switch obj.processingAlgorithm
                        case 1
                            PD = @(PFA, SNR, N) erfc(erfcinv(2*PFA) - sqrt(SNR))/2;
                        case 2
                            PD = @(PFA, SNR, N) marcumq(sqrt(2*SNR), sqrt(-2*log(PFA)));
                        case {3, 4, 5, 6}
                            PD = @(PFA, SNR, N) nan;
                        otherwise
                            error('algorithm %d is not valid for monostatic network', obj.processingAlgorithm);
                    end
            end
        end

        function T = get.threshold(obj)
            switch obj.network.networkMode
                case "multiStatic"
                    switch obj.processingAlgorithm
                        case 1 % Coherent detector - Perfect knowledge of range
                            %   Maximal ratio combining
                            %   Deterministic signal
        
                            % EXACTLY CORRECT FORMULA
                            T = @(PFA, SNR, N) sqrt(SNR)*erfcinv(2*PFA);
        
                        case 2 % Coherent detector - Fully correlated spatial observations
                            %   Spatially fully correlated phase and amplitude
                            %   For weak stochastic signal, it is optimal for any distribution
                            %   Otherwise, it is optimal under AWGN
                            %   For generalized LRT, it is optimal adaptive processing algorithm
        
                            % EXACTLY CORRECT FORMULA
                            T = @(PFA, SNR, N) sqrt(-SNR*log(PFA));
        
                        case 3 % Noncoherent detector - Square-law envelope detector
                            %   Weak stochastic signal
                            %   Spatially independent phase, arbitrarily correlated amplitude
                            %   For weak stochastic signal, it is optimal for any distribution
                            %   Otherwise, it is an approximation of optimal processing algorithm under AWGN
        
                            T = @(PFA, SNR, N) nan;
        
                        case 4 % Noncoherent detector - Linear-law envelope detector
                            %   Strong stochastic signal
                            %   Spatially independent phase, fully correlated amplitude
                            %   For generalized LRT, it is an approximation of optimal adaptive processing algorithm
        
                            T = @(PFA, SNR, N) nan;
        
                        case 5 % Noncoherent detector - Square-law envelope detector
                            %   Neither strong nor weak stochastic signal
                            %   Spatially independent phase and amplitude
        
                            T = @(PFA, SNR, N) nan;
        
                        case 6 % Noncoherent detector - Square-law envelope detector
                            %   Strong stochastic signal
                            %   Spatially independent phase and amplitude
                            %   For generalized LRT, it is optimal adaptive processing algorithm
        
                            T = @(PFA, SNR, N) obj.incompleteGammaInverse(PFA, N);
        
                    end
                    if ~obj.network.passiveNetwork % Single transmitter networks
                    else % Passive networks
                    end
                case "monoStatic"
                    switch obj.processingAlgorithm
                        case 1
                            T = @(PFA, SNR, N) sqrt(SNR)*erfcinv(2*PFA);
                        case 2
                            T = @(PFA, SNR, N) sqrt(-SNR*log(PFA));
                        case {3, 4, 5, 6}
                            T = @(PFA, SNR, N) obj.incompleteGammaInverse(PFA, N);
                        otherwise
                            error('algorithm %d is not valid for monostatic network', obj.processingAlgorithm);
                    end
            end
        end

        function cfg = get.configuration(obj)
            cfg = obj.configuration;
            switch obj.network.networkMode
                case "multiStatic"
                    cfg.PD = obj.ROC(cfg.PFA, obj.outputSNR_lin, obj.network.numberOfActiveBistaticPairs);
                    cfg.threshold = obj.threshold(cfg.PFA, obj.outputSNR_lin, obj.network.numberOfActiveBistaticPairs);
                    switch obj.processingAlgorithm
                        case {1, 4} % linear
                            cfg.threshold = cfg.threshold(1);
                            cfg.threshold_dB = 20*log10(abs(cfg.threshold));
                        case {2, 3, 5, 6} % square
                            cfg.threshold_dB = 10*log10(abs(cfg.threshold));
                    end
                case "monoStatic"
                    cfg.PD = zeros(1, obj.interfaces.numberOfReceivingNodes);
                    cfg.threshold = zeros(1, obj.interfaces.numberOfReceivingNodes);
                    for rxID = 1 : obj.interfaces.numberOfReceivingNodes
                        cfg.PD(rxID) = obj.ROC(cfg.PFA, mean(obj.outputSNR_lin{rxID}, 3), 1);
                        cfg.threshold(rxID) = obj.threshold(cfg.PFA, mean(obj.outputSNR_lin{rxID}, 3), 1);
                    end
                    switch obj.processingAlgorithm
                        case 1 % linear
                            cfg.threshold_dB = 20*log10(abs(cfg.threshold));
                        otherwise % square
                            cfg.threshold_dB = 10*log10(abs(cfg.threshold));
                    end
            end
        end

        function n = get.noisePowersPerSample_W(obj)
            n = 10.^(.1*[obj.network.activeReceivingNodes.noisePowerPerSample_dB]); % kTBN
        end

        function w = get.expectedAmplitudes(obj)
            % (Ntx x Nrx x Nt matrix)
            switch obj.network.networkMode
                case "multiStatic"
                    w = 10.^(.05*obj.interfaces.receivedPowerFromScatterers_dBW); % (Ntx x Nrx x Nt matrix)
                case "monoStatic"
                    w = zeros(1, obj.network.numberOfActiveReceivingNodes, obj.interfaces.numberOfTargets, obj.monteCarlo.numberOfTrialsParallel);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        w(:, rxID, :, :) = 10.^(.05*obj.interfaces.receivedPowerFromScatterers_dBW(obj.network.monoStaticTransmitterIDs(rxID), rxID, :, :));
                    end
            end
        end
        
        function w = get.expectedAmplitudesReferenced(obj)
            % (Ntx x Nrx x Nt matrix)
            w = obj.expectedAmplitudes;
            w = w./w(:, obj.configuration.referenceActiveReceiverID, :);
            w(isnan(w)) = 0;
        end

        function w = get.integrationWeights(obj)
            % (1 x Nrx cell of Ntx x 1 x Nt x Nch matrix)
            % (Ntx x Nrx x Nt matrix)
            switch obj.network.networkMode
                case "multiStatic"
                    switch obj.processingAlgorithm
                        case 1 % Coherent detector - Perfect knowledge of range
                            %   Maximal ratio combining
                            %   Deterministic signal
        
                            w = obj.expectedAmplitudes;
                            w = w./sqrt(obj.noisePowersPerSample_W);
                            w(any(w == 0, 3), :) = 1./sqrt(obj.noisePowersPerSample_W(any(w == 0, 3)));
        
                        case 2 % Coherent detector - Fully correlated spatial observations
                            %   Spatially fully correlated phase and amplitude
                            %   For weak stochastic signal, it is optimal for any distribution
                            %   Otherwise, it is optimal under AWGN
                            %   For generalized LRT, it is optimal adaptive processing algorithm
        
                            w = 1;%obj.integrationWeights;
                            w = w./w(:, obj.configuration.referenceActiveReceiverID, :);
        
                            % w = obj.hypothesizedAmplitudes.*exp(1j*obj.hypothesizedPhasesReferenced)./obj.noisePowersPerSample_W;
                            % w(w == 0) = 1./obj.noisePowersPerSample_W;
        
                        case 3 % Noncoherent detector - Square-law envelope detector
                            %   Weak stochastic signal
                            %   Spatially independent phase, arbitrarily correlated amplitude
                            %   For weak stochastic signal, it is optimal for any distribution
                            %   Otherwise, it is an approximation of optimal processing algorithm under AWGN
        
                            w = obj.expectedAmplitudesReferenced.^2./obj.noisePowersPerSample_W;
                            w(w == 0) = 1./obj.noisePowersPerSample_W;
        
                        case 4 % Noncoherent detector - Linear-law envelope detector
                            %   Strong stochastic signal
                            %   Spatially independent phase, fully correlated amplitude
                            %   For generalized LRT, it is an approximation of optimal adaptive processing algorithm
        
                            w = obj.expectedAmplitudesReferenced./sqrt(obj.noisePowersPerSample_W);
                            w(w == 0) = 1./sqrt(obj.noisePowersPerSample_W);
        
                        case 5 % Noncoherent detector - Square-law envelope detector
                            %   Neither strong nor weak stochastic signal
                            %   Spatially independent phase and amplitude
        
                            w = obj.expectedAmplitudesReferenced.^2./obj.noisePowersPerSample_W/(1 + obj.inputSNR_lin);
                            w(w == 0) = 1./obj.noisePowersPerSample_W;
        
                        case 6 % Noncoherent detector - Square-law envelope detector
                            %   Strong stochastic signal
                            %   Spatially independent phase and amplitude
                            %   For generalized LRT, it is optimal adaptive processing algorithm
        
                            w = ones(obj.network.numberOfActiveTransmittingNodes, obj.network.numberOfActiveReceivingNodes);
                    end
                    if ~obj.network.passiveNetwork % Single transmitter networks
                        switch obj.network.networkCoherency
                            case "short-term coherent"
                                if any(obj.processingAlgorithm, [1, 2])
                                    error('Processing algorithm %d cannot be used for short-term coherent networks', obj.processingAlgorithm);
                                end
                            case "incoherent"
                                if any(obj.processingAlgorithm, [1, 2, 4])
                                    error('Processing algorithm %d cannot be used for incoherent networks', obj.processingAlgorithm);
                                end
                        end
                    else % Passive networks
                        error('Passive networks are not implemented');
                    end
                case "monoStatic"
                    w = ones(1, obj.network.numberOfActiveReceivingNodes);
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

        function A = get.hypothesizedAmplitudes(obj)
            gridScan = obj.gridPointsMesh;
            hypothesizedPositions = permute(cat(4, gridScan.x, gridScan.y, gridScan.z), [4 1 2 3]);
            hypothesizedPositions = reshape(hypothesizedPositions, [3, 1, 1, prod(obj.gridSize)]);
            rangeReceivers = sqrt(sum((hypothesizedPositions - permute([obj.network.activeReceivingNodes.position], [1 3 2])).^2));
            rangeTransmitters = sqrt(sum((hypothesizedPositions - [obj.network.activeTransmittingNodes.position]).^2));
            powerGains = 10*log10([obj.network.activeTransmittingNodes.inputPower_W].') + permute(obj.interfaces.targets.RCS_dbms, [1 3 2]) ...
            + 20*log10([obj.network.activeTransmittingNodes.carrierWavelength].');
            powerLosses = [obj.network.activeReceivingNodes.systemLoss_dB] + 30*log10(4*pi) + 20*log10(rangeReceivers.*rangeTransmitters);
            Pr = powerGains - powerLosses;
            A = permute(10.^(.05*Pr), [2 3 4 1]);
        end

        function A = get.hypothesizedAmplitudesReferenced(obj)
            A = obj.hypothesizedAmplitudes;
            A = A./A(:, obj.configuration.referenceActiveReceiverID, :, :);
            A(isnan(A)) = 0;
        end

        function phi = get.hypothesizedPhases(obj)
            phi = obj.hypothesizedTimeDelays./2./pi./[obj.network.activeTransmittingNodes.carrierFrequency];
        end

        function phi = get.hypothesizedPhasesReferenced(obj)
            phi = obj.hypothesizedPhases;
            phi = phi - phi(:, obj.configuration.referenceActiveReceiverID, :);
        end

        function y = get.signalsMatchFilteredTrials(obj)
            % (1 x Nrx cell of Ns + L - 1 x Nch x Ncmp matrix) Ncmp : number of parallel trials
            mc = obj.monteCarlo;
            t = [obj.network.activeReceivingNodes.samplingInstants]; % Ns x Nrx matrix
            demodulator = exp(1j*2*pi*shiftdim([obj.network.activeTransmittingNodes.carrierFrequency], -1).*t);
            mf = obj.network.matchFilter; % L x Nrx x Ntx matrix
            y = cell(1, obj.network.numberOfActiveReceivingNodes);
            signalsBeamformed = obj.interfaces.signalBeamformed;
            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                switch obj.network.beamformingMode
                    case 'bypass'
                        numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).array.numberOfTotalElements;
                    otherwise
                        numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
                end
                y{rxID} = zeros(size(t, 1) + size(mf, 1) - 1, numberOfTotalChannels, mc.numberOfTrialsParallel);
                s = signalsBeamformed{rxID}; % Ns x Nch x Nmcp matrix
                for chID = 1 : numberOfTotalChannels
                    for mcID = 1 : mc.numberOfTrialsParallel
                        signalChannel = s(:, chID, mcID);
                        signalChannel(isinf(signalChannel)) = 0;
                        y{rxID}(:, chID, mcID) = conv(signalChannel.*demodulator(:, rxID), conj((mf(:, rxID, 1))));
                    end
                end
                y{rxID} = y{rxID}./sqrt(obj.noisePowersPerSample_W(rxID));
            end
        end

        function Y = get.signalsIntegrated(obj)
            % (Ni x Nmcp matrix)
            if isscalar(obj.signalsMatchFiltered) || isempty(obj.signalsMatchFiltered)
                obj.setmatchfilteredsignals;
            end
            mc = obj.monteCarlo;
            visibleZone = ~obj.blindZone;
            switch obj.network.networkMode
                case "multiStatic"
                    Y = zeros(prod(obj.gridSize), 1, mc.numberOfTrialsParallel);
                    if isempty(obj.integrationIndices)
                        fprintf('integration indices are not set');
                        return;
                    end
                    switch obj.processingAlgorithm
                        case 1 % Coherent detector, deterministic signal
                            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                                Y(visibleZone, :, :) = Y(visibleZone, :, :) + obj.integrationWeights(:, rxID, 1).*obj.signalsMatchFiltered{rxID}(obj.integrationIndices(1, rxID, visibleZone), ceil(end/2), :);
                            end
                            Y = real(Y);
                        case 2 % Coherent detector, fully correlated stochastic signal
                            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                                Y(visibleZone, :, :) = Y(visibleZone, :, :) + obj.integrationWeights(:, rxID, 1).*obj.signalsMatchFiltered{rxID}(obj.integrationIndices(1, rxID, visibleZone), ceil(end/2), :);
                            end
                            Y = abs(Y).^2;
                        case {3, 5, 6} % Square-law envelope detectors
                            for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                                for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                                    %%% which channel not implemented
                                    Y(visibleZone, :, :) = Y(visibleZone, :, :) + obj.integrationWeights(txID, rxID, :).*abs(obj.signalsMatchFiltered{rxID}(obj.integrationIndices(txID, rxID, visibleZone), ceil(end/2), :)).^2; 
                                end
                            end
                        case 4 % Linear-law envelope detector
                            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                                Y(visibleZone, :, :) = Y(visibleZone, :, :) + obj.integrationWeights(:, rxID, :).*abs(obj.signalsMatchFiltered{rxID}(obj.integrationIndices(1, rxID, visibleZone), ceil(end/2), :));
                            end
                    end
                    if ~obj.network.passiveNetwork % Single transmitter networks
                    else % Passive networks
                    end
                case "monoStatic"
                    Y = cell(1, obj.network.numberOfActiveReceivingNodes);
                    if isempty(obj.integrationIndices)
                        fprintf('integration indices are not set');
                        return;
                    end
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        switch obj.network.beamformingMode
                            case 'bypass'
                                numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).array.numberOfTotalElements;
                            otherwise
                                numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
                        end
                        monoStaticTXIDs = obj.network.monoStaticTransmitterIDs(rxID);
                        y = zeros(prod(obj.gridSize), numberOfTotalChannels, mc.numberOfTrialsParallel);
                        y(visibleZone, :, :) = obj.integrationWeights(rxID).*obj.signalsMatchFiltered{rxID}(obj.integrationIndices(monoStaticTXIDs, rxID, visibleZone), :, :);
                        switch obj.processingAlgorithm
                            case 1 % Coherent detector, deterministic signal
                                Y{rxID} = real(y);
                            otherwise
                                Y{rxID} = abs(y).^2;
                        end
                    end
            end
        end

        function Y = get.signalsIntegratedDirectly(obj)
            if isscalar(obj.signalsMatchFiltered) || isempty(obj.signalsMatchFiltered)
                obj.setmatchfilteredsignals;
            end
            mc = obj.monteCarlo;
            switch obj.network.networkMode
                case "multiStatic"
                    sizeMF = [obj.network.activeReceivingNodes.numberOfSamplesPerCPI] + obj.network.pulseWidthSample(:, 1).' - 1;
                    Y = zeros([sizeMF, mc.numberOfTrialsParallel]);
                    switch obj.processingAlgorithm
                        case 1 % Coherent detector, deterministic signal
                        case 2 % Coherent detector, fully correlated stochastic signal
                        case {3, 5, 6} % Square-law envelope detectors
                            for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                                for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                                    %%% which channel not implemented
                                    Y = Y + shiftdim(abs(obj.signalsMatchFiltered{rxID}(:, ceil(end/2), :)).^2, 1 - rxID); 
                                end
                            end
                        case 4 % Linear-law envelope detector
                    end
                case "monoStatic"
                    Y = nan;
            end
        end

        function report = get.signalsCompressedTrials(obj)
            % (Ni x Nmcp matrix)
            mc = obj.monteCarlo;
            config = obj.configuration;
            t = [obj.network.activeReceivingNodes.samplingInstants]; % Ns x Nrx matrix
            demodulator = exp(1j*2*pi*shiftdim([obj.network.activeTransmittingNodes.carrierFrequency], -1).*t);
            Ns = size([obj.network.activeReceivingNodes.samplingInstants], 1);
            signalBeamformed = obj.interfaces.signalBeamformed;
            switch obj.network.networkMode
                case 'multiStatic'
                    Nrx = obj.network.numberOfActiveReceivingNodes;
                    measurements = zeros(Ns, Nrx, mc.numberOfTrialsParallel);
                    for rxID = 1 : Nrx
                        s = signalBeamformed{rxID}; % Ns x Nch x Nmcp matrix
                        measurements(:, rxID, :) = s(:, ceil(end/2), :).*demodulator(:, rxID)./sqrt(obj.noisePowersPerSample_W(rxID));
                    end
                case 'monoStatic'
                    Nrx = obj.network.activeReceivingNodes.array.numberOfTotalElements;
                    measurements = zeros(Ns, Nrx, mc.numberOfTrialsParallel);
                    s = signalBeamformed{1}; % Ns x Nch x Nmcp matrix
                    for rxID = 1 : Nrx
                        measurements(:, rxID, :) = s(:, rxID, :).*demodulator(:, 1)./sqrt(obj.noisePowersPerSample_W(1));
                    end
            end
            dict = obj.dictionaryCompression;
            Nc = size(dict, 2); % number of cells
            report = struct( ...
                'cellIDs', [], ...
                'numberOfIterations', [], ...
                'residualHistory', [], ...
                'integratedSignals', [], ...
                'complexAmplitudes', []);
            if isnan(config.numberOfTargets)
                targetResidual = config.threshold;
                numberOfIterations = Ns;
            else
                targetResidual = -Inf;
                numberOfIterations = config.numberOfTargets;
                if numberOfIterations > Ns
                    error('Number of targets cannot be larger than the number of samples');
                end
            end
            for mcID = 1 : mc.numberOfTrialsParallel
                switch obj.detectionAlgorithm
                    case "NCOMP"
                        orthAtomSet = zeros(Ns, numberOfIterations, Nrx);
                        nonorthAtomSet = zeros(Ns, numberOfIterations, Nrx);
                        numOfParallelOpts = Nrx;
                    case "OMP"
                        dict = reshape(permute(dict, [1 3 2]), Ns*Nrx, []);
                        measurements = reshape(measurements, Ns*Nrx, []);
                        orthAtomSet = zeros(Ns*Nrx, numberOfIterations);
                        nonorthAtomSet = zeros(Ns*Nrx, numberOfIterations);
                        numOfParallelOpts = 1;
                    case "CoSaMP"
                        %%% not implemented
                        dict = reshape(permute(dict, [1 3 2]), Ns*Nrx, []);
                        report(mcID) = CoSaMP(dict, measurements(:), config.numberOfTargets);
                end
                cellIDset = zeros(numberOfIterations, 1);
                residualPowerHistory = zeros(numberOfIterations + 1, 1);
                if obj.saveResiduals
                    report(mcID).integratedSignals = nan(Nc, numberOfIterations + 1);
                end
                thresholdReached = false;
                residual = measurements(:, :, mcID);
                for currentIterationID = 1 : numberOfIterations
                    switch obj.network.networkMode
                        case 'multiStatic'
                            switch obj.processingAlgorithm
                                case 6
                                    integratedSignal = sum(abs(pagemtimes(dict, 'ctranspose', permute(residual, [1 3 2]), 'none')).^2, 3); % noncoherent integration
                                    [maxPower, cellID] = max(integratedSignal);
                            end
                        case 'monoStatic'
                            switch obj.processingAlgorithm
                                case 2
                                    switch obj.network.beamformingMode
                                        case 'bypass'
                                            integratedSignal = sum(pagemtimes(dict, 'ctranspose', permute(residual, [1 3 2]), 'none'), 3)./sqrt(obj.network.activeReceivingNodes.array.numberOfTotalElements); % coherent integration
                                        otherwise
                                            integratedSignal = sum(pagemtimes(dict, 'ctranspose', permute(residual, [1 3 2]), 'none'), 3); % coherent integration
                                    end
                                    [maxPower, cellID] = max(abs(integratedSignal).^2);
                            end
                    end
                    % if currentIterationID == 1 && mcID == 1
                    %     fprintf('Max power = %g dB, Threshold = %g dB\n', 10*log10(maxPower), 10*log10(config.threshold));
                    % end
                    if maxPower < targetResidual
                        thresholdReached = true;
                        break;
                    end
                    if obj.saveResiduals
                        report(mcID).integratedSignals(:, currentIterationID) = integratedSignal;
                    end
                    residualPowerHistory(currentIterationID) = sum(abs(residual).^2, 'all')/Ns;
                    % if residualPowerHistory(currentIterationID) < targetResidual || (currentIterationID > 1 && residualPowerHistory(currentIterationID) > residualPowerHistory(currentIterationID - 1))
                    %     thresholdReached = true;
                    %     break;
                    % end
                    if ismember(cellID, cellIDset(1 : currentIterationID - 1))
                        warning('Shouldn''t happen...');
                        thresholdReached = true;
                        break;
                    end
                    cellIDset(currentIterationID) = cellID;
                    currentAtoms = permute(dict(:, cellID, :), [1 3 2]);
                    atomNorms = sqrt(sum(abs(currentAtoms).^2));
                    nonZeroAtoms = atomNorms ~= 0;
                    currentAtoms(:, nonZeroAtoms) = currentAtoms(:, nonZeroAtoms)./atomNorms(nonZeroAtoms);
                    nonorthAtomSet(:, currentIterationID, :) = currentAtoms; 
            
                    % -- Step 2: update residual
                
                    % First, orthogonalize 'atoms' against all previous atoms
                    % We use Modified Gram Schmidt
                    for previousIterationID = 1 : (currentIterationID - 1)
                        for rxID = 1 : numOfParallelOpts
                            currentAtoms(:, rxID) = currentAtoms(:, rxID) - (orthAtomSet(:, previousIterationID, rxID)'*currentAtoms(:, rxID))*orthAtomSet(:, previousIterationID, rxID);
                        end
                    end
                    % Second, normalize:
                    atomNorms = sqrt(sum(abs(currentAtoms).^2));
                    nonZeroAtoms = atomNorms ~= 0;
                    currentAtoms(:, nonZeroAtoms) = currentAtoms(:, nonZeroAtoms)./atomNorms(nonZeroAtoms);
                    orthAtomSet(:, currentIterationID, :) = currentAtoms;
                    % Third, solve least-squares problem
                    for rxID = 1 : numOfParallelOpts
                        % Fourth, update residual:
                        residual(:, rxID) = residual(:, rxID) - orthAtomSet(:, 1 : currentIterationID, rxID)*(orthAtomSet(:, 1 : currentIterationID, rxID)'*residual(:, rxID));
                    end
                end
                if thresholdReached
                    currentIterationID = currentIterationID - 1;
                end
                %%% iterations end
                % For the last iteration, we need to do this without orthogonalizing dictionary
                % so that the complexAmplitudes match what is expected.
                report(mcID).complexAmplitudes = zeros(currentIterationID, numOfParallelOpts);
                for rxID = 1 : numOfParallelOpts
                    if all(any(nonorthAtomSet(:, 1 : currentIterationID, rxID)))
                        report(mcID).complexAmplitudes(:, rxID) = pinv(nonorthAtomSet(:, 1 : currentIterationID, rxID))*measurements(:, rxID, mcID);
                    else
                        report(mcID).complexAmplitudes(:, rxID) = 0;
                    end
                end
                visibleCellIDs = false(Nc, 1);
                visibleCellIDs(cellIDset(1 : currentIterationID)) = true;
                gridCellIDs = false(prod(obj.gridSize), 1);
                gridCellIDs(~obj.blindZone) = visibleCellIDs;
                report(mcID).cellIDs = find(gridCellIDs);
                report(mcID).residualHistory = residualPowerHistory(1 : currentIterationID + 1);
                if obj.saveResiduals
                    report(mcID).integratedSignals = report(mcID).integratedSignals(:, 1 : currentIterationID + 1) ;
                end
                report(mcID).numberOfIterations = currentIterationID + 1;
            end
        end

        function T = get.thresholdCFAR(obj)
            % (1 x Nrx cell of Ns + L - 1 x Nch x Nmcp matrix)
            T = cell(1, obj.network.numberOfActiveReceivingNodes);
            config = obj.configuration;
            mc = obj.monteCarlo;
            NTC = config.numberOfTrainingCells;
            NGC = config.numberOfGuardCells;
            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                switch obj.network.beamformingMode
                    case 'bypass'
                        numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).array.numberOfTotalElements;
                    otherwise
                        numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
                end
                txID = obj.network.monoStaticTransmitterIDs(rxID);
                L = obj.network.pulseWidthSample(rxID, txID) - 2;
                N = obj.network.activeReceivingNodes(rxID).numberOfSamplesPerCPI;
                samples = -L : N;
                numberOfSamples = length(samples);
                T{rxID} = zeros(numberOfSamples, numberOfTotalChannels, mc.numberOfTrialsParallel);
                for sampleID = 1 : numberOfSamples % OSCFAR
                    trainingStart = sampleID + NGC + 1;
                    upperTraining = min(trainingStart : trainingStart + NTC - 1, numberOfSamples);
                    trainingStart = sampleID - NGC - 1;
                    lowerTraining = max(trainingStart - NTC + 1 : trainingStart, 1);
                    sampleIdx = [lowerTraining, upperTraining];
                    Nt = length(sampleIdx);
                    alpha = Nt.*(config.PFA.^(-1/Nt) - 1);
                    T{rxID}(sampleID, :, :) = alpha.*mean(abs(obj.signalsMatchFiltered{rxID}(sampleIdx, :, :)).^2);
                end
            end
        end

        function detection = get.detectionFromMatchFiltration(obj)
            % (Nrx x Nmcp cell of 1 x Nd vector)
            if isscalar(obj.signalsMatchFiltered) || isempty(obj.signalsMatchFiltered)
                obj.setmatchfilteredsignals;
            end
            config = obj.configuration;
            mc = obj.monteCarlo;
            detection = struct( ...
                'power', [], ...
                'range', [], ...
                'elevation', [], ...
                'azimuth', [], ...
                'numberOfDetections', [], ...
                'position', []);
            detection = repmat(detection, obj.network.numberOfActiveReceivingNodes, mc.numberOfTrialsParallel);
            switch obj.network.networkMode
                case "monoStatic"
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        Ts = obj.network.activeReceivingNodes(rxID).samplingPeriod;
                        N = obj.network.activeReceivingNodes(rxID).numberOfSamplesPerCPI;
                        switch obj.network.beamformingMode
                            case 'bypass'
                                numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).array.numberOfTotalElements;
                            otherwise
                                numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
                        end
                        if isempty(obj.network.receivingNodes(rxID).directionFinder.table)
                            obj.network.setdirectionfinders;
                        end
                        DF = obj.network.receivingNodes(rxID).directionFinder.table;
                        azimuthDF = obj.network.receivingNodes(rxID).directionFinder.azimuth;
                        azimuthResolution = obj.network.receivingNodes(rxID).directionFinder.azimuthResolution;
                        elevationSteering = obj.network.receivingNodes(rxID).array.steeringElevation;
                        azimuthSteering = obj.network.receivingNodes(rxID).array.steeringAzimuth;
                        sizeDF = size(DF, 1);
                        rotationMat = obj.network.receivingNodes(rxID).array.rotationMatrix;
                        origin = obj.network.receivingNodes(rxID).position;
                        txID = obj.network.monoStaticTransmitterIDs(rxID);
                        L = obj.network.pulseWidthSample(rxID, txID) - 2;
                        Ns = N + L + 1;
                        T = config.threshold(rxID)*ones(Ns*numberOfTotalChannels, mc.numberOfTrialsParallel);
                        obs = obj.signalsMatchFiltered{rxID};
                        absObs = reshape(abs(obs), [], mc.numberOfTrialsParallel);
                        switch obj.processingAlgorithm
                            case 1 % Coherent detector, deterministic signal
                                thObs = real(obj.signalsMatchFiltered{rxID});
                            otherwise % Square-law envelope detectors
                                if config.CFAR
                                    T = reshape(obj.thresholdCFAR{rxID}, [], mc.numberOfTrialsParallel);
                                end
                                thObs = abs(obj.signalsMatchFiltered{rxID}).^2;
                        end
                        thObs = reshape(thObs, [], mc.numberOfTrialsParallel);
                        for mcID = 1 : mc.numberOfTrialsParallel
                            switch obj.network.activeTransmittingNodes(txID).transmissionType
                                case "continuous"
                                    error('not implemented');
                                case "pulsed"
                                    switch obj.detectionAlgorithm
                                        case "thresholding"
                                            detectionIndex = find(thObs(:, mcID) > T(:, mcID));
                                        case "peak"
                                            detectionIndex = find(thObs(:, mcID) > T(:, mcID));
                                            [~, peakIndex] = max(thObs(detectionIndex, mcID));
                                            detectionIndex = detectionIndex(peakIndex);
                                        otherwise
                                            detectionIndex = [];
                                    end
                            end
                            detectionSamples = mod(detectionIndex - 1, Ns) + 1;
                            if config.removeBoundaryDetectionMF
                                boundaryDetections = detectionSamples ~= 1 & detectionSamples ~= Ns;
                                detectionIndex = detectionIndex(boundaryDetections);
                                detectionSamples = detectionSamples(boundaryDetections);
                            end
                            numberOfDetections = length(detectionIndex);
                            if numberOfDetections
                                detectionPower = 20*log10(absObs(detectionIndex, mcID));
                                idx = 1;
                                while idx < numberOfDetections
                                    group = abs(detectionSamples - detectionSamples(idx)) < config.groupingCloseness;
                                    groupPowers = detectionPower(group);
                                    groupSamples = detectionSamples(group);
                                    groupIndices = detectionIndex(group);
                                    [powerRepresenter, groupRepresenterIndex] = max(groupPowers);
                                    indexRepresenter = groupIndices(groupRepresenterIndex);
                                    sampleRepresenter = groupSamples(groupRepresenterIndex);
                                    detectionPower = [powerRepresenter; detectionPower(~group)];
                                    detectionIndex = [indexRepresenter; detectionIndex(~group)];
                                    detectionSamples = [sampleRepresenter; detectionSamples(~group)];
                                    idx = idx + 1 - sum(group(idx - 1 : -1 : 1));
                                    numberOfDetections = length(detectionIndex);
                                end
                                if config.interpolationMF
                                    detectionSamplesInterpolated = zeros(numberOfDetections, 1);
                                    interpolatables = detectionSamples ~= 1 & detectionSamples ~= Ns;
                                    rightNgbh = absObs(detectionIndex(interpolatables) + 1, mcID);
                                    leftNgbh = absObs(detectionIndex(interpolatables) - 1, mcID);
                                    center = absObs(detectionIndex(interpolatables), mcID);
                                    detectionSamplesInterpolated(interpolatables) = detectionSamples(interpolatables) + (leftNgbh - rightNgbh)./(2*center - (rightNgbh + leftNgbh))/2;
                                else
                                    detectionSamplesInterpolated = detectionSamples;
                                end
                                spectrum = abs(DF*obs(detectionSamples, :, mcID)');
                                [maxVal, idxDF] = max(spectrum);
                                if config.removeBoundaryDetectionDF
                                    boundaryDetections = idxDF == 1 | idxDF == length(azimuthDF);
                                    maxVal = maxVal(~boundaryDetections);
                                    idxDF = idxDF(~boundaryDetections);
                                    detectionPower = detectionPower(~boundaryDetections);
                                    detectionSamplesInterpolated = detectionSamplesInterpolated(~boundaryDetections);
                                    numberOfDetections = length(idxDF);
                                end
                                if numberOfDetections
                                    elevation = zeros(1, length(idxDF));
                                    azimuth = azimuthDF(idxDF);
                                    if config.interpolationDF
                                        interpolatables = idxDF ~= 1 & idxDF ~= length(azimuthDF);
                                        specVec = sizeDF*(0 : numberOfDetections - 1);
                                        rightNgbh = spectrum(idxDF(interpolatables) + 1 + specVec(interpolatables));
                                        leftNgbh = spectrum(idxDF(interpolatables) - 1 + specVec(interpolatables));
                                        corrections = (leftNgbh - rightNgbh)./(2*maxVal(interpolatables) - (rightNgbh + leftNgbh))/2;
                                        azimuth(interpolatables) = azimuth(interpolatables) + azimuthResolution.*corrections;
                                    end
                                    range = obj.speedOfLight.*Ts*(detectionSamplesInterpolated.' - L - 1)/2;
                                    elevation = asind(sind(elevation) + sind(elevationSteering));
                                    azimuth = asind(sind(azimuth) + sind(azimuthSteering));
                                    imaginaryDetections = imag(azimuth) | imag(elevation) | imag(range);
                                    range = range(~imaginaryDetections);
                                    elevation = elevation(~imaginaryDetections);
                                    azimuth = azimuth(~imaginaryDetections);
                                    numberOfDetections = length(range);
                                    if numberOfDetections
                                        [x, y, z] = sph2cart(pi*azimuth/180, pi*elevation/180, range);
                                        position = [x; y; z];
                                        detection(rxID, mcID).power = detectionPower;
                                        detection(rxID, mcID).range = range;
                                        detection(rxID, mcID).elevation = elevation;
                                        detection(rxID, mcID).azimuth = azimuth;
                                        detection(rxID, mcID).numberOfDetections = numberOfDetections;
                                        detection(rxID, mcID).position = rotationMat*position + origin;
                                    end
                                end
                            end
                        end
                    end
                case "multiStatic"
                    warning('network mode must be monostatic');
            end
        end

        function idx = get.hypothesisTestingResults(obj)
            switch obj.detectionAlgorithm
                case {"thresholding", "peak"}
                    if isscalar(obj.signalsMatchFiltered) || isempty(obj.signalsMatchFiltered)
                        obj.setmatchfilteredsignals;
                    end
                    z = obj.signalsIntegrated;
                case {"NCOMP", "OMP"}
                    if isempty(obj.compressionReport(1).complexAmplitudes)
                        obj.applycompression;
                    end
                case "CoSaMP"
            end
            config = obj.configuration;
            mc = obj.monteCarlo;
            switch obj.network.networkMode
                case "multiStatic"
                    idx = cell(1, mc.numberOfTrialsParallel);
                    for mcID = 1 : mc.numberOfTrialsParallel
                        switch obj.detectionAlgorithm
                            case "thresholding"
                                idx{mcID} = find(z(:, :, mcID) > config.threshold);
                            case "peak"
                                idx{mcID} = find(z(:, :, mcID) > config.threshold);
                                [~, peakIndex] = max(z(idx{mcID}, :, mcID));
                                idx{mcID} = idx{mcID}(peakIndex);
                            case {"NCOMP", "OMP"}
                                idx{mcID} = obj.compressionReport(mcID).cellIDs;
                            case "CoSaMP"
                        end
                    end
                    if ~obj.network.passiveNetwork % Single transmitter networks
                    else % Passive networks
                    end
                case "monoStatic"
                    idx = cell(1, obj.network.numberOfActiveReceivingNodes);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        T = config.threshold(rxID);
                        switch obj.network.beamformingMode
                            case 'bypass'
                                switch obj.detectionAlgorithm
                                    case {"thresholding", "peak"}
                                        numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).array.numberOfTotalElements;
                                    otherwise
                                        numberOfTotalChannels = 1;
                                end
                            otherwise
                                numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
                        end
                        idx{rxID} = cell(numberOfTotalChannels, mc.numberOfTrialsParallel);
                        for mcID = 1 : mc.numberOfTrialsParallel
                            switch obj.detectionAlgorithm
                                case "thresholding"
                                    for chID = 1 : numberOfTotalChannels
                                        idx{rxID}{chID, mcID} = find(z{rxID}(:, chID, mcID) > T);
                                    end
                                case "peak"
                                    for chID = 1 : numberOfTotalChannels
                                        idx{rxID}{chID, mcID} = find(z{rxID}(:, chID, mcID) > T);
                                        [~, peakIndex] = max(z{rxID}(idx{rxID}{chID, mcID}, chID, mcID));
                                        idx{rxID}{chID, mcID} = idx{rxID}{chID, mcID}(peakIndex);
                                    end
                                case {"NCOMP", "OMP"}
                                    idx{rxID}{1, mcID} = obj.compressionReport(mcID).cellIDs;
                                case "CoSaMP"
                            end
                        end
                    end
            end
        end

        function pos = get.estimatedPositions(obj)
            mc = obj.monteCarlo;
            gridScan = obj.gridPointsMesh;
            idx = obj.hypothesisTestingResults;
            switch obj.network.networkMode
                case "multiStatic"
                    pos = cell(1, mc.numberOfTrialsParallel);
                    for mcID = 1 : mc.numberOfTrialsParallel
                        pos{mcID} = [gridScan.x(idx{mcID}) gridScan.y(idx{mcID}) gridScan.z(idx{mcID})].';
                    end
                case "monoStatic"
                    pos = cell(1, obj.network.numberOfActiveReceivingNodes);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        switch obj.network.beamformingMode
                            case 'bypass'
                                numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).array.numberOfTotalElements;
                            otherwise
                                numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
                        end
                        pos{rxID} = cell(numberOfTotalChannels, mc.numberOfTrialsParallel);
                        for mcID = 1 : mc.numberOfTrialsParallel
                            for chID = 1 : numberOfTotalChannels
                                pos{rxID}{chID, mcID} = [gridScan.x(idx{rxID}{chID, mcID}) gridScan.y(idx{rxID}{chID, mcID}) gridScan.z(idx{rxID}{chID, mcID})].';
                            end
                        end
                    end
            end
        end

        function detection = get.detectionFromIntegration(obj)
            % (Nrx x Nmcp cell of 1 x Nd vector)
            if isscalar(obj.signalsMatchFiltered) || isempty(obj.signalsMatchFiltered)
                obj.setmatchfilteredsignals;
            end
            config = obj.configuration;
            mc = obj.monteCarlo;
            gridScan = obj.gridPointsMesh;
            detection = struct( ...
                'power', [], ...
                'x', [], ...
                'y', [], ...
                'z', [], ...
                'numberOfDetections', [], ...
                'position', []);
            detection = repmat(detection, 1, mc.numberOfTrialsParallel);
            if isempty(obj.integrationIndices)
                fprintf('integration indices are not set');
                return;
            end
            Y = zeros(prod(obj.gridSize), 1, mc.numberOfTrialsParallel);
            switch obj.processingAlgorithm
                case 1 % Coherent detector, deterministic signal
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        Y = Y + obj.integrationWeights(:, rxID, :).*obj.signalsMatchFiltered{rxID}(obj.integrationIndices(1, rxID, :), 1, :);
                    end
                    Y = real(Y);
                case 2 % Coherent detector, fully correlated stochastic signal
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        Y = Y + obj.integrationWeights(:, rxID, :).*obj.signalsMatchFiltered{rxID}(obj.integrationIndices(1, rxID, :), 1, :);
                    end
                    Y = abs(Y).^2;
                case {3, 5, 6} % Square-law envelope detectors
                    for txID = 1 %: obj.network.numberOfActiveTransmittingNodes
                        for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                            %%% which channel not implemented
                            Y = Y + obj.integrationWeights(txID, rxID, :).*abs(obj.signalsMatchFiltered{rxID}(obj.integrationIndices(txID, rxID, :), ceil(end/2), :)).^2;
                        end
                    end
                case 4 % Linear-law envelope detector
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        Y = Y + obj.integrationWeights(:, rxID, :).*abs(obj.signalsMatchFiltered{rxID}(obj.integrationIndices(1, rxID, :), 1, :));
                    end
            end
            if ~obj.network.passiveNetwork % Single transmitter networks
            else % Passive networks
            end

            for mcID = 1 : mc.numberOfTrialsParallel
                switch obj.detectionAlgorithm
                    case "thresholding"
                        idx = find(Y(:, :, mcID) > config.threshold);
                    case "peak"
                        idx = find(Y(:, :, mcID) > config.threshold);
                        [~, peakIndex] = max(Y(idx, :, mcID));
                        idx = idx(peakIndex);
                    otherwise
                        idx = [];
                end
                x = gridScan.x(idx);
                y = gridScan.y(idx);
                z = gridScan.z(idx);
                position = [x y z].';
                detection(mcID).power = Y(idx, :, mcID);
                detection(mcID).x = x;
                detection(mcID).y = y;
                detection(mcID).z = z;
                detection(mcID).numberOfDetections = numel(idx);
                detection(mcID).position = position;
            end
        end

        %%%

        function err = get.positionError(obj)
            mc = obj.monteCarlo;
            switch obj.network.networkMode
                case "multiStatic"
                    err = cell(1, mc.numberOfTrialsParallel);
                    for mcID = 1 : mc.numberOfTrialsParallel
                        idx = obj.estimatedPositions{mcID};
                        numberOfEstimations = size(idx, 2);
                        err{mcID} = zeros(3, numberOfEstimations, obj.interfaces.numberOfTargets);
                        if ~isempty(obj.interfaces.targets)
                            for targetID = 1 : obj.interfaces.numberOfTargets
                                err{mcID}(:, :, targetID) = obj.interfaces.targets.position(:, targetID) - idx;
                            end
                        end
                    end
                case "monoStatic"
                    err = cell(1, obj.network.numberOfActiveReceivingNodes);
                    gridScan = obj.gridPointsMesh;
                    idx = obj.hypothesisTestingResults;
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        switch obj.network.beamformingMode
                            case 'bypass'
                                numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).array.numberOfTotalElements;
                            otherwise
                                numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
                        end
                        err{rxID} = cell(numberOfTotalChannels, mc.numberOfTrialsParallel);
                        for mcID = 1 : mc.numberOfTrialsParallel
                            for chID = 1 : numberOfTotalChannels
                                indices = idx{rxID}{chID, mcID};
                                estimatedPosition = [gridScan.x(indices) gridScan.y(indices) gridScan.z(indices)].';
                                if ~isempty(estimatedPosition)
                                    numberOfEstimations = size(estimatedPosition, 2);
                                    err{rxID}{chID, mcID} = zeros(3, numberOfEstimations, obj.interfaces.numberOfTargets);
                                    if ~isempty(obj.interfaces.targets)
                                        for targetID = 1 : obj.interfaces.numberOfTargets
                                            err{rxID}{chID, mcID}(:, :, targetID) = permute(obj.interfaces.targets.position(:, targetID), [1 3 2]) - estimatedPosition;
                                        end
                                    end
                                end
                            end
                        end
                    end
            end
        end

        function err = get.positionErrorTotal(obj)
            mc = obj.monteCarlo;
            posErr = obj.positionError;
            switch obj.network.networkMode
                case "multiStatic"
                    err = cell(1, mc.numberOfTrialsParallel);
                    for mcID = 1 : mc.numberOfTrialsParallel
                        err{mcID} = sqrt(sum(posErr{mcID}.^2));
                    end
                case "monoStatic"
                    err = cell(1, obj.network.numberOfActiveReceivingNodes);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        switch obj.network.beamformingMode
                            case 'bypass'
                                numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).array.numberOfTotalElements;
                            otherwise
                                numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
                        end
                        err{rxID} = cell(numberOfTotalChannels, mc.numberOfTrialsParallel);
                        for mcID = 1 : mc.numberOfTrialsParallel
                            for chID = 1 : numberOfTotalChannels
                                err{rxID}{chID, mcID} = shiftdim(sqrt(sum(posErr{rxID}{chID, mcID}.^2)), 1);
                            end
                        end
                    end
            end
        end

        function C = get.covarianceMatrixInitial(obj)
            switch obj.network.networkMode
                case "multiStatic"
                case "monoStatic"
                    lambda = [obj.network.transmittingNodes.carrierWavelength];
                    arraysRX = [obj.network.receivingNodes.array];
                    N = [arraysRX.numberOfElements];
                    d = [arraysRX.spacing];
                    azimuthBeamWidth = 70*lambda./(N(1, :).*d(1, :));
                    elevationBeamWidth = 70*lambda./(N(2, :).*d(2, :));
                    rangeResolution = obj.speedOfLight.*[obj.network.receivingNodes.samplingPeriod]/2;
                    C = zeros(3, 3, obj.network.numberOfActiveReceivingNodes);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        C(:, :, rxID) = diag([rangeResolution(rxID), elevationBeamWidth(rxID), azimuthBeamWidth(rxID)]); % R, theta, phi
                    end
            end
        end

        %%% set methods

        function setintegrationindices(obj)
            obj.integrationIndices = ones(obj.network.numberOfActiveTransmittingNodes, obj.network.numberOfActiveReceivingNodes, prod(obj.gridSize));
            obj.integrationWeights = ones(obj.network.numberOfActiveTransmittingNodes, obj.network.numberOfActiveReceivingNodes, prod(obj.gridSize));
            obj.blindZone = true(1, prod(obj.gridSize));
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
                    timeDifference = abs(obj.hypothesizedTimeDelays(txID, rxID, :) - timeDelays);
                    obj.blindZone = obj.blindZone & permute(all(timeDifference > Ts, 4), [1 3 2]);
                    [~, obj.integrationIndices(txID, rxID, :)] = min(timeDifference, [], 4);
                    % obj.integrationWeights(txID, rxID, :) = obj.hypothesizedAmplitudes(txID, rxID, :).*exp(1j*obj.hypothesizedPhases(txID, rxID, :));
                end
            end
            %%% add back of array to blind zone
            gridScan = obj.gridPointsMesh;
            hypothesizedPositions = reshape(permute(cat(4, gridScan.x, gridScan.y, gridScan.z), [4 1 2 3]), [3 prod(obj.gridSize)]);
            targets = target('position', hypothesizedPositions);
            int = interface('network', obj.network, 'targets', targets);
            backOfArray = all(int.transmitBackOfArrayTarget, 3) & all(int.receiveBackOfArrayTarget, 3);
            obj.blindZone(backOfArray) = true;
        end

        function setmatchfilteredsignals(obj)
            obj.signalsMatchFiltered = obj.signalsMatchFilteredTrials;
        end

        function setdictionary(obj)
            Ns = size([obj.network.activeReceivingNodes.samplingInstants], 1);
            gridScan = obj.gridPointsMesh;
            hypothesizedPositions = reshape(permute(cat(4, gridScan.x, gridScan.y, gridScan.z), [4, 1, 2, 3]), [3, prod(obj.gridSize)]);
            hypothesizedPositions = hypothesizedPositions(:, ~obj.blindZone);
            targets = target( ...
                "position", hypothesizedPositions, ...
                "meanRCS_dbms", 0);
            int = interface( ...
                'network', obj.network, ...
                'targets', targets);
            int.configure( ...
                'noise', 0, ...
                'directPath', 0, ...
                'pathLoss', 0, ...
                'spatialCoherency', 'deterministic');
            signals = int.signalReceivedFromScatterers; % Ns x 1 x Nt x Ntx x M
            switch obj.network.networkMode
                case 'multiStatic'
                    dictionary = zeros(Ns, obj.network.numberOfActiveReceivingNodes, size(hypothesizedPositions, 2));
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        dictionary(:, rxID, :) = signals{rxID};
                        % dictionary(:, rxID, :) = sum(signals{rxID}, 4); %%% not implemented
                    end
                case 'monoStatic'
                    dictionary = zeros(Ns, obj.network.activeReceivingNodes.array.numberOfTotalElements, size(hypothesizedPositions, 2));
                    for elementID = 1 : obj.network.activeReceivingNodes.array.numberOfTotalElements
                        dictionary(:, elementID, :) = signals{1}(:, :, :, :, elementID);
                    end
            end
            obj.dictionaryCompression = permute(dictionary, [1 3 2]);
        end

        function applycompression(obj, options)
            arguments
                obj
                options.saveResiduals (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.saveResiduals, [0, 1])} = 0
            end
            if isempty(obj.dictionaryCompression)
                obj.setdictionary;
            end
            obj.saveResiduals = options.saveResiduals;
            obj.compressionReport = obj.signalsCompressedTrials;
        end

        function configure(obj, options)
            arguments
                obj
                options.PFA (1, 1) double {mustBeNonnegative} = obj.configuration.PFA
                options.numberOfTargets (1, 1) double = obj.configuration.numberOfTargets
                options.removeBoundaryDetectionDF (1, 1) double {mustBeNonnegative} = obj.configuration.removeBoundaryDetectionDF
                options.removeBoundaryDetectionMF (1, 1) double {mustBeNonnegative} = obj.configuration.removeBoundaryDetectionMF
                options.interpolationDF (1, 1) double {mustBeNonnegative} = obj.configuration.interpolationDF
                options.interpolationMF (1, 1) double {mustBeNonnegative} = obj.configuration.interpolationMF
                options.groupingCloseness (1, 1) double {mustBePositive} = obj.configuration.groupingCloseness
                options.numberOfTrainingCells (1, 1) {mustBePositive, mustBeInteger} = obj.configuration.numberOfTrainingCells
                options.numberOfGuardCells (1, 1) {mustBeNonnegative, mustBeInteger} = obj.configuration.numberOfGuardCells
                options.CFAR (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.CFAR, [0, 1])} = obj.configuration.CFAR
                options.processingAlgorithm (1, 1) {mustBeInteger, mustBeInRange(options.processingAlgorithm, 1, 6)} = obj.processingAlgorithm
                options.detectionAlgorithm (1, 1) string {mustBeMember(options.detectionAlgorithm, ["thresholding", "peak", "CoSaMP", "OMP", "NCOMP"])} = obj.detectionAlgorithm
                options.numberOfTrials (1, 1) {mustBeNonnegative, mustBeInteger} = obj.monteCarlo.numberOfTrials
                options.numberOfTrialsParallel (1, 1) {mustBeNonnegative, mustBeInteger} = obj.monteCarlo.numberOfTrialsParallel
                options.seed (1, 1) {mustBeNonnegative, mustBeInteger, mustBeLessThan(options.seed, 4294967296)} = 0
                options.seedShuffle (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.seedShuffle, [0, 1])} = obj.seedShuffle
            end
            cfg = obj.configuration;
            cfg.PFA = options.PFA;
            cfg.numberOfTargets = options.numberOfTargets;
            cfg.removeBoundaryDetectionDF = options.removeBoundaryDetectionDF;
            cfg.removeBoundaryDetectionMF = options.removeBoundaryDetectionMF;
            cfg.interpolationDF = options.interpolationDF;
            cfg.interpolationMF = options.interpolationMF;
            cfg.groupingCloseness = options.groupingCloseness;
            cfg.numberOfTrainingCells = options.numberOfTrainingCells;
            cfg.numberOfGuardCells = options.numberOfGuardCells;
            cfg.CFAR = options.CFAR;
            obj.configuration = cfg;
            if any(obj.processingAlgorithm == [1 2]) && ~strcmpi(obj.network.networkCoherency, "coherent")
                warning('processing algorithm %d requires coherent network', options.processingAlgorithm);
            else
                obj.processingAlgorithm = options.processingAlgorithm;
            end
            obj.detectionAlgorithm = options.detectionAlgorithm;
            switch obj.detectionAlgorithm %%% not implemented fully
                case 'OMP'
                    switch obj.network.networkMode
                        case 'multiStatic'
                            obj.processingAlgorithm = 6;
                        case 'monoStatic'
                            obj.processingAlgorithm = 2;
                    end
                case 'NCOMP'
                    switch obj.network.networkMode
                        case 'multiStatic'
                            obj.processingAlgorithm = 6;
                        case 'monoStatic'
                            obj.processingAlgorithm = 2;
                    end
            end
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

        function cellID = neighbours(obj, cellID, options)
            arguments
                obj
                cellID (1, 1) {mustBePositive, mustBeInteger}
                options.offset (1, 1) {mustBePositive, mustBeInteger} = 1
            end
            [xInd, yInd, zInd] = ind2sub(obj.gridSize, cellID);
            xNeigbours = validneighbours(xInd, 1);
            yNeigbours = validneighbours(yInd, 2);
            zNeigbours = validneighbours(zInd, 3);
            [xNeigbours, yNeigbours, zNeigbours] = meshgrid(xNeigbours, yNeigbours, zNeigbours);
            cellID = sub2ind(obj.gridSize, xNeigbours(:), yNeigbours(:), zNeigbours(:));
            function ngbhrs = validneighbours(center, index)
                ngbhrs = [center, center + (1 : options.offset), center - (1 : options.offset)];
                ngbhrs = ngbhrs(ngbhrs > 0 & ngbhrs <= obj.gridSize(index));
            end
        end

        %%% monte carlo simulation

        function [PD, PFA] = simulatedetection(obj)
            arguments
                obj
            end
            mc = obj.monteCarlo;
            gridScan = obj.gridPointsMesh;
            dims = find(size(gridScan.x) ~= 1);
            dimensions = {"x", "y", "z"};
            dimensions = dimensions(dims);
            numberOfCells = prod(obj.gridSize);
            probabilities = zeros(obj.gridSize([2, 1, 3]));
            numberOfTotalTrials = mc.numberOfTrials*mc.numberOfTrialsParallel;
            for mcID = 1 : mc.numberOfTrials
                switch obj.detectionAlgorithm
                    case {"thresholding", "peak"}
                        obj.setmatchfilteredsignals;
                    case {"NCOMP", "OMP"}
                        obj.applycompression;
                    case "CoSaMP"
                end
                idx = obj.hypothesisTestingResults;
                for mcpID = 1 : mc.numberOfTrialsParallel
                    probabilities(idx{mcpID}) = probabilities(idx{mcpID}) + 1/numberOfTotalTrials;
                end
                if ~mod(mcID, 10)
                    fprintf('trial = %d\n', mcID);
                end
            end
            if ~isempty(obj.interfaces.targets)
                targetCellIDs = cell(1, obj.interfaces.numberOfTargets);
                cellPositions = cat(3, gridScan.x, gridScan.y, gridScan.z);
                cellPositions = reshape(cellPositions(:, :, dims), [], numel(dims)).';
                targetPositions = obj.interfaces.targets.position(dims, :);
                for targetID = 1 : obj.interfaces.numberOfTargets
                    targetCellIDs{targetID} = find(sum(abs(targetPositions(:, targetID) - cellPositions).^2) < sum(abs(obj.gridResolution(dims)).^2));
                end
                targetCellIDs = unique(cell2mat(targetCellIDs));
            else
                targetCellIDs = [];
            end
            noiseCellIDs = setdiff(1 : numberOfCells, targetCellIDs);
            PD = nan(obj.gridSize([2, 1, 3]));
            PFA = nan(obj.gridSize([2, 1, 3]));
            PD(targetCellIDs) = probabilities(targetCellIDs);
            PFA(noiseCellIDs) = probabilities(noiseCellIDs);
            xLabel = dimensions{1} + " (km)";
            yLabel = dimensions{2} + " (km)";
            x1 = obj.gridPoints{1}/1e3;
            x2 = obj.gridPoints{2}/1e3;
            figure; img = imagesc(x1, x2, PD);
            targetCells = false(obj.gridSize([2, 1, 3]));
            targetCells(targetCellIDs) = true;
            set(img, 'AlphaData', targetCells);
            delete(datatip(img, 2, 2));
            img.DataTipTemplate.DataTipRows(1).Label = "x";
            img.DataTipTemplate.DataTipRows(1).Value = gridScan.x;
            img.DataTipTemplate.DataTipRows(2).Label = "y";
            img.DataTipTemplate.DataTipRows(2).Value = gridScan.y;
            img.DataTipTemplate.DataTipRows(3).Label = "PD";
            img.DataTipTemplate.DataTipRows(3).Value = PD;
            colorbar; colormap('spring');
            hold on;
            if ~isempty(obj.interfaces.targets)
                plot3(targetPositions(1, :)/1e3, targetPositions(2, :)/1e3, ones(1, size(targetPositions, 2)), 'ok', 'LineWidth', 1);
            end
            grid on; grid minor; view(0, 90);
            xlabel(xLabel); ylabel(yLabel); zlabel('P_D');
            title('Probability of detection'); hold off;
            figure; img = imagesc(x1, x2, PFA);
            noiseCells = true(obj.gridSize([2, 1, 3]));
            noiseCells(targetCellIDs) = false;
            set(img, 'AlphaData', noiseCells);
            delete(datatip(img, 2, 2));
            img.DataTipTemplate.DataTipRows(1).Label = "x";
            img.DataTipTemplate.DataTipRows(1).Value = gridScan.x;
            img.DataTipTemplate.DataTipRows(2).Label = "y";
            img.DataTipTemplate.DataTipRows(2).Value = gridScan.y;
            img.DataTipTemplate.DataTipRows(3).Label = "PFA";
            img.DataTipTemplate.DataTipRows(3).Value = PFA;
            colorbar; colormap('spring');
            grid on; grid minor; view(0, 90);
            xlabel(xLabel); ylabel(yLabel); zlabel('P_{FA}');
            title(sprintf('Probability of false alarm %g', mean(PFA, 'all', 'omitnan'))); hold off;
        end

        function simulatecoverage(obj, options)
            arguments
                obj
                options.meanRCS_dbms (1, 1) double = 0
                options.onCellCenters (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.onCellCenters, [0, 1])} = 0; 
            end
            mc = obj.monteCarlo;
            numberOfTotalTrials = mc.numberOfTrials*mc.numberOfTrialsParallel;
            originalTargets = obj.interfaces.targets;
            originalNumberOfTargets = obj.configuration.numberOfTargets;
            cleanup = onCleanup(@() cleanupFunction(obj, originalTargets, originalNumberOfTargets));
            obj.configure("numberOfTargets", nan);
            gridScan = obj.gridPointsMesh;
            dimensions = find(size(gridScan.x) ~= 1);
            visibleZone = ~obj.blindZone;
            cellPositions = cat(3, gridScan.x, gridScan.y, gridScan.z);
            cellPositions = reshape(cellPositions(:, :, dimensions), [], numel(dimensions)).';
            cellPositions = cellPositions(:, visibleZone);
            width = zeros(3, 1);
            width(dimensions) = obj.gridResolution(dimensions);
            obj.coverageSimulationReport.PD = zeros(obj.gridSize([2, 1, 3]));
            obj.coverageSimulationReport.thresholds = zeros(obj.gridSize([2, 1, 3]));
            obj.coverageSimulationReport.SNRs = zeros(obj.gridSize([2, 1, 3]));
            visibleCellIDs = 1 : size(cellPositions, 2); %% ayar
            numberOfCells = length(visibleCellIDs);
            targetCellIDs = find(visibleZone);
            targetCellIDs = targetCellIDs(visibleCellIDs);
            obj.coverageSimulationReport.targetCellID = targetCellIDs;
            targetPositions = zeros(3, numberOfCells);
            targetPositions(dimensions, :) = cellPositions(dimensions, visibleCellIDs);
            for targetID = 1 : numberOfCells
                targetCellID = targetCellIDs(targetID);
                targetNeighbourCellIDs = obj.neighbours(targetCellID);
                obj.interfaces.settargets(target( ...
                    "position", targetPositions(:, targetID), ...
                    "meanRCS_dbms", options.meanRCS_dbms));
                obj.coverageSimulationReport.thresholds(targetCellID) = obj.configuration.threshold_dB;
                switch obj.network.networkMode
                    case 'multiStatic'
                        obj.coverageSimulationReport.SNRs(targetCellID) = obj.outputSNR_dB;
                    case 'monoStatic'
                        obj.coverageSimulationReport.SNRs(targetCellID) = obj.outputSNR_dB{1};
                end
                for mcID = 1 : mc.numberOfTrials
                    if ~options.onCellCenters
                        obj.interfaces.settargetpositions("width", width);
                    end
                    switch obj.detectionAlgorithm
                        case {"thresholding", "peak"}
                            obj.setmatchfilteredsignals;
                        case {"NCOMP", "OMP"}
                            obj.applycompression;
                        case "CoSaMP"
                    end
                    idx = obj.hypothesisTestingResults;
                    switch obj.network.networkMode
                        case 'multiStatic'
                            for mcpID = 1 : mc.numberOfTrialsParallel
                                if any(ismember(idx{mcpID}, targetNeighbourCellIDs))
                                    obj.coverageSimulationReport.PD(targetCellID) = obj.coverageSimulationReport.PD(targetCellID) + 1/numberOfTotalTrials;
                                end
                            end
                        case 'monoStatic'
                            for mcpID = 1 : mc.numberOfTrialsParallel
                                if any(ismember(idx{1}{mcpID}, targetNeighbourCellIDs))
                                    obj.coverageSimulationReport.PD(targetCellID) = obj.coverageSimulationReport.PD(targetCellID) + 1/numberOfTotalTrials;
                                end
                            end
                    end
                    if ~mod(mcID, 10)
                        fprintf('trial = %d/%d\n', mcID, mc.numberOfTrials);
                    end
                end
                if ~mod(targetID, 10)
                    fprintf('target = %d/%d\n', targetID, numberOfCells);
                end
            end
            function cleanupFunction(obj, originalTargets, originalNumberOfTargets)
                obj.interfaces.targets = originalTargets;
                obj.interfaces.settargetpositions("width", 0);
                obj.configure("numberOfTargets", originalNumberOfTargets);
            end
        end

        function visualizecoveragesimulation(obj)
            gridScan = obj.gridPointsMesh;
            dimensions = {"x", "y", "z"};
            dimensions = dimensions(size(gridScan.x) ~= 1);
            xLabel = dimensions{1} + " (km)";
            yLabel = dimensions{2} + " (km)";
            x1 = obj.gridPoints{1}/1e3;
            x2 = obj.gridPoints{2}/1e3;
            visibleZone = reshape(~obj.blindZone, obj.gridSize([2, 1, 3]));
            figure; img = imagesc(x1, x2, obj.coverageSimulationReport.PD);
            colorbar; colormap('default'); clim([0 1]);
            ax = gca; set(ax, 'Ydir', 'Normal');
            set(img, 'AlphaData', visibleZone);
            delete(datatip(img, 2, 2));
            grid on; grid minor;
            xlabel(xLabel); ylabel(yLabel); zlabel('p_D');
            title('Probability of detection'); hold off;
            img.DataTipTemplate.DataTipRows(1).Label = "x";
            img.DataTipTemplate.DataTipRows(1).Value = gridScan.x;
            img.DataTipTemplate.DataTipRows(2).Label = "y";
            img.DataTipTemplate.DataTipRows(2).Value = gridScan.y;
            img.DataTipTemplate.DataTipRows(3).Label = "PD";
            img.DataTipTemplate.DataTipRows(3).Value = obj.coverageSimulationReport.PD;

            figure; img = imagesc(x1, x2, obj.coverageSimulationReport.SNRs);
            colorbar; colormap('default');
            ax = gca; set(ax, 'Ydir', 'Normal');
            set(img, 'AlphaData', visibleZone);
            delete(datatip(img, 2, 2));
            grid off; grid on; grid minor;
            xlabel(xLabel); ylabel(yLabel); zlabel('SNR (dB)');
            title('SNR'); hold off;
            img.DataTipTemplate.DataTipRows(1).Label = "x";
            img.DataTipTemplate.DataTipRows(1).Value = gridScan.x;
            img.DataTipTemplate.DataTipRows(2).Label = "y";
            img.DataTipTemplate.DataTipRows(2).Value = gridScan.y;
            img.DataTipTemplate.DataTipRows(3).Label = "SNR";
            img.DataTipTemplate.DataTipRows(3).Value = obj.coverageSimulationReport.SNRs;

            figure;
            plot(obj.coverageSimulationReport.SNRs(:), obj.coverageSimulationReport.PD(:), '.');
            grid off; grid on; grid minor;
            xlabel('SNR'); ylabel('p_D');
            title('p_D vs SNR');
        end

        %%%

        function pb_call(obj, varargin)
            S = varargin{3};
            obj.isStep = 0;
            if obj.isWait
                uiresume(S.f);
                obj.isWait = 0;
                set(S.pb, 'String', char(9208));
            else
                obj.isWait = 1;
                set(S.pb, 'String', char(9205));
            end
        end

        function sb_call(obj, varargin)
            S = varargin{3};
            uiresume(S.f);
            obj.isStep = 1;
            obj.isWait = 1;
            set(S.pb, 'String', char(9205));
        end

        function rb_call(obj, varargin)
            S = varargin{3};
            uiresume(S.f);
            obj.isRestart = 1;
        end

        function eb_call(obj, varargin)
            S = varargin{3};
            uiresume(S.f);
            if obj.isSimulation
                obj.isSimulation = 0;
            else
                close(S.f);
            end
        end

        function detection = startsimulation(obj, options)
            arguments
                obj
                options.duration (1, 1) double {mustBePositive} = 1
            end
            stepTime = obj.network.beamTime;
            numberOfSteps = ceil(options.duration ./stepTime);
            currentNumberOfTrials = obj.monteCarlo.numberOfTrials;
            currentNumberOfTrialsParallel = obj.monteCarlo.numberOfTrialsParallel;
            obj.configure("numberOfTrials", 1, "numberOfTrialsParallel", 1);
            obj.resetsimulation;
            figID = 19234;
            h = figure(figID);
            if exist('h', 'var') && ishandle(h)
                close(19234);
            end
            obj.isWait = 0;
            obj.isStep = 0;
            obj.isRestart = 0;
            obj.isSimulation = 1;
            S.f = figure(figID); hold on; view(0, 90);
            set(S.f, 'WindowState', 'maximized');
            S.pb = uicontrol( ...
                'Style', 'pushbutton', ...
                'Units', 'normalized', ...
                'Position', [0.92 0.85 0.05 0.05], ...
                'FontSize', 14);
            S.sb = uicontrol( ...
                'Style', 'pushbutton', ...
                'Units', 'normalized', ...
                'Position', [0.92 0.78 0.05 0.05], ...
                'FontSize', 14);
            S.rb = uicontrol( ...
                'Style', 'pushbutton', ...
                'Units', 'normalized', ...
                'Position', [0.92 0.71 0.05 0.05], ...
                'FontSize', 14);
            S.eb = uicontrol( ...
                'Style', 'pushbutton', ...
                'Units', 'normalized', ...
                'Position', [0.92 0.64 0.05 0.05], ...
                'FontSize', 14);
            set(S.pb, 'Callback', {@obj.pb_call, S});
            set(S.sb, 'Callback', {@obj.sb_call, S});
            set(S.rb, 'Callback', {@obj.rb_call, S});
            set(S.eb, 'Callback', {@obj.eb_call, S});
            set(S.pb, 'String', char(9208));
            set(S.sb, 'String', char(8618));
            set(S.rb, 'String', char(9198));
            set(S.eb, 'String', char(10006));
            posRX = [obj.network.activeReceivingNodes.position]/1e3;
            posTX = [obj.network.activeTransmittingNodes.position]/1e3;
            posRXrep = repelem(posRX, 1, obj.interfaces.numberOfTransmittingNodes);
            posTXrep = repmat(posTX, 1, obj.interfaces.numberOfReceivingNodes);

            numberOfActiveTransmittingNodes = obj.network.numberOfActiveTransmittingNodes;
            numberOfActiveReceivingNodes = obj.network.numberOfActiveReceivingNodes;

            colors = {'m', 'c'};
            phi = linspace(-180, 180, 1000);
            d = max(sqrt(sum(obj.network.boundary.^2)));
            interfacesTX = interface.empty;
            interfacesRX = interface.empty;
            visualizationNormalizationTX = zeros(1, numberOfActiveTransmittingNodes);
            for txID = 1 : numberOfActiveTransmittingNodes
                targetsTemp = target( ...
                    'position', [ ...
                    d*cosd(phi) + posTX(1, txID)*1e3; ...
                    d*sind(phi) + posTX(2, txID)*1e3; ...
                    zeros(1, 1000) + posTX(3, txID)*1e3]);
                interfacesTX(txID) = obj.network.getinterface(targetsTemp);
                visualizationNormalizationTX(txID) = max(interfacesTX(txID).distanceTX(txID, :, :), [], 'all')/1e3;
            end
            visualizationNormalizationRX = zeros(1, numberOfActiveReceivingNodes);
            for rxID = 1 : numberOfActiveReceivingNodes
                targetsTemp = target( ...
                    'position', [ ...
                    d*cosd(phi) + posRX(1, rxID)*1e3; ...
                    d*sind(phi) + posRX(2, rxID)*1e3; ...
                    zeros(1, 1000) + posRX(3, rxID)*1e3]);
                interfacesRX(rxID) = obj.network.getinterface(targetsTemp);
                visualizationNormalizationRX(rxID) = max(interfacesRX(rxID).distanceRX(:, rxID, :), [], 'all')/1e3;
            end

            plot3(posRX(1, :), posRX(2, :), posRX(3, :), 'vb', 'LineWidth', 2, 'MarkerSize', 5);
            plot3(posTX(1, :), posTX(2, :), posTX(3, :), 'vr', 'LineWidth', 2, 'MarkerSize', 5);
            line([posRXrep(1, :); posTXrep(1, :)], [posRXrep(2, :); posTXrep(2, :)], [posRXrep(3, :); posTXrep(3, :)], 'lineStyle', '--', 'Color', 'k');

            ax = gca;
            lines = ax.Children;

            grid off; grid on; grid minor; title('scenario');
            xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
            xlim(obj.network.boundaryListened(1, :)/1e3);
            ylim(obj.network.boundaryListened(2, :)/1e3);
            zlim(obj.network.boundaryListened(3, :)/1e3);

            switch obj.network.networkMode
                case "multiStatic"
                    detection = struct( ...
                        'power', [], ...
                        'x', [], ...
                        'y', [], ...
                        'z', [], ...
                        'numberOfDetections', [], ...
                        'position', []);
                    detection = repmat(detection, 1, numberOfSteps);
                case "monoStatic"
                    detection = struct( ...
                        'power', [], ...
                        'range', [], ...
                        'elevation', [], ...
                        'azimuth', [], ...
                        'numberOfDetections', [], ...
                        'position', []);
                    detection = repmat({repmat(detection, 1, numberOfSteps)}, obj.network.numberOfActiveReceivingNodes, 1);
            end
            numberOfLinesInitial = length(ax.Children);
            stepID = 0;
            start = true;
            while obj.isSimulation
                if stepID == numberOfSteps
                    set(S.pb, 'String', char(9205));
                    uiwait(S.f);
                end
                if obj.isRestart && stepID
                    obj.isRestart = 0;
                    obj.resetsimulation;
                    stepID = 0;
                    if ~obj.isStep
                        obj.isWait = 0;
                        set(S.pb, 'String', char(9208));
                    end
                    continue
                end
                if stepID == numberOfSteps
                    continue
                end
                figure(figID);
                lineID = numberOfLinesInitial + 1;
                if ~start
                    if ~isempty(obj.interfaces.targets)
                        obj.interfaces.targets.step(stepTime);
                    end
                    for rxID = 1 : numberOfActiveReceivingNodes
                        obj.network.activeReceivingNodes(rxID).array.step(stepTime);
                    end
                    for txID = 1 : numberOfActiveTransmittingNodes
                        obj.network.transmittingNodes(txID).array.step(stepTime);
                    end
                end
                obj.setmatchfilteredsignals;
                switch obj.network.networkMode
                    case "multiStatic"
                        detections = obj.detectionFromIntegration;
                    case "monoStatic"
                        detections = obj.detectionFromMatchFiltration;
                end
                if ~start && sum([detections.numberOfDetections]) == 0
                    continue;
                end
                if (obj.isWait || obj.isStep) && stepID
                    uiwait(S.f);
                    if obj.isRestart
                        obj.isRestart = 0;
                        obj.resetsimulation;
                        stepID = 0;
                        continue
                    end
                end

                if ~isempty(obj.interfaces.targets)
                    x = obj.interfaces.targets.position(1, :)/1e3;
                    y = obj.interfaces.targets.position(2, :)/1e3;
                    z = obj.interfaces.targets.position(3, :)/1e3;
    
                    if start
                        if obj.interfaces.numberOfTargets < 11
                            plot3(x, y, z, '+k', 'LineWidth', 2, 'MarkerSize', 10);
                        elseif obj.interfaces.numberOfTargets < 21
                            plot3(x, y, z, '*k', 'LineWidth', 1, 'MarkerSize', 10);
                        else
                            plot3(x, y, z, '.k');
                        end
                    else
                        lines(lineID).XData = x;
                        lines(lineID).YData = y;
                        lines(lineID).ZData = z;
                        lineID = lineID + 1;
                    end
                end

                switch obj.network.networkMode
                    case "monoStatic"
                        for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                            detection{rxID}(stepID + 1) = detections(rxID);
                            pos = detection{rxID}(stepID + 1).position;
                            if isempty(pos)
                                if start
                                    % dummy to create line
                                    if obj.network.numberOfActiveReceivingNodes == 2
                                        plot3(nan, nan, nan, ['.', colors{rxID}], 'LineWidth', 1, 'MarkerSize', 10);
                                    else
                                        plot3(nan, nan, nan, '.c', 'LineWidth', 1, 'MarkerSize', 10);
                                    end
                                    continue;
                                else
                                    lines(lineID).XData = [];
                                    lines(lineID).YData = [];
                                    lines(lineID).ZData = [];
                                    lineID = lineID + 1;
                                    continue
                                end
                            end
                            xEstimation = pos(1, :)/1e3;
                            yEstimation = pos(2, :)/1e3;
                            zEstimation = pos(3, :)/1e3;
                            if start
                                if obj.network.numberOfActiveReceivingNodes == 2
                                    plot3(xEstimation, yEstimation, zEstimation, ['.', colors{rxID}], 'LineWidth', 1, 'MarkerSize', 10);
                                else
                                    plot3(xEstimation, yEstimation, zEstimation, '.c', 'LineWidth', 1, 'MarkerSize', 10);
                                end
                            else
                                lines(lineID).XData = xEstimation;
                                lines(lineID).YData = yEstimation;
                                lines(lineID).ZData = zEstimation;
                                lineID = lineID + 1;
                            end
                        end
                    case "multiStatic"
                        detection(stepID + 1) = detections;
                        pos = detection(stepID + 1).position;
                        if isempty(pos)
                            if start
                                % dummy to create line
                                plot3(nan, nan, nan, '.r', 'LineWidth', 1, 'MarkerSize', 10);
                                continue;
                            else
                                lines(lineID).XData = [];
                                lines(lineID).YData = [];
                                lines(lineID).ZData = [];
                                continue
                            end
                        end
                        xEstimation = pos(1, :)/1e3;
                        yEstimation = pos(2, :)/1e3;
                        zEstimation = pos(3, :)/1e3;
                        if start
                            plot3(xEstimation, yEstimation, zEstimation, '.r', 'LineWidth', 1, 'MarkerSize', 10);
                        else
                            lines(lineID).XData = xEstimation;
                            lines(lineID).YData = yEstimation;
                            lines(lineID).ZData = zEstimation;
                        end
                end

                if (~strcmpi(obj.network.surveillanceMode, "staticBeam") && ~start) || start
                    for txID = 1 : numberOfActiveTransmittingNodes
                        G = visualizationNormalizationTX(txID).*permute(abs(interfacesTX(txID).transmittedBeam(txID, :, :).^2./interfacesTX(txID).network.activeTransmittingNodes(txID).array.numberOfTotalElements), [1 3 2]);
                        u = permute(interfacesTX(txID).unitDirectionTX, [1 3 2]);
                        if start
                            plot3(G.*u(1, :, txID) + posTX(1, txID), G.*u(2, :, txID) + posTX(2, txID), u(3, :, txID) + posTX(3, txID), 'r');
                        else
                            lines(lineID).XData = G.*u(1, :, txID) + posTX(1, txID);
                            lines(lineID).YData = G.*u(2, :, txID) + posTX(2, txID);
                            lines(lineID).ZData = u(3, :, txID) + posTX(3, txID);
                            lineID = lineID + 1;
                        end
                    end
                    for rxID = 1 : numberOfActiveReceivingNodes
                        switch obj.network.networkMode
                            case "monoStatic"
                                txID = obj.network.monoStaticTransmitterIDs(rxID);
                                if isnan(txID)
                                    break;
                                end
                            case "multiStatic"
                                txID = 1;
                        end
                        G = visualizationNormalizationRX(rxID).*permute(abs(interfacesRX(rxID).receivedBeamSpaceObservations{rxID}(txID, :, :, :).^2./interfacesRX(rxID).network.activeReceivingNodes(rxID).array.numberOfTotalElements), [3 4 1 2]);
                        u = permute(interfacesRX(rxID).unitDirectionRX, [3 1 2]);
                        if start
                            if obj.network.numberOfActiveReceivingNodes == 2
                                plot3(G.*u(:, 1, rxID) + posRX(1, rxID), G.*u(:, 2, rxID) + posRX(2, rxID), u(:, 3, rxID) + posRX(3, rxID), colors{rxID}, 'LineWidth', 2);
                            else
                                plot3(G.*u(:, 1, rxID) + posRX(1, rxID), G.*u(:, 2, rxID) + posRX(2, rxID), u(:, 3, rxID) + posRX(3, rxID), 'b', 'LineWidth', 2);
                            end
                        else
                            switch obj.network.beamformingMode
                                case 'bypass'
                                    numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).array.numberOfTotalElements;
                                otherwise
                                    numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
                            end
                            for chID = 1 : numberOfTotalChannels
                                lines(lineID).XData = G(:, chID).*u(:, 1, rxID) + posRX(1, rxID);
                                lines(lineID).YData = G(:, chID).*u(:, 2, rxID) + posRX(2, rxID);
                                lines(lineID).ZData = u(:, 3, rxID) + posRX(3, rxID);
                                lineID = lineID + 1;
                            end
                        end
                    end
                    for txID = 1 : numberOfActiveTransmittingNodes
                        n = visualizationNormalizationTX(txID).*obj.network.activeTransmittingNodes(txID).array.normalVector/2;
                        if start
                            quiver3(posTX(1, txID), posTX(2, txID), posTX(3, txID), n(1), n(2), n(3), 'Color', 'k');
                        else
                            lines(lineID).XData = posTX(1, txID);
                            lines(lineID).YData = posTX(2, txID);
                            lines(lineID).ZData = posTX(3, txID);
                            lines(lineID).UData = n(1);
                            lines(lineID).VData = n(2);
                            lines(lineID).WData = n(3);
                            lineID = lineID + 1;
                        end
                    end
                    for rxID = 1 : numberOfActiveReceivingNodes
                        n = visualizationNormalizationRX(rxID).*obj.network.activeReceivingNodes(rxID).array.normalVector/2;
                        if start
                            quiver3(posRX(1, rxID), posRX(2, rxID), posRX(3, rxID), n(1), n(2), n(3), 'Color', 'k');
                        else
                            lines(lineID).XData = posRX(1, rxID);
                            lines(lineID).YData = posRX(2, rxID);
                            lines(lineID).ZData = posRX(3, rxID);
                            lines(lineID).UData = n(1);
                            lines(lineID).VData = n(2);
                            lines(lineID).WData = n(3);
                            lineID = lineID + 1;
                        end
                    end
                end
                % for txID = 1 : numberOfActiveTransmittingNodes
                %     s = 1.5*visualizationNormalizationTX(txID).*obj.network.activeTransmittingNodes(txID).array.steeringUnitDirection;
                %     if start
                %         quiver3(posTX(1, txID), posTX(2, txID), posTX(3, txID), s(1), s(2), s(3), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
                %     else
                %         lines(lineID).XData = posTX(1, txID);
                %         lines(lineID).YData = posTX(2, txID);
                %         lines(lineID).ZData = posTX(3, txID);
                %         lines(lineID).UData = s(1);
                %         lines(lineID).VData = s(2);
                %         lines(lineID).WData = s(3);
                %         lineID = lineID + 1;
                %     end
                % end
                % for rxID = 1 : numberOfActiveReceivingNodes
                %     s = 1.5*visualizationNormalizationRX(rxID).*obj.network.activeReceivingNodes(rxID).array.steeringUnitDirection;
                %     if start
                %         if obj.network.numberOfActiveReceivingNodes == 2
                %             quiver3(posRX(1, rxID), posRX(2, rxID), posRX(3, rxID), s(1), s(2), s(3), 'Color', colors{rxID}, 'LineStyle', '--', 'LineWidth', 2);
                %         else
                %             quiver3(posRX(1, rxID), posRX(2, rxID), posRX(3, rxID), s(1), s(2), s(3), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
                %         end
                %     else
                %         lines(lineID).XData = posRX(1, rxID);
                %         lines(lineID).YData = posRX(2, rxID);
                %         lines(lineID).ZData = posRX(3, rxID);
                %         lines(lineID).UData = s(1);
                %         lines(lineID).VData = s(2);
                %         lines(lineID).WData = s(3);
                %         lineID = lineID + 1;
                %     end
                % end
                
                if start
                    lines = flipud(ax.Children);
                    start = false;
                end
                if stepID < numberOfSteps
                    stepID = stepID + 1;
                end
                drawnow;
            end
            obj.configure( ...
                "numberOfTrials", currentNumberOfTrials, ...
                "numberOfTrialsParallel", currentNumberOfTrialsParallel);
        end

        function resetsimulation(obj)
            if ~isempty(obj.interfaces.targets)
                obj.interfaces.targets.reset;
            end
            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                obj.network.activeReceivingNodes(rxID).array.reset;
            end
            for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                obj.network.transmittingNodes(txID).array.reset;
            end
        end
        
        %%% visualization methods

        function visualizedictionary(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.cellInterested double {mustBeInteger, mustBePositive} = []
                options.plot (1, 1) string {mustBeMember(options.plot, ["magnitude", "phase"])} = "magnitude"
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            if isempty(options.cellInterested)
                dictionary = obj.dictionaryCompression(:, :, 1);
                switch options.plot
                    case "magnitude"
                        imagesc(abs(dictionary)); colorbar;
                        xlabel('cells'); ylabel('sample');
                        title('magnitude of dictionary');
                    case "phase"
                        imagesc(180*angle(dictionary)/pi); colorbar;
                        xlabel('cells'); ylabel('sample');
                        title('phase of dictionary');
                end
            elseif obj.network.numberOfActiveReceivingNodes > 1
                dictionary = squeeze(obj.dictionaryCompression(:, options.cellInterested, :));
                switch options.plot
                    case "magnitude"
                        imagesc(abs(dictionary)); colorbar;
                        xlabel('rxID'); ylabel('sample');
                        title('magnitude of dictionary');
                    case "phase"
                        imagesc(180*angle(dictionary)/pi); colorbar;
                        xlabel('rxID'); ylabel('sample');
                        title('phase of dictionary');
                end
            else
                dictionary = squeeze(obj.dictionaryCompression(:, options.cellInterested, :));
                switch options.plot
                    case "magnitude"
                        plot(abs(dictionary), 'b');
                        xlabel('sample'); ylabel('magnitude');
                        title('magnitude of dictionary');
                    case "phase"
                        plot(180*angle(dictionary)/pi, '*r');
                        xlabel('sample'); ylabel('phase (degree)');
                        title('phase of dictionary');
                end
            end
        end

        function visualizefilteredsignals(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.receivingNodeIDs (1, :) double {mustBeInteger, mustBeNonnegative} = 1 : obj.network.numberOfActiveReceivingNodes
                options.trialID (1, 1) double {mustBePositive, mustBeInteger} = 1
            end
            config = obj.configuration;
            mustBeInRange(options.receivingNodeIDs, 1, obj.network.numberOfReceivingNodes);
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            if isscalar(obj.signalsMatchFiltered) || isempty(obj.signalsMatchFiltered)
                obj.setmatchfilteredsignals;
            end
            for rxID = options.receivingNodeIDs
                s = obj.signalsMatchFiltered{rxID}(:, :, options.trialID);
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
                            L = obj.network.pulseWidthSample(rxID, txID) - 2;
                            N = obj.network.activeReceivingNodes(rxID).numberOfSamplesPerCPI;
                    end
                    t = (-L : N)*Ts;
                    if ~isscalar(options.receivingNodeIDs)
                        plot(t, 20*log10(abs(s(:, ceil(end/2)))));
                    else
                        plot(t, 20*log10(abs(s)));
                    end
                    hold on;
                    if strcmpi(obj.network.networkMode, "monoStatic") && obj.processingAlgorithm == 2 && config.CFAR
                        T = obj.thresholdCFAR{rxID}(:, :, options.trialID);
                        if ~isscalar(options.receivingNodeIDs)
                            plot(t, 10*log10(T(:, ceil(end/2))), 'LineStyle', '--', 'LineWidth', 2);
                        else
                            plot(t, 10*log10(T), 'LineStyle', '--', 'LineWidth', 2);
                        end
                    end
                end
            end
            if strcmpi(obj.network.networkMode, "monoStatic")
                if obj.processingAlgorithm ~= 2 || ~config.CFAR
                    yline(config.threshold_dB, 'LineStyle', '--', 'LineWidth', 2);
                end
            end
            grid off; grid on; grid minor;
            xlabel('time (s)'); ylabel('power (dB)');
            title('filtered signal');
            if ~isscalar(options.receivingNodeIDs) && ~(strcmpi(obj.network.networkMode, "monoStatic") && obj.processingAlgorithm == 2 && config.CFAR)
                leg = legend(num2str(options.receivingNodeIDs.'), 'Location', 'best');
                title(leg, 'RX ID');
            end
            drawnow; hold off;
        end

        function visualizehypothesizedtimedelays(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.transmittingNodeID (1, 1) double {mustBeInteger, mustBePositive} = 1
                options.receivingNodeID (1, 1) double {mustBeInteger, mustBePositive} = 1
                options.dimension (1, 1) string {mustBeMember(options.dimension, ["x-y", "y-z", "z-x", "x-y-z"])} = "x-y"
            end
            mustBeInRange(options.transmittingNodeID, 1, obj.network.numberOfActiveTransmittingNodes);
            mustBeInRange(options.receivingNodeID, 1, obj.network.numberOfActiveReceivingNodes);
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            gridScan = obj.gridPointsMesh;
            tau = obj.hypothesizedTimeDelays(options.transmittingNodeID, options.receivingNodeID, :)*1e6;
            switch options.dimension
                case "x-y"
                    tau = reshape(squeeze(tau), [obj.gridSize]);
                    x = gridScan.x(:, :, 1); x = x(:)/1e3;
                    y = gridScan.y(:, :, 1); y = y(:)/1e3;
                    t = tau(:, :, 1);
                    plot3(x, y, t(:), '.'); hold on;
                    if size(tau, 3) ~= 1
                        for zID = unique(ceil(linspace(2, size(tau, 3), 7)))
                            t = tau(:, :, zID);
                            plot3(x, y, t(:), '.');
                        end
                    end
                    grid off; grid on; grid minor;
                    xlabel('x (km)'); ylabel('y (km');
                    zlabel('time delay (\mus)');
                    title('hypothesized time delays');
                case "y-z"
                    tau = reshape(squeeze(tau), [obj.gridSize]).';
                    y = gridScan.y(:, 1, :); y = y(:)/1e3;
                    z = gridScan.z(:, 1, :); z = z(:)/1e3;
                    t = tau(:, 1, :);
                    plot3(y, z, t(:), '.'); hold on;
                    if size(tau, 2) ~= 1
                        for xID = unique(ceil(linspace(2, size(tau, 2), 7)))
                            t = tau(:, xID, :);
                            plot3(y, z, t(:), '.');
                        end
                    end
                    grid off; grid on; grid minor;
                    xlabel('y (km)'); ylabel('z (km)');
                    zlabel('time delay (\mus)');
                    title('hypothesized time delays');
                case "z-x"
                    tau = reshape(squeeze(tau), [obj.gridSize]).';
                    z = gridScan.z(1, :, :); z = z(:)/1e3;
                    x = gridScan.x(1, :, :); x = x(:)/1e3;
                    t = tau(1, :, :);
                    plot3(z, x, t(:), '.'); hold on;
                    if size(tau, 1) ~= 1
                        for yID = unique(ceil(linspace(2, size(tau, 1), 7)))
                            t = tau(yID, :, :);
                            plot3(z, x, t(:), '.');
                        end
                    end
                    grid off; grid on; grid minor;
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
                    grid off; grid on; grid minor;
                    xlabel('x (meter)'); ylabel('y (meter)'); zlabel('z (meter)');
                    title('hypothesized locations');
            end
        end

        function visualizeintegrationindices(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.receivingNodeIDs (1, :) double {mustBeInteger, mustBeNonnegative} = 1 : obj.network.numberOfActiveReceivingNodes
            end
            mustBeInRange(options.receivingNodeIDs, 1, obj.network.numberOfReceivingNodes);
            if isempty(obj.integrationIndices)
                fprintf('integration indices are not set');
                return;
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            for rxID = options.receivingNodeIDs
                switch obj.network.networkMode
                    case "multiStatic"
                        plot(squeeze(obj.integrationIndices(:, rxID, :)).', '.');
                    case "monoStatic"
                        plot(squeeze(obj.integrationIndices(obj.network.monoStaticTransmitterIDs(rxID), rxID, :)).', '.');
                end
                hold on;
            end
            xlabel('grid point'); ylabel('index');
            xlim tight; ylim tight;
            grid off; grid on; grid minor;
            drawnow; hold off;
        end

        function visualizeintegratedsignals(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.scale (1, 1) string {mustBeMember(options.scale, ["linear", "power", "dB", "SNR"])} = "dB"
                options.plotMode (1, 1) string {mustBeMember(options.plotMode, ["image", "mesh"])} = "image"
                options.trialID (1, 1) double {mustBePositive, mustBeInteger} = 1
                options.monoStaticNetworkRXID (1, 1) {mustBePositive, mustBeInteger} = 1
                options.monoStaticNetworkCHID (1, 1) {mustBePositive, mustBeInteger} = 1
                options.iterationID (1, :) {mustBeNonnegative, mustBeInteger} = 0
            end
            if isnan(obj.integrationWeights)
                return;
            end
            config = obj.configuration;
            options.iterationID = unique(min(options.iterationID, max([obj.compressionReport.numberOfIterations])));
            switch obj.network.networkMode
                case "multiStatic"
                    if options.iterationID 
                        pos = obj.estimatedPositions{options.trialID};
                        if ~isempty(pos)
                            pos = pos(:, unique(min(size(pos, 2), options.iterationID)));
                        end
                    else
                        pos = obj.estimatedPositions{options.trialID};
                    end
                case "monoStatic"
                    pos = obj.estimatedPositions{options.monoStaticNetworkRXID}{options.monoStaticNetworkCHID, options.trialID};
            end
            if any(options.iterationID > 0) && ~isempty(obj.compressionReport(options.trialID).integratedSignals)
                z = nan(prod(obj.gridSize), 1);
                if isscalar(options.iterationID)
                    z(~obj.blindZone) = obj.compressionReport(options.trialID).integratedSignals(:, options.iterationID);
                else
                    for plotID = 1 : numel(options.iterationID)
                        obj.visualizeintegratedsignals( ...
                            "figureID", options.figureID + plotID - 1, ...
                            "scale", options.scale, ...
                            "plotMode", options.plotMode, ...
                            "iterationID", options.iterationID(plotID), ...
                            "monoStaticNetworkCHID", options.monoStaticNetworkCHID, ...
                            "monoStaticNetworkRXID", options.monoStaticNetworkRXID, ...
                            "trialID", options.trialID);
                    end
                    return;
                end
            else
                z = obj.signalsIntegrated;
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            gridScan = obj.gridPointsMesh;
            dims = find(size(gridScan.x) ~= 1);
            dimensions = {"x", "y", "z"};
            dimensions = dimensions(dims);
            switch length(dims)
                case 1
                    xLabel = dimensions{1} + " (km)";
                    x1 = gridScan.(dimensions{1})/1e3;
                    switch obj.network.networkMode
                        case "multiStatic"
                            switch obj.processingAlgorithm
                                case {3, 5, 6}
                                    plot(x1, 10*log10(z(:, :, options.trialID)));
                                otherwise
                                    plot(x1, 20*log10(z(:, :, options.trialID)));
                            end
                        case "monoStatic"
                            switch obj.processingAlgorithm
                                case {3, 5, 6}
                                    plot(x1, 10*log10(z{options.monoStaticNetworkRXID}(:, :, options.trialID)));
                                otherwise
                                    plot(x1, 20*log10(z{options.monoStaticNetworkRXID}(:, :, options.trialID)));
                            end
                    end
                    grid off; grid on; grid minor;
                    xlabel(xLabel); ylabel('power (dB)');
                    title('integrated signal');
                case 2
                    zLabel = 'z_{processed}';
                    switch obj.network.networkMode
                        case "multiStatic"
                            Y = reshape(z(:, :, options.trialID), obj.gridSize([2, 1, 3]));
                        case "monoStatic"
                            if iscell(z)
                                Y = reshape(z{options.monoStaticNetworkRXID}(:, options.monoStaticNetworkCHID, options.trialID), obj.gridSize([2, 1, 3]));
                            else
                                Y = reshape(z(:, :, options.trialID), obj.gridSize([2, 1, 3]));
                            end
                    end
                    switch options.scale
                        case "linear"
                            zLabel = [zLabel, ', linear'];
                            switch obj.processingAlgorithm
                                case {3, 5, 6}
                                    Y = sqrt(Y);
                            end
                        case "power"
                            zLabel = [zLabel, ', power'];
                            switch obj.processingAlgorithm
                                case {3, 5, 6}
                                otherwise
                                    Y = Y.^2;
                            end
                        case "dB"
                            zLabel = [zLabel, ', power (dB)'];
                            switch obj.processingAlgorithm
                                case {3, 5, 6}
                                    Y = 10*log10(abs(Y));
                                otherwise
                                    Y = 20*log10(abs(Y));
                            end
                        case "SNR"
                            zLabel = [zLabel, ', SNR (dB)'];
                            switch obj.processingAlgorithm
                                case {3, 5, 6}
                                    Y = 10*log10(abs(Y)) - unique([obj.interfaces.network.activeReceivingNodes.noisePower]);
                                otherwise
                                    Y = 20*log10(abs(Y)) - unique([obj.interfaces.network.activeReceivingNodes.noisePower]);
                            end
                    end
                    maxY = max(Y, [], 'all');
                    xLabel = dimensions{1} + " (km)";
                    yLabel = dimensions{2} + " (km)";
                    x1 = gridScan.(dimensions{1})/1e3;
                    x2 = gridScan.(dimensions{2})/1e3;
                    posRX = [obj.network.activeReceivingNodes.position]/1e3;
                    posTX = [obj.network.activeTransmittingNodes.position]/1e3;
                    posRX = posRX(dims, :); posTX = posTX(dims, :);
                    switch options.plotMode
                        case "mesh"
                            hold on;
                            plot3(posRX(1, :), posRX(2, :), repmat(maxY, [1, size(posRX, 2)]), 'vb', 'LineWidth', 2);
                            plot3(posTX(1, :), posTX(2, :), repmat(maxY, [1, size(posTX, 2)]), 'vr', 'LineWidth', 2);
                            estimations = pos(dims, :)/1e3;
                            plot3(estimations(1, :), estimations(2, :), repmat(maxY, [1, size(estimations, 2)]), 'om', 'LineWidth', 2);
                            if ~isempty(obj.interfaces.targets)
                                x = obj.interfaces.targets.position(dims, :)/1e3;
                                plot3(x(1, :), x(2, :), repmat(maxY, [1, size(x, 2)]), '+k', 'LineWidth', 2);
                            end
                            plot3(x2(:), x1(:), config.threshold_dB(options.monoStaticNetworkRXID)*ones(1, numel(x2)), 'k');
                            m = mesh(x2, x1, Y);
                            set(m, 'AlphaData', ~isinf(Y) & ~isnan(Y));
                            m.FaceColor = 'flat';
                        case "image"
                            hold on;
                            img = imagesc(obj.gridPoints{1}/1e3, obj.gridPoints{2}/1e3, Y);
                            set(img, 'AlphaData', ~isinf(Y) & ~isnan(Y));
                            delete(datatip(img, 2, 2));
                            img.DataTipTemplate.DataTipRows(1).Label = "x";
                            img.DataTipTemplate.DataTipRows(1).Value = x1;
                            img.DataTipTemplate.DataTipRows(2).Label = "y";
                            img.DataTipTemplate.DataTipRows(2).Value = x2;
                            img.DataTipTemplate.DataTipRows(3).Label = "power";
                            img.DataTipTemplate.DataTipRows(3).Value = Y;
                            plot(posRX(1, :), posRX(2, :), 'vb', 'LineWidth', 2);
                            plot(posTX(1, :), posTX(2, :), 'vr', 'LineWidth', 2);
                            estimations = pos(dims, :)/1e3;
                            plot(estimations(1, :), estimations(2, :), 'om', 'LineWidth', 2);
                            if ~isempty(obj.interfaces.targets)
                                x = obj.interfaces.targets.position(dims, :)/1e3;
                                plot(x(1, :), x(2, :), '+k', 'LineWidth', 2);
                            end
                    end
                    colorbar; clim([0 inf]);
                    xlim tight; ylim tight;
                    grid on; grid minor; view(0, 90);
                    xlabel(xLabel); ylabel(yLabel); zlabel(zLabel);
                    switch obj.network.networkMode
                        case "multiStatic"
                            title(sprintf('integrated signal SNR = %.3g dB, PFA = %.3e, PD = %.3g', max(obj.outputSNR_dB), config.PFA, max(config.PD))); hold off;
                        case "monoStatic"
                            title(sprintf('integrated signal SNR = %.3g dB, PFA = %.3e, PD = %.3g', obj.outputSNR_dB{options.monoStaticNetworkRXID}, config.PFA, config.PD(options.monoStaticNetworkRXID))); hold off;
                    end
                case 3
            end
            if isempty(pos)
                if ~isempty(obj.interfaces.targets)
                    legend('RX', 'TX', 'targets', 'Location', 'best');
                else
                    legend('RX', 'TX', 'Location', 'best');
                end
            else
                if ~isempty(obj.interfaces.targets)
                    legend('RX', 'TX', 'estimations', 'targets', 'Location', 'best');
                else
                    legend('RX', 'TX', 'estimations', 'Location', 'best');
                end
            end
            hold off; drawnow;
        end

        function visualizecompression(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.trialID (1, 1) double {mustBePositive, mustBeInteger} = 1
                options.monoStaticNetworkRXID (1, 1) {mustBePositive, mustBeInteger} = 1
                options.monoStaticNetworkCHID (1, 1) {mustBePositive, mustBeInteger} = 1
            end
            config = obj.configuration;
            switch obj.network.networkMode
                case "multiStatic"
                    pos = obj.estimatedPositions{options.trialID};
                case "monoStatic"
                    pos = obj.estimatedPositions{options.monoStaticNetworkRXID}{options.monoStaticNetworkCHID, options.trialID};
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
            z = abs(obj.compressionReport(options.trialID).complexAmplitudes);
            Y = sum(z.^2, 2);
            switch sum(dims)
                case 1
                    xLabel = dimensions{1} + " (km)";
                    x1 = gridScan.(dimensions{1})/1e3;
                    plot(x1, 20*log10(z));
                    xlim tight; ylim tight; grid off; grid on; grid minor;
                    xlabel(xLabel); ylabel('power (dB)');
                    title('compressed signal');
                case 2
                    Y = 10*log10(Y);
                    xLabel = dimensions{1} + " (km)";
                    yLabel = dimensions{2} + " (km)";
                    posRX = [obj.network.activeReceivingNodes.position]/1e3;
                    posTX = [obj.network.activeTransmittingNodes.position]/1e3;
                    posRX = posRX(dims, :); posTX = posTX(dims, :);
                    estimations = pos(dims, :)/1e3;
                    plot(estimations(1, :), estimations(2, :), 'om', 'LineWidth', 2); hold on;
                    plot(posRX(1, :), posRX(2, :), 'vb', 'LineWidth', 2);
                    plot(posTX(1, :), posTX(2, :), 'vr', 'LineWidth', 2);
                    if ~isempty(obj.interfaces.targets)
                        posTargets = obj.interfaces.targets.position(dims, :)/1e3;
                        plot(posTargets(1, :), posTargets(2, :), '+k', 'LineWidth', 2);
                    end
                    text(estimations(1, :), estimations(2, :), num2str(Y), ...
                        'FontSize', 12, ...
                        'Color', 'blue', ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'bottom');
                    xlim tight; ylim tight; grid off; grid on; grid minor;
                    xlabel(xLabel); ylabel(yLabel);
                    switch obj.network.networkMode
                        case "multiStatic"
                            title(sprintf('compressed signal SNR = %.3g dB, PFA = %.3e, PD = %.3g', max(obj.outputSNR_dB(:, :, 1)), config.PFA, max(config.PD))); hold off;
                        case "monoStatic"
                            title(sprintf('compressed signal SNR = %.3g dB, PFA = %.3e, PD = %.3g', obj.outputSNR_dB{options.monoStaticNetworkRXID}, config.PFA, config.PD(options.monoStaticNetworkRXID))); hold off;
                    end
                case 3
            end
            if isempty(pos)
                legend('RX', 'TX', 'targets', 'Location', 'best');
            else
                legend('RX', 'TX', 'estimations', 'targets', 'Location', 'best');
            end
            hold off; drawnow;
            figure; plot(10*log10(obj.compressionReport(options.trialID).residualHistory), 'b');
            hold on; yline(obj.configuration.threshold_dB, 'r--');
            grid off; grid on; grid minor;
            xlim tight; ylim tight;
            xlabel('iteration ID'); ylabel('residual power (dB)');
            title('residual power vs iteration ID of the compressing sensing algorithm');
            legend('residual power', 'threshold', 'Location', 'best');
        end

        function visualizedirectlyintegratedsignals(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.cuttingIndices double = [nan, nan, ceil(([obj.network.activeReceivingNodes(3 : end).numberOfSamplesPerCPI] + obj.network.pulseWidthSample(3 : end, 1).' - 1)/2)]
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            cuttingRxIDs = find(~isnan(options.cuttingIndices));
            visualizedRxIDs = setdiff(1 : obj.network.numberOfActiveReceivingNodes, cuttingRxIDs);
            Y = obj.signalsIntegratedDirectly;
            for rxID = cuttingRxIDs
                Y = indexer(Y, rxID, options.cuttingIndices(rxID));
            end
            txID = 1;
            t = cell(1, obj.network.numberOfActiveReceivingNodes);
            for rxID = visualizedRxIDs
                switch obj.network.activeTransmittingNodes(txID).transmissionType
                    case "continuous"
                        error('not implemented');
                    case "pulsed"
                        L = obj.network.pulseWidthSample(rxID, txID) - 2;
                        N = obj.network.activeReceivingNodes(rxID).numberOfSamplesPerCPI;
                end
                t{rxID} = (-L : N)*obj.network.activeReceivingNodes(rxID).samplingPeriod*1e6;
            end
            mesh(t{visualizedRxIDs(1)}, t{visualizedRxIDs(2)}, 10*log10(squeeze(Y)));
            xlabel(['\tau_{', num2str(visualizedRxIDs(1)), '} (\mus)']); ylabel(['\tau_{', num2str(visualizedRxIDs(2)), '} (\mus)']); zlabel('power (dB)');
            xlim tight; ylim tight; zlim tight;
            grid on; grid minor; colorbar;
        end

        function visualizeestimation(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.trialID (1, 1) double {mustBePositive, mustBeInteger} = 1
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            posRX = [obj.network.activeReceivingNodes.position]/1e3;
            posTX = [obj.network.activeTransmittingNodes.position]/1e3;
            plot3(posRX(1, :), posRX(2, :), posRX(3, :), 'vb', 'LineWidth', 2, 'MarkerSize', 10);
            hold on; plot3(posTX(1, :), posTX(2, :), posTX(3, :), 'vr', 'LineWidth', 2, 'MarkerSize', 10);
            if ~isempty(obj.interfaces.targets)
                x = obj.interfaces.targets.position(1, :)/1e3;
                y = obj.interfaces.targets.position(2, :)/1e3;
                z = obj.interfaces.targets.position(3, :)/1e3;
                plot3(x, y, z, '*k', 'LineWidth', 1, 'MarkerSize', 10);
            end
            switch obj.network.networkMode
                case "multiStatic"
                    est = obj.estimatedPositions;
                    pos = cell(1, obj.network.numberOfActiveReceivingNodes);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        pos{rxID} = [est{rxID}{:, options.trialID}];
                    end
                    pos = [pos{:}];
                    xEstimation = pos(1, :)/1e3;
                    yEstimation = pos(2, :)/1e3;
                    zEstimation = pos(3, :)/1e3;
                    if size(pos, 2) < 9
                        plot3(xEstimation, yEstimation, zEstimation, '+c', 'LineWidth', 2, 'MarkerSize', 10);
                    else
                        plot3(xEstimation, yEstimation, zEstimation, '.c', 'LineWidth', 1, 'MarkerSize', 10);
                    end
                case "monoStatic"
                    colors = {'m', 'c'};
                    est = obj.detectionFromMatchFiltration;
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        pos = est(rxID, options.trialID).position;
                        if isempty(pos)
                            continue
                        end
                        xEstimation = pos(1, :)/1e3;
                        yEstimation = pos(2, :)/1e3;
                        zEstimation = pos(3, :)/1e3;
                        if obj.network.numberOfActiveReceivingNodes == 2
                            if size(pos, 2) < 9
                                plot3(xEstimation, yEstimation, zEstimation, ['+', colors{rxID}], 'LineWidth', 2, 'MarkerSize', 10);
                            else
                                plot3(xEstimation, yEstimation, zEstimation, ['.', colors{rxID}], 'LineWidth', 1, 'MarkerSize', 10);
                            end
                        else
                            if size(pos, 2) < 9
                                plot3(xEstimation, yEstimation, zEstimation, '+c', 'LineWidth', 2, 'MarkerSize', 10);
                            else
                                plot3(xEstimation, yEstimation, zEstimation, '.c', 'LineWidth', 1, 'MarkerSize', 10);
                            end
                        end
                    end
            end
            text(posRX(1, :), posRX(2, :), posRX(3, :), num2str((1 : obj.network.numberOfActiveReceivingNodes).'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            text(posTX(1, :), posTX(2, :), posTX(3, :), num2str((1 : obj.network.numberOfActiveTransmittingNodes).'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            if size(posRX, 2) < 11
                text(x, y, z, num2str((1 : obj.interfaces.numberOfTargets).'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            end
            posRX = repelem(posRX, 1, obj.interfaces.network.numberOfActiveTransmittingNodes);
            posTX = repmat(posTX, 1, obj.interfaces.network.numberOfActiveReceivingNodes);
            line([posRX(1, :); posTX(1, :)], [posRX(2, :); posTX(2, :)], [posRX(3, :); posTX(3, :)], 'lineStyle', '--', 'Color', 'k');
            grid off; grid on; grid minor; title('scenario');
            xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
            legend('RX', 'TX', 'targets', 'estimations', 'Location', 'best');
            view(0, 90); hold off; drawnow;
        end

        function visualizereceiveroperatingcharacteristics(obj, options)
            arguments
                obj
                options.snr_dB (1, :) double = [0 3 10 13] % dB
                options.numberOfReceivers (1, :) double = [1 : 5 10]
                options.scanAlgorithms (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.scanAlgorithms, [0, 1])} = 0
            end
            nop = 1000;
            if options.scanAlgorithms
                PFA = logspace(-10, 0, nop);
                PD = zeros(nop, 6);
                currentAlgorithm = obj.processingAlgorithm;
                for k = 1 : 6
                    obj.configure("processingAlgorithm", k);
                    PD(:, k) = obj.ROC(PFA, 10.^(.1*options.snr_dB(1)), obj.network.numberOfActiveReceivingNodes);
                end
                obj.configure("processingAlgorithm", currentAlgorithm);
                figure; semilogx(PFA, PD);
                grid on; grid minor;
                xlabel('PFA'); ylabel('PD');
                leg = legend(num2str((1 : 6).'), 'Location', 'best');
                title(leg, 'algorithm');
                title("Receiver operating characteristic (ROC) curve for \color{blue}SNR_{out} = " + num2str(options.snr_dB(1)) + " dB");
            else
                PFA = logspace(-10, 0, nop);
                if ~isscalar(options.snr_dB)
                    PD = zeros(nop, length(options.snr_dB));
                    for k = 1 : length(options.snr_dB)
                        PD(:, k) = obj.ROC(PFA, 10.^(.1*options.snr_dB(k)), obj.network.numberOfActiveReceivingNodes);
                    end
                    figure; semilogx(PFA, PD);
                    grid on; grid minor;
                    xlabel('PFA'); ylabel('PD');
                    leg = legend(num2str(options.snr_dB.'), 'Location', 'best');
                    switch obj.processingAlgorithm
                        case 6
                            title(leg, 'mean SNR_{out}');
                        otherwise
                            title(leg, 'SNR_{out}');
                    end
                    title("Receiver operating characteristic (ROC) curve of algorithm \color{blue}L_{" + num2str(obj.processingAlgorithm) + "}");
                else
                    PD = zeros(nop, length(options.numberOfReceivers));
                    for k = 1 : length(options.numberOfReceivers)
                        PD(:, k) = obj.ROC(PFA, 10.^(.1*options.snr_dB), options.numberOfReceivers(k));
                    end
                    figure; semilogx(PFA, PD);
                    grid on; grid minor;
                    xlabel('PFA'); ylabel('PD');
                    leg = legend(num2str(options.numberOfReceivers.', '%d'), 'Location', 'best');
                    title(leg, 'N_{RX}');
                    title("Receiver operating characteristic (ROC) curve of algorithm \color{blue}L_{" + num2str(obj.processingAlgorithm) + "}");
                end
            end
        end

        function visualizedetectioncharacteristics(obj, options)
            arguments
                obj
                options.PFA (1, :) double = [1e-4 1e-6 1e-8]
                options.numberOfReceivers (1, :) double = [1 : 5 10]
            end
            nop = 1000;
            SNR = linspace(1, 33, nop);
            if ~isscalar(options.PFA)
                PD = zeros(nop, length(options.PFA));
                for k = 1 : length(options.PFA)
                    PD(:, k) = obj.ROC(options.PFA(k), 10.^(.1*SNR), obj.network.numberOfActiveReceivingNodes);
                end
                legStr = num2str(options.PFA.', '%.2e');
                legTitle = 'PFA';
            else
                PD = zeros(nop, length(options.numberOfReceivers));
                for n = 1 : length(options.numberOfReceivers)
                    PD(:, n) = obj.ROC(options.PFA, 10.^(.1*SNR), options.numberOfReceivers(n));
                end
                legStr = num2str(options.numberOfReceivers.', '%d');
                legTitle = 'N_{RX}';
            end
            figure; plot(SNR, PD);
            grid on; grid minor;
            switch obj.processingAlgorithm
                case 6
                    xlabel('mean SNR_{out} (dB)');
                otherwise
                    xlabel('SNR_{out} (dB)');
            end
            ylabel('PD');
            leg = legend(legStr, 'Location', 'best');
            title(leg, legTitle);
            title("Detection characteristic curve of algorithm \color{blue}L_{" + num2str(obj.processingAlgorithm) + "}");
        end
    end

    methods (Static)
        function T = incompleteGammaInverse(PFA, N)
            switch N
                case 1
                    T = -log(PFA);
                otherwise
                    n = (0 : N - 1).';
                    T = zeros(1, length(PFA));
                    for pfaIdx = 1 : length(PFA)
                        if PFA(pfaIdx) == 1
                            T(pfaIdx) = 0;
                        else
                            thresholdSet = -log(PFA(pfaIdx)).*logspace(0, sqrt(N - 1), 1e4);
                            T(pfaIdx) = thresholdSet(islocalmin(abs(PFA(pfaIdx) - exp(-thresholdSet).*sum(thresholdSet.^n./factorial(n), 1))));
                        end
                    end
            end
        end
    end
end

% function Y = indexer(Y, loc, idx)
%     dimY = numel(size(Y));
%     expr = repmat(':,', 1, dimY);
%     expr = expr(1 : end - 1);
%     expr = [expr(1 : 2*loc - 2), 'idx', expr(2*loc : end)];
%     Y = eval(['Y(', expr, ')']);
% end