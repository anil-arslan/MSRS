classdef spu < handle
    %spu Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = private, GetAccess = public)
        interfaces (1, :) interface = interface.empty()
        monteCarlo (1, 1) struct = struct( ...
            'seed', 0, ...
            'numberOfTrials', 1, ...
            'numberOfTrialsParallel', 1)
        configuration (1, 1) struct = struct( ...
            'PFA', 1e-6, ...
            'PD', 1, ...
            'threshold', 0, ...
            'threshold_dB', -inf);
        configurationCompression (1, 1) struct = struct( ...
            'numberOfTargets', nan, ...
            'neighbourOffset', 1)
        configurationMonostatic (1, 1) struct = struct( ...
            'removeBoundaryDetectionDF', 1, ...
            'removeBoundaryDetectionMF', 0, ...
            'interpolationDF', 1, ...
            'interpolationMF', 0, ...
            'groupingCloseness', 2, ... 
            'CFAR', 1, ... % For coherent monostatic L2
            'numberOfTrainingCells', 20, ...
            'numberOfGuardCells', 5);
        gridResolution (3, 1) double {mustBeNonnegative} = 100*ones(3, 1)
        processingAlgorithm (1, 1) {mustBeInteger, mustBeInRange(processingAlgorithm, 1, 6)} = 1
        % 1: Coherent detector - Perfect knowledge of range
        %   Maximal ratio combining
        %   Deterministic signal
        % 2: Coherent detector - Fully correlated spatial observations
        %   Uniformly most powerful test (UMP)
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
        detectionAlgorithm (1, 1) string {mustBeMember(detectionAlgorithm, ["thresholding", "peak", "CoSaMP", "OMP", "OMPjoint"])} = "thresholding"
    end

    properties (Dependent)
        network (1, 1) radarNetwork
        gridPoints cell
        gridSize (1, 3) double {mustBeNonnegative}
        gridPointsMesh struct
        hypothesizedInterface interface
    end

    properties (Dependent)
        inputSNR_dB double % dB power, % (Ntx x Nrx x Nt x Nmcp matrix)
        inputSNR_lin dobule % linear scale power, % (Ntx x Nrx x Nt x Nmcp matrix)
        outputSNR_dB dobule % dB power, (1 x Nmonorx x Nt x Nmcp matrix)
        outputSNR_lin dobule % linear scale power, (1 x Nmonorx x Nt x Nmcp matrix)
        threshold function_handle
        ROC function_handle
        noisePowersPerSample_W double {mustBePositive} % linear scale noise power
        integrationWeights double % (Ntx x Nrx x Ni matrix)
    end

    properties (SetAccess = private, GetAccess = public)
        integrationIndices double {mustBeInteger, mustBePositive} % (Ntx x Nrx x Ni vector)
        expectedCellObservations double % (Ntx x Nrx x Ni matrix)
        blindZone logical {mustBeNumericOrLogical, mustBeMember(blindZone, [0, 1])} = []; 
    end

    properties (Dependent, Hidden)
        signalsMatchFilteredTrials cell % (1 x Nrx cell of Ns + L - 1 x Nrxch x Ntxch x Nmcp matrix) Ncmp : number of parallel trials
        signalsCompressedTrials double % (Ns_i x Ncmp matrix)
    end

    properties (SetAccess = private, GetAccess = public)
        signalsMatchFiltered cell = {0} % (1 x Nrx cell of Ns + L - 1 x Nrxch x Ntxch x Nmcp matrix) Ncmp : number of parallel trials
        dictionaryCompression = [] % (Ns_i x Ni x Nrx x Ntx matrix)
    end

    properties (Dependent)
        signalsIntegrated double % (Ni x 1 x Ncmp matrix)
        signalsIntegratedDirectly double % (Ns_i x Ncmp matrix)

        thresholdCFAR cell % (1 x Nrx cell of Ns + L - 1 x Nrxch x Nmcp matrix)
        detectionFromMatchFiltration struct % (Nrx x Nmcp struct of 1 x Nd vectors) Nd : number of detections

        hypothesisTestingResults cell % (1 x Nrx cell of Nrxch x Nmcp cell of Nd x 1 vector) Nd : number of detections
        estimatedPositions cell % (1 x Nrx cell of Nrxch x Nmcp cell of 3 x Nd matrix)
        detectionFromIntegration struct % (Nrx x Nmcp struct of 1 x Nd vectors) Nd : number of detections

        positionError cell % (1 x Nrx cell of Nrxch x Nmcp cell of 3 x Nd x Nt matrix)
        positionErrorTotal cell % (1 x Nrx cell of Nrxch x Nmcp cell of Nd x Nt matrix)

        covarianceMatrixInitial (3, 3) double
    end

    properties (SetAccess = private, GetAccess = public)
        compressionReport struct = struct( ...
            'cellIDs', [], ...
            'numberOfIterations', [], ...
            'integratedSignals', [], ...
            'complexAmplitudes', [], ...
            'residualPowerHistory', [], ...
            'maximumPowerHistory', []);
        coverageSimulationReport struct = struct( ...
            'targetCellIDs', [], ...
            'PD', [], ...
            'SNRs', [], ...
            'SNRsMean', [])
    end

    properties (Access = private)
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
            % (Ntx x Nrx x Nt x Nmcp matrix)
            SNRin = 10.^(.1*obj.inputSNR_dB);
        end

        function SNRin = get.inputSNR_dB(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            SNRin = obj.interfaces.inputSNR_dB;
        end

        function SNRout = get.outputSNR_lin(obj)
            % (1 x Nmonorx x Nt x Nmcp matrix)
            switch obj.network.networkMode
                case "multiStatic"
                    SNRout = sum(10.^(.1*(obj.network.processingGain_dB.' + obj.inputSNR_dB)), [1 2]); % (1 x 1 x Nt x Nmcp matrix)
                case "monoStatic"
                    SNRout = zeros([1 size(obj.inputSNR_dB, [2 3 4])]);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        txID = obj.network.monoStaticTransmitterIDs(rxID);
                        SNRout(:, rxID, :, :) = 10.^(.1*(obj.network.processingGain_dB(rxID, txID) + obj.inputSNR_dB(txID, rxID, :, :)));
                    end
            end
        end

        function SNRout = get.outputSNR_dB(obj)
            % (1 x Nmonorx x Nt x Nmcp matrix)
            SNRout = 10*log10(obj.outputSNR_lin);
        end

        function PD = get.ROC(obj)
            if obj.network.passiveNetwork
                % Passive networks
                error('not implemented');
            else
                switch obj.processingAlgorithm
                    case 1 % Coherent detector - Perfect knowledge of range
                        %   Maximal ratio combining
                        %   Deterministic signal

                        % EXACTLY CORRECT FORMULA
                        PD = @(PFA, SNR, N) erfc(erfcinv(2*PFA) - sqrt(SNR))/2;

                    case 2 % Coherent detector - Fully correlated spatial observations
                        %   Uniformly most powerful test (UMP)
                        %   Spatially fully correlated phase and amplitude
                        %   For weak stochastic signal, it is optimal for any distribution
                        %   Otherwise, it is optimal under AWGN
                        %   For generalized LRT, it is optimal adaptive processing algorithm

                        % EXACTLY CORRECT FORMULA
                        switch obj.interfaces.swerling
                            case 0 % only absolute phase fluctuates
                                % If the received signal has the same spatial
                                % probability distributions
                                PD = @(PFA, SNR, N) marcumq(sqrt(2*SNR), sqrt(-2*log(PFA)));
                            case 1 % both absolute phase and absolute amplitude fluctuates
                                % If the received signal has different but fully
                                % correlated spatial probability distributions
                                PD = @(PFA, SNR, N) PFA.^(1./(1 + SNR));
                        end

                    case 3 % Noncoherent detector - Square-law envelope detector
                        %   Weak stochastic signal
                        %   Spatially independent phase, arbitrarily correlated amplitude
                        %   For weak stochastic signal, it is optimal for any distribution
                        %   Otherwise, it is an approximation of optimal processing algorithm under AWGN

                        PD = @(PFA, SNR, N) nan;
                        % not implemented
                        switch obj.network.networkMode
                            case "multiStatic"
                                switch obj.network.numberOfActiveBistaticPairs
                                    case 1
                                        PD = @(PFA, SNR, N) PFA.^(1./(1 + SNR)); % same as swerling 1
                                    otherwise % check again
                                        PD = @(PFA, SNR, N) exp(-obj.threshold(PFA, N)./(1 + SNR/N)).*sum((obj.threshold(PFA, N)./(1 + SNR/N)).^(shiftdim(0 : N - 1, -1))./factorial(shiftdim(0 : N - 1, -1)), 3);
                                end
                            case "monoStatic"
                                PD = @(PFA, SNR, N) nan;
                        end

                    case 4 % Noncoherent detector - Linear-law envelope detector
                        %   Strong stochastic signal
                        %   Spatially independent phase, fully correlated amplitude
                        %   For generalized LRT, it is an approximation of optimal adaptive processing algorithm

                        PD = @(PFA, SNR, N) nan;
                        % not implemented

                    case 5 % Noncoherent detector - Square-law envelope detector
                        %   Neither strong nor weak stochastic signal
                        %   Spatially independent phase and amplitude

                        PD = @(PFA, SNR, N) nan;
                        % not implemented
                        switch obj.network.networkMode
                            case "multiStatic"
                                switch obj.network.numberOfActiveBistaticPairs
                                    case 1
                                        PD = @(PFA, SNR, N) PFA.^(1./(1 + SNR)); % same as swerling 1
                                    otherwise % check again
                                        PD = @(PFA, SNR, N) exp(-obj.threshold(PFA, N)./(1 + SNR/N)).*sum((obj.threshold(PFA, N)./(1 + SNR/N)).^(shiftdim(0 : N - 1, -3))./factorial(shiftdim(0 : N - 1, -3)), 5);
                                end
                            case "monoStatic"
                                PD = @(PFA, SNR, N) nan;
                        end

                    case 6 % Noncoherent detector - Square-law envelope detector
                        %   Strong stochastic signal
                        %   Spatially independent phase and amplitude
                        %   For generalized LRT, it is optimal adaptive processing algorithm
                        switch obj.network.networkMode
                            case "multiStatic"
                                switch obj.network.numberOfActiveBistaticPairs
                                    case 1
                                        PD = @(PFA, SNR, N) PFA.^(1./(1 + SNR)); % same as swerling 1
                                    otherwise % check again
                                        PD = @(PFA, SNR, N) exp(-obj.threshold(PFA, N)./(1 + SNR/N)).*sum((obj.threshold(PFA, N)./(1 + SNR/N)).^(shiftdim(0 : N - 1, -3))./factorial(shiftdim(0 : N - 1, -3)), 5);
                                end
                            case "monoStatic"
                                PD = @(PFA, SNR, N) nan;
                        end
                end
            end
        end

        function T = get.threshold(obj)
            if obj.network.passiveNetwork
                % Passive networks
                error('not implemented');
            else
                switch obj.processingAlgorithm
                    case 1 % Coherent detector - Perfect knowledge of range
                        %   Maximal ratio combining
                        %   Deterministic signal

                        % EXACTLY CORRECT FORMULA
                        T = @(PFA, N) erfcinv(2*PFA).^2;

                    case 2 % Coherent detector - Fully correlated spatial observations
                        %   Spatially fully correlated phase and amplitude
                        %   For weak stochastic signal, it is optimal for any distribution
                        %   Otherwise, it is optimal under AWGN
                        %   For generalized LRT, it is optimal adaptive processing algorithm

                        % EXACTLY CORRECT FORMULA
                        T = @(PFA, N) -log(PFA);

                    case 3 % Noncoherent detector - Square-law envelope detector
                        %   Weak stochastic signal
                        %   Spatially independent phase, arbitrarily correlated amplitude
                        %   For weak stochastic signal, it is optimal for any distribution
                        %   Otherwise, it is an approximation of optimal processing algorithm under AWGN

                        T = @(PFA, N) gammaincinv(1 - PFA, N);

                    case 4 % Noncoherent detector - Linear-law envelope detector
                        %   Strong stochastic signal
                        %   Spatially independent phase, fully correlated amplitude
                        %   For generalized LRT, it is an approximation of optimal adaptive processing algorithm

                        T = @(PFA, N) nan;

                    case 5 % Noncoherent detector - Square-law envelope detector
                        %   Neither strong nor weak stochastic signal
                        %   Spatially independent phase and amplitude

                        % not implemented/for same SNR formula changes
                        T = @(PFA, N) gammaincinv(1 - PFA, N);

                    case 6 % Noncoherent detector - Square-law envelope detector
                        %   Strong stochastic signal
                        %   Spatially independent phase and amplitude
                        %   For generalized LRT, it is optimal adaptive processing algorithm

                        % check again
                        T = @(PFA, N) gammaincinv(PFA, N, 'upper');

                end
            end
        end

        function cfg = get.configuration(obj)
            cfg = obj.configuration;
            switch obj.network.networkMode
                case "multiStatic"
                    cfg.PD = obj.ROC(cfg.PFA, obj.outputSNR_lin, obj.network.numberOfActiveBistaticPairs);
                    cfg.threshold = obj.threshold(cfg.PFA, obj.network.numberOfActiveBistaticPairs);
                case "monoStatic"
                    cfg.PD = obj.ROC(cfg.PFA, obj.outputSNR_lin, 1);
                    cfg.threshold = obj.threshold(cfg.PFA, 1);
            end
            cfg.threshold_dB = 10*log10(abs(cfg.threshold));
        end

        function n = get.noisePowersPerSample_W(obj)
            n = 10.^(.1*[obj.network.activeReceivingNodes.noisePowerPerSample_dB]); % kTBN
        end

        function w = get.integrationWeights(obj)
            % (Ntx x Nrx x Ni matrix)
            switch obj.processingAlgorithm
                case 1 % Coherent detector - Perfect knowledge of range
                    %   Maximal ratio combining
                    %   Deterministic signal

                    % w = obj.expectedCellObservations;
                    % not implemented
                    w = ones(obj.network.numberOfActiveTransmittingNodes, obj.network.numberOfActiveReceivingNodes, prod(obj.gridSize));

                case 2 % Coherent detector - Fully correlated spatial observations
                    %   Spatially fully correlated phase and amplitude
                    %   For weak stochastic signal, it is optimal for any distribution
                    %   Otherwise, it is optimal under AWGN
                    %   For generalized LRT, it is optimal adaptive processing algorithm

                    % checked/ not implemented
                    w = obj.expectedCellObservations;

                case 3 % Noncoherent detector - Square-law envelope detector
                    %   Weak stochastic signal
                    %   Spatially independent phase, arbitrarily correlated amplitude
                    %   For weak stochastic signal, it is optimal for any distribution
                    %   Otherwise, it is an approximation of optimal processing algorithm under AWGN

                    % not checked
                    w = abs(obj.expectedCellObservations).^2;
                    w = w./vecnorm(vecnorm(w, 2, 1), 2, 2).*obj.network.numberOfActiveBistaticPairs;

                case 4 % Noncoherent detector - Linear-law envelope detector
                    %   Strong stochastic signal
                    %   Spatially independent phase, fully correlated amplitude
                    %   For generalized LRT, it is an approximation of optimal adaptive processing algorithm

                    % not checked
                    w = abs(obj.expectedCellObservations);

                case 5 % Noncoherent detector - Square-law envelope detector
                    %   Neither strong nor weak stochastic signal
                    %   Spatially independent phase and amplitude

                    % not checked
                    w = 1./(1 + 1./abs(obj.expectedCellObservations).^2);
                    w = w./vecnorm(vecnorm(w, 2, 1), 2, 2).*obj.network.numberOfActiveBistaticPairs;

                case 6 % Noncoherent detector - Square-law envelope detector
                    %   Strong stochastic signal
                    %   Spatially independent phase and amplitude
                    %   For generalized LRT, it is optimal adaptive processing algorithm

                    w = ones(obj.network.numberOfActiveTransmittingNodes, obj.network.numberOfActiveReceivingNodes, prod(obj.gridSize));
            end
            switch obj.network.networkCoherency
                case "short-term coherent"
                    if any(obj.processingAlgorithm, [1 2])
                        error('Processing algorithm %d cannot be used for short-term coherent networks', obj.processingAlgorithm);
                    end
                case "incoherent"
                    if any(obj.processingAlgorithm, [1 2 4])
                        error('Processing algorithm %d cannot be used for incoherent networks', obj.processingAlgorithm);
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

        function y = get.signalsMatchFilteredTrials(obj)
            % (1 x Nrx cell of Ns + L - 1 x Nrxch x Ntxch Ncmp matrix) Ncmp : number of parallel trials
            mc = obj.monteCarlo;
            matchFilters = obj.network.matchFilter; % L x Nrx x Ntx matrix
            %%% carrier demodulation not implemented
            % t = [obj.network.activeReceivingNodes.samplingInstants]; % Ns x Nrx matrix
            % demodulator = exp(1j*2*pi*shiftdim([obj.network.activeTransmittingNodes.carrierFrequency], -1).*t);
            Ns = size([obj.network.activeReceivingNodes.samplingInstants], 1);
            signalsBeamformed = obj.interfaces.signalBeamformed;
            y = cell(1, obj.network.numberOfActiveReceivingNodes);
            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
                switch obj.network.networkMode
                    case 'multiStatic'
                        numberOfTransmitters = obj.network.numberOfActiveTransmittingNodes;
                    case 'monoStatic'
                        numberOfTransmitters = 1;
                end
                y{rxID} = zeros(Ns + size(matchFilters, 1) - 1, numberOfTotalChannels, numberOfTransmitters, mc.numberOfTrialsParallel);
                s = signalsBeamformed{rxID}; % Ns x Nrxch x Ntxch x Nmcp matrix
                for txID = 1 : numberOfTransmitters
                    switch obj.network.networkMode
                        case 'multiStatic'
                            matchFilter = conj(matchFilters(:, rxID, txID));
                        case 'monoStatic'
                            matchFilter = conj(matchFilters(:, rxID, obj.network.monoStaticTransmitterIDs(rxID)));
                    end
                    for chID = 1 : numberOfTotalChannels
                        for mcID = 1 : mc.numberOfTrialsParallel
                            signalChannel = s(:, chID, txID, mcID);
                            signalChannel(isinf(signalChannel)) = 0;
                            y{rxID}(:, chID, txID, mcID) = conv(signalChannel, matchFilter);
                        end
                    end
                end
                y{rxID} = y{rxID}./sqrt(obj.noisePowersPerSample_W(rxID));
            end
        end

        function report = get.signalsCompressedTrials(obj)
            % (Ni x Nmcp matrix)
            if ~any(strcmpi(obj.detectionAlgorithm, ["OMP", "OMPjoint", "CoSaMP"]))
                error('Select a compression algorithm. Current algorithm is "%s"', obj.detectionAlgorithm);
            end
            mc = obj.monteCarlo;
            config = obj.configuration;
            %%% carrier demodulation not implemented
            % t = [obj.network.activeReceivingNodes.samplingInstants]; % Ns x Nrx matrix
            % demodulator = exp(1j*2*pi*shiftdim([obj.network.activeTransmittingNodes.carrierFrequency], -1).*t);
            Ns = size([obj.network.activeReceivingNodes.samplingInstants], 1);
            signalBeamformed = obj.interfaces.signalBeamformed;
            Nfus = obj.network.numberOfParallelSignalFusions;
            report = struct( ...
                'cellIDs', [], ...
                'numberOfIterations', [], ...
                'integratedSignals', [], ...
                'complexAmplitudes', [], ...
                'residualPowerHistory', [], ...
                'maximumPowerHistory', []);
            report = repmat(report, Nfus, mc.numberOfTrialsParallel);
            for fusionID = 1 : Nfus
                switch obj.network.networkMode
                    case 'multiStatic'
                        Nrxch = obj.network.numberOfActiveReceivingNodes;
                        Ntxrx = obj.network.numberOfActiveTransmittingNodes;
                        measurements = zeros(Ns, Nrxch, Ntxrx, mc.numberOfTrialsParallel);
                        for rxID = 1 : Nrxch
                            s = signalBeamformed{rxID}; % Ns x Nrxch x Ntxch x Nmcp matrix
                            %%% multi channel not implemented
                            measurements(:, rxID, :, :) = s(:, ceil(end/2), :, :)./sqrt(obj.noisePowersPerSample_W(rxID));
                        end
                        dict = obj.dictionaryCompression; % (Ns_i x Ni x Nrxch x Ntxch matrix)
                        w = obj.integrationWeights; % (Ntx x Nrx x Ni matrix)
                        w = permute(w, [3 4 2 1]);
                    case 'monoStatic'
                        Nrxch = obj.network.activeReceivingNodes(fusionID).numberOfTotalChannels;
                        Ntxrx = 1;
                        measurements = signalBeamformed{fusionID}./sqrt(obj.noisePowersPerSample_W(fusionID)); % Ns x Nrxch x 1 x Nmcp matrix
                        dict = obj.dictionaryCompression{fusionID}; % (Ns_i x Ni x Nrxch matrix)
                        w = 1;
                end
                Nc = size(dict, 2); % number of cells
                if isnan(obj.configurationCompression.numberOfTargets)
                    targetResidual = config.threshold;
                    numberOfIterations = Ns;
                else
                    targetResidual = -Inf;
                    numberOfIterations = obj.configurationCompression.numberOfTargets;
                    if numberOfIterations > Ns
                        error('Number of targets cannot be larger than the number of samples');
                    end
                end
                switch obj.detectionAlgorithm
                    case "OMPjoint"
                        numOfParallelOpts = Nrxch*Ntxrx;
                        numOfSamples = Ns;
                    case "OMP"
                        % check
                        numOfParallelOpts = 1;
                        numOfSamples = Ns*Nrxch*Ntxrx;
                        dict = reshape(permute(dict, [1 3 4 2]), numOfSamples, []); % Ns_i*Nrxch*Ntxch x Ni
                        measurements = reshape(measurements, numOfSamples, 1, 1, []); % Ns_i*Nrxch*Ntxch x 1 x 1 x Nmcp
                    case "CoSaMP"
                        error('not implemented');
                        %%% not implemented
                        % report = CoSaMP(dict, measurements(:), obj.configurationCompression.numberOfTargets);
                end
                for mcID = 1 : mc.numberOfTrialsParallel
                    orthAtomSet = zeros(numOfSamples, numberOfIterations, numOfParallelOpts);
                    nonorthAtomSet = zeros(numOfSamples, numberOfIterations, numOfParallelOpts);
                    cellIDset = zeros(numberOfIterations, 1);
                    residualPowerHistory = zeros(numberOfIterations + 1, 1);
                    maximumPowerHistory = zeros(numberOfIterations + 1, 1);
                    if obj.saveResiduals
                        report(fusionID, mcID).integratedSignals = nan(Nc, numberOfIterations + 1);
                    end
                    thresholdReached = false;
                    residual = permute(measurements(:, :, :, mcID), [1 4 2 3]); % (Ns x 1 x Nrxch x Ntxch) or (Ns*Nrxch*Ntxch x 1)
                    for currentIterationID = 1 : numberOfIterations
                        switch obj.processingAlgorithm
                            case 1
                                integratedSignal = real(sum(w.*pagemtimes(dict, 'ctranspose', residual, 'none'), [3 4])); % coherent integration
                                integratedSignal = sign(integratedSignal).*integratedSignal.^2;
                            case 2
                                integratedSignal = abs(sum(w.*pagemtimes(dict, 'ctranspose', residual, 'none'), [3 4])).^2; % coherent integration
                                switch obj.network.networkMode
                                    case 'monoStatic'
                                        switch obj.network.activeReceivingNodes(fusionID).beamformingMode
                                            case 'bypass'
                                                integratedSignal = integratedSignal./Nrxch;
                                            otherwise
                                        end
                                end
                            case {3, 5}
                                integratedSignal = sum(w.*abs(pagemtimes(dict, 'ctranspose', residual, 'none')).^2, [3 4]); % noncoherent integration
                            case 6
                                integratedSignal = sum(abs(pagemtimes(dict, 'ctranspose', residual, 'none')).^2, [3 4]); % noncoherent integration
                        end
                        [maximumPowerHistory(currentIterationID), cellID] = max(integratedSignal, [], 1);
                        residualPowerHistory(currentIterationID) = sum(abs(residual).^2, 'all');
                        if obj.saveResiduals
                            report(fusionID, mcID).integratedSignals(:, currentIterationID) = integratedSignal; % Ni x 1
                        end
                        % if currentIterationID == 1 && mcID == 1
                        %     fprintf('Max power = %g dB, Threshold = %g dB\n', 10*log10(maxPower), 10*log10(config.threshold));
                        % end
                        % residual power max powerdan daha yuksek geliyor.
                        % if residualPowerHistory(currentIterationID) < targetResidual || (currentIterationID > 1 && residualPowerHistory(currentIterationID) > residualPowerHistory(currentIterationID - 1))
                        %     thresholdReached = true;
                        %     break;
                        % end
                        if maximumPowerHistory(currentIterationID) < targetResidual
                            thresholdReached = true;
                            break;
                        end
                        if ismember(cellID, cellIDset(1 : currentIterationID - 1))
                            warning('same cellID with previous iteration');
                            thresholdReached = true;
                            break;
                        end
                        cellIDset(currentIterationID) = cellID;
                        % dict: (Ns_i*Nrxch*Ntxch x Ni) or (Ns_i x Ni x Nrxch x Ntxch matrix)
                        currentAtoms = reshape(dict(:, cellID, :, :), numOfSamples, numOfParallelOpts);
                        atomNorms = sqrt(sum(abs(currentAtoms).^2, 1));
                        nonZeroAtoms = atomNorms ~= 0;
                        currentAtoms(:, nonZeroAtoms) = currentAtoms(:, nonZeroAtoms)./atomNorms(nonZeroAtoms);
                        nonorthAtomSet(:, currentIterationID, :) = currentAtoms; 
                
                        % -- Step 2: update residual
                    
                        % First, orthogonalize 'atoms' against all previous atoms
                        % We use Modified Gram Schmidt
                        for previousIterationID = 1 : (currentIterationID - 1)
                            for chID = 1 : numOfParallelOpts
                                currentAtoms(:, chID) = currentAtoms(:, chID) - (orthAtomSet(:, previousIterationID, chID)'*currentAtoms(:, chID))*orthAtomSet(:, previousIterationID, chID);
                            end
                        end
                        % Second, normalize:
                        atomNorms = sqrt(sum(abs(currentAtoms).^2, 1));
                        nonZeroAtoms = atomNorms ~= 0;
                        currentAtoms(:, nonZeroAtoms) = currentAtoms(:, nonZeroAtoms)./atomNorms(nonZeroAtoms);
                        orthAtomSet(:, currentIterationID, :) = currentAtoms;
                        % Third, solve least-squares problem
                        for chID = 1 : numOfParallelOpts
                            % Fourth, update residual:
                            residual(:, chID) = residual(:, chID) - orthAtomSet(:, 1 : currentIterationID, chID)*(orthAtomSet(:, 1 : currentIterationID, chID)'*residual(:, chID));
                        end
                    end
                    if thresholdReached
                        currentIterationID = currentIterationID - 1;
                    end
                    %%% iterations end
                    % For the last iteration, we need to do this without orthogonalizing dictionary
                    % so that the complexAmplitudes match what is expected.
                    if currentIterationID
                        report(fusionID, mcID).complexAmplitudes = zeros(currentIterationID, numOfParallelOpts);
                        for chID = 1 : numOfParallelOpts
                            if all(any(nonorthAtomSet(:, 1 : currentIterationID, chID)))
                                measurement = measurements(:, :, :, mcID);
                                report(fusionID, mcID).complexAmplitudes(:, chID) = pinv(nonorthAtomSet(:, 1 : currentIterationID, chID))*measurement(:, chID);
                            else
                                report(fusionID, mcID).complexAmplitudes(:, chID) = 0;
                            end
                        end
                    else
                        report(fusionID, mcID).complexAmplitudes = zeros(1, numOfParallelOpts);
                    end
                    visibleCellIDs = false(Nc, 1);
                    visibleCellIDs(cellIDset(1 : currentIterationID)) = true;
                    gridCellIDs = false(prod(obj.gridSize), 1);
                    gridCellIDs(~obj.blindZone) = visibleCellIDs;
                    report(fusionID, mcID).cellIDs = find(gridCellIDs);
                    report(fusionID, mcID).residualPowerHistory = residualPowerHistory(1 : currentIterationID + 1);
                    report(fusionID, mcID).maximumPowerHistory = maximumPowerHistory(1 : currentIterationID + 1);
                    if obj.saveResiduals
                        report(fusionID, mcID).integratedSignals = report(fusionID, mcID).integratedSignals(:, 1 : currentIterationID + 1) ;
                    end
                    report(fusionID, mcID).numberOfIterations = currentIterationID + 1;
                end
            end
        end

        function Y = get.signalsIntegrated(obj)
            % (Ni x Nmcp matrix)
            if isempty(obj.integrationIndices)
                fprintf('integration indices are not set\n');
                Y = [];
                return;
            end
            if isscalar(obj.signalsMatchFiltered{1})
                obj.setmatchfilteredsignals;
            end
            mc = obj.monteCarlo;
            visibleZone = ~obj.blindZone;
            switch obj.network.networkMode
                case "multiStatic"
                    w = permute(obj.integrationWeights, [3 1 2]);
                    if obj.network.passiveNetwork % Passive networks
                    else 
                        Y = zeros(prod(obj.gridSize), 1, 1, mc.numberOfTrialsParallel);
                        switch obj.processingAlgorithm
                            case 1 % Coherent detector, deterministic signal
                                for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                                        %%% which channel not implemented
                                        Y(visibleZone, :, :, :) = Y(visibleZone, :, :, :) + w(:, txID, rxID).*obj.signalsMatchFiltered{rxID}(obj.integrationIndices(txID, rxID, :), ceil(end/2), txID, :);
                                    end
                                end
                                Y = real(Y);
                                Y = sign(Y).*Y.^2;
                            case 2 % Coherent detector, fully correlated stochastic signal
                                for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                                        %%% which channel not implemented
                                        Y(visibleZone, :, :, :) = Y(visibleZone, :, :, :) + w(:, txID, rxID).*obj.signalsMatchFiltered{rxID}(obj.integrationIndices(txID, rxID, :), ceil(end/2), txID, :);
                                    end
                                end
                                Y = abs(Y).^2;
                            case {3, 5} % Square-law envelope detectors
                                for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                                        %%% which channel not implemented
                                        Y(visibleZone, :, :, :) = Y(visibleZone, :, :, :) + w(:, txID, rxID).*abs(obj.signalsMatchFiltered{rxID}(obj.integrationIndices(txID, rxID, :), ceil(end/2), txID, :)).^2; 
                                    end
                                end
                            case 4 % Linear-law envelope detector
                                for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                                        %%% which channel not implemented
                                        Y(visibleZone, :, :, :) = Y(visibleZone, :, :, :) + w(:, txID, rxID).*abs(obj.signalsMatchFiltered{rxID}(obj.integrationIndices(txID, rxID, :), ceil(end/2), txID, :));
                                    end
                                end
                            case 6 % Square-law envelope detectors
                                for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                                        %%% which channel not implemented
                                        Y(visibleZone, :, :, :) = Y(visibleZone, :, :, :) + abs(obj.signalsMatchFiltered{rxID}(obj.integrationIndices(txID, rxID, :), ceil(end/2), txID, :)).^2; 
                                    end
                                end
                        end
                        Y = permute(Y, [1 2 4 3]);
                    end
                case "monoStatic"
                    Y = cell(1, obj.network.numberOfActiveReceivingNodes);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
                        monoStaticTXIDs = obj.network.monoStaticTransmitterIDs(rxID);
                        y = zeros(prod(obj.gridSize), numberOfTotalChannels, 1, mc.numberOfTrialsParallel);
                        y(visibleZone, :, :, :) = obj.signalsMatchFiltered{rxID}(obj.integrationIndices(monoStaticTXIDs, rxID, :), :, :, :);
                        switch obj.processingAlgorithm
                            case 1 % Coherent detector, deterministic signal
                                Y{rxID} = real(y);
                                Y{rxID} = sign(Y{rxID}).*Y{rxID}.^2;
                            otherwise
                                Y{rxID} = abs(y).^2;
                        end
                        Y{rxID} = permute(Y{rxID}, [1 2 4 3]);
                    end
            end
        end

        function Y = get.signalsIntegratedDirectly(obj)
            if isempty(obj.signalsMatchFiltered)
                Y = [];
                fprintf('"signalsMatchFiltered" is empty\n');
                return;
            end
            if isscalar(obj.signalsMatchFiltered{1})
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
                                    Y = Y + shiftdim(abs(obj.signalsMatchFiltered{rxID}(:, ceil(end/2), :, :)).^2, 1 - rxID); 
                                end
                            end
                        case 4 % Linear-law envelope detector
                    end
                case "monoStatic"
                    Y = nan;
            end
        end

        function T = get.thresholdCFAR(obj)
            % (1 x Nrx cell of Ns + L - 1 x Nrxch x Nmcp matrix)
            T = cell(1, obj.network.numberOfActiveReceivingNodes);
            config = obj.configuration;
            mc = obj.monteCarlo;
            NTC = obj.configurationMonostatic.numberOfTrainingCells;
            NGC = obj.configurationMonostatic.numberOfGuardCells;
            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
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
                    T{rxID}(sampleID, :, :) = alpha.*mean(abs(obj.signalsMatchFiltered{rxID}(sampleIdx, :, :, :)).^2);
                end
            end
        end

        function detection = get.detectionFromMatchFiltration(obj)
            % (Nrx x Nmcp cell of 1 x Nd vector)
            detection = struct( ...
                'power', [], ...
                'range', [], ...
                'elevation', [], ...
                'azimuth', [], ...
                'numberOfDetections', [], ...
                'position', []);
            if isempty(obj.signalsMatchFiltered)
                fprintf('"signalsMatchFiltered" is empty\n');
                return;
            end
            if isscalar(obj.signalsMatchFiltered{1})
                obj.setmatchfilteredsignals;
            end
            config = obj.configuration;
            mc = obj.monteCarlo;
            detection = repmat(detection, obj.network.numberOfActiveReceivingNodes, mc.numberOfTrialsParallel);
            switch obj.network.networkMode
                case "monoStatic"
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        Ts = obj.network.activeReceivingNodes(rxID).samplingPeriod;
                        N = obj.network.activeReceivingNodes(rxID).numberOfSamplesPerCPI;
                        numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
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
                        T = config.threshold*ones(Ns*numberOfTotalChannels, mc.numberOfTrialsParallel);
                        obs = obj.signalsMatchFiltered{rxID};
                        absObs = reshape(abs(obs), [], mc.numberOfTrialsParallel);
                        switch obj.processingAlgorithm
                            case 1 % Coherent detector, deterministic signal
                                thObs = real(obj.signalsMatchFiltered{rxID});
                                thObs = sign(thObs).*thObs.^2;
                            otherwise % Square-law envelope detectors
                                if obj.configurationMonostatic.CFAR
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
                            if obj.configurationMonostatic.removeBoundaryDetectionMF
                                boundaryDetections = detectionSamples ~= 1 & detectionSamples ~= Ns;
                                detectionIndex = detectionIndex(boundaryDetections);
                                detectionSamples = detectionSamples(boundaryDetections);
                            end
                            numberOfDetections = length(detectionIndex);
                            if numberOfDetections
                                detectionPower = 20*log10(absObs(detectionIndex, mcID));
                                idx = 1;
                                while idx < numberOfDetections
                                    group = abs(detectionSamples - detectionSamples(idx)) < obj.configurationMonostatic.groupingCloseness;
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
                                if obj.configurationMonostatic.interpolationMF
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
                                if obj.configurationMonostatic.removeBoundaryDetectionDF
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
                                    if obj.configurationMonostatic.interpolationDF
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
            end
        end

        function idx = get.hypothesisTestingResults(obj)
            switch obj.detectionAlgorithm
                case {"thresholding", "peak"}
                    if isempty(obj.signalsMatchFiltered)
                        idx = {};
                        fprintf('"signalsMatchFiltered" is empty\n');
                        return;
                    end
                    if isscalar(obj.signalsMatchFiltered{1})
                        obj.setmatchfilteredsignals;
                    end
                case {"OMPjoint", "OMP"}
                    if isempty(obj.compressionReport) || isempty(obj.compressionReport(1).complexAmplitudes)
                        obj.applycompression;
                    end
                case "CoSaMP"
                    error("not implemented");
            end
            config = obj.configuration;
            mc = obj.monteCarlo;
            switch obj.network.networkMode
                case "multiStatic"
                    idx = cell(1, mc.numberOfTrialsParallel);
                    for mcID = 1 : mc.numberOfTrialsParallel
                        switch obj.detectionAlgorithm
                            case "thresholding"
                                idx{mcID} = find(obj.signalsIntegrated(:, :, mcID) > config.threshold);
                            case "peak"
                                z = obj.signalsIntegrated;
                                idx{mcID} = find(z(:, :, mcID) > config.threshold);
                                [~, peakIndex] = max(z(idx{mcID}, :, mcID));
                                idx{mcID} = idx{mcID}(peakIndex);
                            case {"OMPjoint", "OMP"}
                                idx{mcID} = obj.compressionReport(1, mcID).cellIDs;
                            case "CoSaMP"
                                error("not implemented");
                        end
                    end
                    if ~obj.network.passiveNetwork % Single transmitter networks
                    else % Passive networks
                    end
                case "monoStatic"
                    idx = cell(1, obj.network.numberOfActiveReceivingNodes);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        T = config.threshold;
                        numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
                        switch obj.detectionAlgorithm
                            case "thresholding"
                                z = obj.signalsIntegrated;
                                idx{rxID} = cell(numberOfTotalChannels, mc.numberOfTrialsParallel);
                                for mcID = 1 : mc.numberOfTrialsParallel
                                    for chID = 1 : numberOfTotalChannels
                                        idx{rxID}{chID, mcID} = find(z{rxID}(:, chID, mcID) > T);
                                    end
                                end
                            case "peak"
                                z = obj.signalsIntegrated;
                                idx{rxID} = cell(numberOfTotalChannels, mc.numberOfTrialsParallel);
                                for mcID = 1 : mc.numberOfTrialsParallel
                                    for chID = 1 : numberOfTotalChannels
                                        idx{rxID}{chID, mcID} = find(z{rxID}(:, chID, mcID) > T);
                                        [~, peakIndex] = max(z{rxID}(idx{rxID}{chID, mcID}, chID, mcID));
                                        idx{rxID}{chID, mcID} = idx{rxID}{chID, mcID}(peakIndex);
                                    end
                                end
                            case {"OMPjoint", "OMP"}
                                idx{rxID} = cell(1, mc.numberOfTrialsParallel);
                                for mcID = 1 : mc.numberOfTrialsParallel
                                    idx{rxID}{mcID} = obj.compressionReport(rxID, mcID).cellIDs;
                                end
                            case "CoSaMP"
                                error("not implemented");
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
                        switch obj.detectionAlgorithm
                            case {"thresholding", "peak"}
                                numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
                                pos{rxID} = cell(numberOfTotalChannels, mc.numberOfTrialsParallel);
                                for mcID = 1 : mc.numberOfTrialsParallel
                                    for chID = 1 : numberOfTotalChannels
                                        pos{rxID}{chID, mcID} = [gridScan.x(idx{rxID}{chID, mcID}) gridScan.y(idx{rxID}{chID, mcID}) gridScan.z(idx{rxID}{chID, mcID})].';
                                    end
                                end
                            case {"OMPjoint", "OMP"}
                                pos{rxID} = cell(1, mc.numberOfTrialsParallel);
                                for mcID = 1 : mc.numberOfTrialsParallel
                                    pos{rxID}{mcID} = [gridScan.x(idx{rxID}{mcID}) gridScan.y(idx{rxID}{mcID}) gridScan.z(idx{rxID}{mcID})].';
                                end
                        end
                    end
            end
        end

        function detection = get.detectionFromIntegration(obj)
            % (Nrx x Nmcp cell of 1 x Nd vector)
            detection = struct( ...
                'power', [], ...
                'x', [], ...
                'y', [], ...
                'z', [], ...
                'numberOfDetections', [], ...
                'position', []);
            if isempty(obj.integrationIndices)
                fprintf('integration indices are not set\n');
                return;
            end
            if isscalar(obj.signalsMatchFiltered{1})
                obj.setmatchfilteredsignals;
            end
            mc = obj.monteCarlo;
            gridScan = obj.gridPointsMesh;
            detection = repmat(detection, 1, mc.numberOfTrialsParallel);
            Y = obj.signalsIntegrated;
            if isempty(Y)
                fprintf('"signalsIntegrated" is empty\n');
                return;
            end
            idx = obj.hypothesisTestingResults;
            switch obj.network.networkMode
                case 'multiStatic'
                    for mcID = 1 : mc.numberOfTrialsParallel
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
                case 'monoStatic'
                    detection = repmat(detection, obj.network.numberOfActiveReceivingNodes, 1);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        for chID = 1 : obj.network.receivingNodes(rxID).numberOfTotalChannels
                            for mcID = 1 : mc.numberOfTrialsParallel
                                index = idx{rxID}{chID, mcID};
                                x = gridScan.x(index);
                                y = gridScan.y(index);
                                z = gridScan.z(index);
                                position = [x y z].';
                                detection(mcID).power{chID} = Y{rxID}(index, chID, mcID);
                                detection(mcID).x{chID} = x;
                                detection(mcID).y{chID} = y;
                                detection(mcID).z{chID} = z;
                                detection(mcID).numberOfDetections{chID} = numel(index);
                                detection(mcID).position{chID} = position;
                            end
                        end
                    end
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
                        numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
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
                        numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
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

        function int = get.hypothesizedInterface(obj)
            gridScan = obj.gridPointsMesh;
            hypothesizedPositions = reshape(permute(cat(4, gridScan.x, gridScan.y, gridScan.z), [4 1 2 3]), [3 prod(obj.gridSize)]);
            targets = target( ...
                "position", hypothesizedPositions, ...
                "meanRCS_dbsm", 0);
            int = interface( ...
                'network', obj.network, ...
                'targets', targets);
            int.configure( ...
                'noise', 0, ...
                'directPath', 0, ...
                'pathLoss', 0, ...
                'spatialCoherency', 'deterministic');
        end

        %%% set methods

        function setblindzone(obj)
            %%% add back of array to blind zone
            interfaceGrid = obj.hypothesizedInterface;
            if ~strcmpi(obj.network.surveillanceMode, "rotating")
                backOfArray = all(interfaceGrid.transmitBackOfArrayTarget, 3) & all(interfaceGrid.receiveBackOfArrayTarget, 3);
            else
                backOfArray = false(prod(obj.gridSize), 1);
            end
            currentSynchronization = obj.network.fractionalDelayMode;
            cleanup = onCleanup(@() resetsynchronization(obj, currentSynchronization));
            obj.network.settingsnetwork('fractionalDelayMode', 'off');
            uncoveredZone = squeeze(~any(interfaceGrid.waveformReceivedFromScatterers, [1 2 4]));
            obj.blindZone = backOfArray | uncoveredZone;
            function resetsynchronization(obj, currentSynchronization)
                obj.network.settingsnetwork('fractionalDelayMode', currentSynchronization);
            end
        end

        function setintegrationindices(obj)
            if isempty(obj.blindZone)
                obj.setblindzone;
            end
            visibleZone = ~obj.blindZone;
            numberOfVisibleCells = nnz(visibleZone);
            hypothesizedTimeDelays = obj.hypothesizedInterface.timeDelay(:, :, visibleZone); % (Ntx x Nrx x Ni matrix)
            obj.expectedCellObservations = 10.^(.5*obj.hypothesizedInterface.receivedPowerFromScatterers_dBW(:, :, visibleZone)); % (Ntx x Nrx x Ni matrix)
            % obj.expectedCellObservations = obj.expectedCellObservations.*exp(1j*phase); % not implemented
            obj.integrationIndices = ones(obj.network.numberOfActiveTransmittingNodes, obj.network.numberOfActiveReceivingNodes, numberOfVisibleCells);
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
                    timeDifference = abs(hypothesizedTimeDelays(txID, rxID, :) - timeDelays); % (Ntx x Nrx x Ni matrix)
                    [~, obj.integrationIndices(txID, rxID, :)] = min(timeDifference, [], 4);
                end
            end
        end

        function setdictionary(obj)
            % (Ns_i x Ni x Nrx x Ntx matrix)
            if isempty(obj.blindZone)
                obj.setblindzone;
            end
            visibleZone = ~obj.blindZone;
            numberOfVisibleCells = nnz(visibleZone);
            currentSynchronization = obj.network.fractionalDelayMode;
            cleanup = onCleanup(@() resetsynchronization(obj, currentSynchronization));
            obj.network.settingsnetwork('fractionalDelayMode', 'off');
            signals = obj.hypothesizedInterface.signalReceivedFromScatterers; % Ns x 1 x Nt x Ntx x M
            % pos = hypothesizedPositions(:, (all(abs(squeeze(signals{1}(:, :, :, :, 1))) == 0)));
            % figure; plot(pos(1, :), pos(2, :), '.');
            Ns = size([obj.network.activeReceivingNodes.samplingInstants], 1);
            switch obj.network.networkMode
                case 'multiStatic'
                    dict = zeros(Ns, obj.network.numberOfActiveReceivingNodes, numberOfVisibleCells, obj.network.numberOfActiveTransmittingNodes);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        dict(:, rxID, :, :) = signals{rxID}(:, :, visibleZone, :, :);
                    end
                    obj.dictionaryCompression = permute(dict, [1 3 2 4]); % (Ns_i x Ni x Nrx x Ntx matrix)
                case 'monoStatic'
                    dict = cell(1, obj.network.numberOfActiveReceivingNodes);
                    for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                        dict{rxID} = zeros(Ns, obj.network.activeReceivingNodes(rxID).numberOfTotalChannels, numberOfVisibleCells);
                        for elementID = 1 : obj.network.activeReceivingNodes(rxID).numberOfTotalChannels
                            dict{rxID}(:, elementID, :) = signals{rxID}(:, :, visibleZone, obj.network.monoStaticTransmitterIDs(rxID), elementID);
                        end
                        dict{rxID} = permute(dict{rxID}, [1 3 2 4]); % (Ns_i x Ni x Nrxch matrix)
                    end
                    obj.dictionaryCompression = dict;
            end
            function resetsynchronization(obj, currentSynchronization)
                obj.network.settingsnetwork('fractionalDelayMode', currentSynchronization);
            end
        end

        function setmatchfilteredsignals(obj)
            % (1 x Nrx cell of Ns + L - 1 x Nrxch x Ntxch Ncmp matrix) Ncmp : number of parallel trials
            obj.signalsMatchFiltered = obj.signalsMatchFilteredTrials;
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
                options.numberOfTargets (1, 1) double = obj.configurationCompression.numberOfTargets
                options.neighbourOffset (1, 1) double {mustBeNonnegative} = obj.configurationCompression.neighbourOffset
                options.processingAlgorithm (1, 1) {mustBeInteger, mustBeInRange(options.processingAlgorithm, 1, 6)} = obj.processingAlgorithm
                options.detectionAlgorithm (1, 1) string {mustBeMember(options.detectionAlgorithm, ["thresholding", "peak", "CoSaMP", "OMP", "OMPjoint"])} = obj.detectionAlgorithm
                options.numberOfTrials (1, 1) {mustBeNonnegative, mustBeInteger} = obj.monteCarlo.numberOfTrials
                options.numberOfTrialsParallel (1, 1) {mustBeNonnegative, mustBeInteger} = obj.monteCarlo.numberOfTrialsParallel
                options.seed (1, 1) {mustBeNonnegative, mustBeInteger, mustBeLessThan(options.seed, 4294967296)} = 0
                options.seedShuffle (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.seedShuffle, [0, 1])} = obj.seedShuffle
            end
            obj.configuration.PFA = options.PFA;
            obj.configurationCompression.numberOfTargets = options.numberOfTargets;
            obj.configurationCompression.neighbourOffset = options.neighbourOffset;
            obj.processingAlgorithm = options.processingAlgorithm;
            switch obj.network.networkCoherency
                case "short-term coherent"
                    if any(obj.processingAlgorithm, [1 2])
                        warning('Processing algorithm %d cannot be used for short-term coherent networks', obj.processingAlgorithm);
                    end
                case "incoherent"
                    if any(obj.processingAlgorithm, [1 2 4])
                        warning('Processing algorithm %d cannot be used for incoherent networks', obj.processingAlgorithm);
                    end
            end
            obj.detectionAlgorithm = options.detectionAlgorithm;
            mc = obj.monteCarlo;
            mc.numberOfTrials = options.numberOfTrials;
            mc.numberOfTrialsParallel = options.numberOfTrialsParallel;
            obj.interfaces.numberOfTrialsParallel = options.numberOfTrialsParallel;
            mc.seed = options.seed;
            obj.monteCarlo = mc;
            obj.seedShuffle = options.seedShuffle;
            rng(obj.monteCarlo.seed);
        end

        function configuremonostatic(obj, options)
            arguments
                obj
                options.removeBoundaryDetectionDF (1, 1) double {mustBeNonnegative} = obj.configuration.removeBoundaryDetectionDF
                options.removeBoundaryDetectionMF (1, 1) double {mustBeNonnegative} = obj.configuration.removeBoundaryDetectionMF
                options.interpolationDF (1, 1) double {mustBeNonnegative} = obj.configuration.interpolationDF
                options.interpolationMF (1, 1) double {mustBeNonnegative} = obj.configuration.interpolationMF
                options.groupingCloseness (1, 1) double {mustBePositive} = obj.configuration.groupingCloseness
                options.numberOfTrainingCells (1, 1) {mustBePositive, mustBeInteger} = obj.configuration.numberOfTrainingCells
                options.numberOfGuardCells (1, 1) {mustBeNonnegative, mustBeInteger} = obj.configuration.numberOfGuardCells
                options.CFAR (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.CFAR, [0, 1])} = obj.configuration.CFAR
            end
            obj.configurationMonostatic.removeBoundaryDetectionDF = options.removeBoundaryDetectionDF;
            obj.configurationMonostatic.removeBoundaryDetectionMF = options.removeBoundaryDetectionMF;
            obj.configurationMonostatic.interpolationDF = options.interpolationDF;
            obj.configurationMonostatic.interpolationMF = options.interpolationMF;
            obj.configurationMonostatic.groupingCloseness = options.groupingCloseness;
            obj.configurationMonostatic.numberOfTrainingCells = options.numberOfTrainingCells;
            obj.configurationMonostatic.numberOfGuardCells = options.numberOfGuardCells;
            obj.configurationMonostatic.CFAR = options.CFAR;
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
                    case {"OMPjoint", "OMP"}
                        obj.applycompression;
                    case "CoSaMP"
                        error("not implemented");
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
                options.meanRCS_dbsm (1, 1) double = 0
                options.onCellCenters (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.onCellCenters, [0, 1])} = 0; 
            end
            mc = obj.monteCarlo;
            numberOfTotalTrials = mc.numberOfTrials*mc.numberOfTrialsParallel;
            originalTargets = obj.interfaces.targets;
            originalNumberOfTargets = obj.configurationCompression.numberOfTargets;
            previousCompressionReport = obj.compressionReport;
            cleanup = onCleanup(@() cleanupFunction(obj, originalTargets, originalNumberOfTargets, previousCompressionReport));
            obj.configure("numberOfTargets", nan);
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
                'PD', [], ...
                'PFA', [], ...
                'SNRs', [], ...
                'SNRsMean', []);
            obj.coverageSimulationReport.PD = zeros(prod(obj.gridSize), 1);
            obj.coverageSimulationReport.SNRs = -inf(prod(obj.gridSize), 1);
            obj.coverageSimulationReport.SNRsMean = zeros(prod(obj.gridSize), 1);
            visibleCellIDs = 1 : size(cellPositions, 2); %% ayar
            numberOfCells = length(visibleCellIDs);
            targetCellIDs = find(visibleZone);
            targetCellIDs = targetCellIDs(visibleCellIDs);
            obj.coverageSimulationReport.targetCellIDs = targetCellIDs;
            targetPositions = zeros(3, numberOfCells);
            targetPositions(dimensions, :) = cellPositions(dimensions, visibleCellIDs);
            for targetID = 1 : numberOfCells
                targetCellID = targetCellIDs(targetID);
                targetNeighbourCellIDs = obj.neighbourssquarewindow(targetCellID, "offset", obj.configurationCompression.neighbourOffset);
                obj.interfaces.settargets(target( ...
                    "position", targetPositions(:, targetID), ...
                    "meanRCS_dbsm", options.meanRCS_dbsm));
                obj.coverageSimulationReport.SNRs(targetCellID) = max(10*log10(mean(obj.outputSNR_lin, 4)), [], 2);
                for mcID = 1 : mc.numberOfTrials
                    if ~options.onCellCenters
                        obj.interfaces.settargetpositions("width", width);
                    end
                    switch obj.detectionAlgorithm
                        case {"thresholding", "peak"}
                            obj.setmatchfilteredsignals;
                        case {"OMPjoint", "OMP"}
                            obj.applycompression;
                        case "CoSaMP"
                            error("not implemented");
                    end
                    idx = obj.hypothesisTestingResults;
                    switch obj.network.networkMode
                        case 'multiStatic'
                            for mcpID = 1 : mc.numberOfTrialsParallel
                                targetIndex = ismember(idx{mcpID}, targetNeighbourCellIDs);
                                if any(targetIndex)
                                    obj.coverageSimulationReport.PD(targetCellID) = obj.coverageSimulationReport.PD(targetCellID) + 1/numberOfTotalTrials;
                                    obj.coverageSimulationReport.SNRsMean(targetCellID) = obj.coverageSimulationReport.SNRsMean(targetCellID) + mean(obj.compressionReport(1, mcpID).maximumPowerHistory(find(targetIndex)));
                                end
                            end
                        case 'monoStatic'
                            for mcpID = 1 : mc.numberOfTrialsParallel
                                SNRmax = -inf;
                                for fusionID = 1 : obj.network.numberOfParallelSignalFusions
                                    targetIndex = ismember(idx{fusionID}{mcpID}, targetNeighbourCellIDs);
                                    if any(targetIndex)
                                        SNRmonostatic = mean(obj.compressionReport(fusionID, mcpID).maximumPowerHistory(find(targetIndex)));
                                        if SNRmonostatic > SNRmax
                                            SNRmax = SNRmonostatic;
                                        end
                                    end
                                end
                                if ~isinf(SNRmax)
                                    obj.coverageSimulationReport.PD(targetCellID) = obj.coverageSimulationReport.PD(targetCellID) + 1/numberOfTotalTrials;
                                    obj.coverageSimulationReport.SNRsMean(targetCellID) = obj.coverageSimulationReport.SNRsMean(targetCellID) + SNRmax;
                                end
                            end
                    end
                    if ~mod(mcID, 100)
                        fprintf('trial = %d/%d\n', mcID, mc.numberOfTrials);
                    end
                end
                obj.coverageSimulationReport.SNRsMean(targetCellID) = 10*log10(obj.coverageSimulationReport.SNRsMean(targetCellID)./obj.coverageSimulationReport.PD(targetCellID)/numberOfTotalTrials);
                if ~mod(targetID, 10)
                    fprintf('target = %d/%d\n', targetID, numberOfCells);
                end
            end
            obj.coverageSimulationReport.SNRsMean(~logical(obj.coverageSimulationReport.PD)) = -inf;
            function cleanupFunction(obj, originalTargets, originalNumberOfTargets, previousCompressionReport)
                obj.compressionReport = previousCompressionReport;
                obj.interfaces.targets = originalTargets;
                obj.interfaces.settargetpositions("width", 0);
                obj.configure("numberOfTargets", originalNumberOfTargets);
            end
        end

        function visualizecoveragesimulation(obj, options)
            arguments
                obj
                options.saveFigures (1, 1) {mustBeNumericOrLogical, mustBeMember(options.saveFigures, [0, 1])} = false
                options.monoStaticNetworkRXID (1, 1) {mustBePositive, mustBeInteger} = 1
                options.saveDirectory (1, :) {mustBeText} = ''
            end
            if isempty(obj.coverageSimulationReport.targetCellIDs)
                fprintf('coverage simulation had not executed\n');
                return;
            end
            if ~isfolder(options.saveDirectory)
                mkdir(options.saveDirectory);
            end
            config = obj.configuration;
            gridScan = obj.gridPointsMesh;
            dimensions = {"x", "y", "z"};
            dimensions = dimensions(size(gridScan.x) ~= 1);
            xLabel = dimensions{1} + " (km)";
            yLabel = dimensions{2} + " (km)";
            x1 = obj.gridPoints{1}/1e3;
            x2 = obj.gridPoints{2}/1e3;
            visibleZone = reshape(~obj.blindZone, obj.gridSize([2 1 3]));

            PD = reshape(obj.coverageSimulationReport.PD(:, options.monoStaticNetworkRXID), obj.gridSize([2 1 3]));
            SNRs = reshape(obj.coverageSimulationReport.SNRs(:, options.monoStaticNetworkRXID), obj.gridSize([2 1 3]));
            SNRsMean = reshape(obj.coverageSimulationReport.SNRsMean(:, options.monoStaticNetworkRXID), obj.gridSize([2 1 3]));

            fig = figure;
            img = imagesc(x1, x2, PD);
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
            img.DataTipTemplate.DataTipRows(3).Value = PD;
            savePath = fullfile(options.saveDirectory, 'PD.fig');
            hold off; drawnow;
            if ~exist(savePath, 'file')
                savefig(fig, savePath);
            end

            fig = figure;
            img = imagesc(x1, x2, SNRs);
            colorbar; colormap('default');
            ax = gca; set(ax, 'Ydir', 'Normal');
            set(img, 'AlphaData', visibleZone);
            delete(datatip(img, 2, 2));
            grid off; grid on; grid minor;
            xlabel(xLabel); ylabel(yLabel);
            title('modeled SNR averaged over trials'); hold off;
            img.DataTipTemplate.DataTipRows(1).Label = "x";
            img.DataTipTemplate.DataTipRows(1).Value = gridScan.x;
            img.DataTipTemplate.DataTipRows(2).Label = "y";
            img.DataTipTemplate.DataTipRows(2).Value = gridScan.y;
            img.DataTipTemplate.DataTipRows(3).Label = "mean SNR";
            img.DataTipTemplate.DataTipRows(3).Value = SNRs;
            hold off; drawnow;
            savePath = fullfile(options.saveDirectory, 'SNRmodeled.fig');
            if ~exist(savePath, 'file')
                savefig(fig, savePath);
            end

            fig = figure;
            img = imagesc(x1, x2, SNRsMean);
            colorbar; colormap('default');
            ax = gca; set(ax, 'Ydir', 'Normal');
            set(img, 'AlphaData', visibleZone);
            delete(datatip(img, 2, 2));
            grid off; grid on; grid minor;
            xlabel(xLabel); ylabel(yLabel);
            title('realized SNR with straddle loss averaged over trials'); hold off;
            img.DataTipTemplate.DataTipRows(1).Label = "x";
            img.DataTipTemplate.DataTipRows(1).Value = gridScan.x;
            img.DataTipTemplate.DataTipRows(2).Label = "y";
            img.DataTipTemplate.DataTipRows(2).Value = gridScan.y;
            img.DataTipTemplate.DataTipRows(3).Label = "mean SNR";
            img.DataTipTemplate.DataTipRows(3).Value = SNRsMean;
            hold off; drawnow;
            savePath = fullfile(options.saveDirectory, 'SNRrealized.fig');
            if ~exist(savePath, 'file')
                savefig(fig, savePath);
            end

            fig = figure; hold on;
            plot(SNRs(:), SNRsMean(:), '.b');
            minSNR = min(SNRs(:), SNRsMean(:));
            maxSNR = max(SNRs(:), SNRsMean(:));
            plot([minSNR, maxSNR], [minSNR, maxSNR]);
            grid off; grid on; grid minor;
            xlabel('modeled SNR_{out} (dB)'); ylabel('realized SNR_{out} (dB)');
            title('modeled SNR vs realized SNR with straddle loss both averaged over trials');
            hold off; drawnow;
            savePath = fullfile(options.saveDirectory, 'SNRcomparison.fig');
            if ~exist(savePath, 'file')
                savefig(fig, savePath);
            end

            fig = figure; hold on;
            plot(SNRs(:), PD(:), '.b');
            plot(SNRs(:), obj.ROC(config.PFA, 10.^(.1*SNRs(:)), obj.network.numberOfActiveBistaticPairs), '.r');
            grid off; grid on; grid minor;
            xlabel('modeled SNR_{out} (dB)'); ylabel('p_D');
            title('p_D vs modeled SNR averaged over trials');
            legend('simulation', 'theoretical', 'Location', 'best');
            hold off; drawnow;
            savePath = fullfile(options.saveDirectory, 'ROCmodeled.fig');
            if ~exist(savePath, 'file')
                savefig(fig, savePath);
            end

            fig = figure; hold on;
            plot(SNRsMean(:), PD(:), '.b');
            plot(SNRsMean(:), obj.ROC(config.PFA, 10.^(.1*SNRsMean(:)), obj.network.numberOfActiveBistaticPairs), '.r');
            grid off; grid on; grid minor;
            xlabel('realized SNR_{out} (dB)'); ylabel('p_D');
            title('p_D vs realized SNR averaged over trials');
            legend('simulation', 'theoretical', 'Location', 'best');
            hold off; drawnow;
            savePath = fullfile(options.saveDirectory, 'ROCrealized.fig');
            if ~exist(savePath, 'file')
                savefig(fig, savePath);
            end

            figure;
            PD(~visibleZone) = nan;
            contour(x1, x2, PD, 0 : .1 : 1, 'LineWidth', 2, 'ShowText', 'on');
            grid on; grid minor;
            xlabel(xLabel); ylabel(yLabel); zlabel('p_D');
            title('Probability of detection');
            hold off; drawnow;
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
            cleanup = onCleanup(@() recoverstate(obj, currentNumberOfTrials, currentNumberOfTrialsParallel));
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
                            numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
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
            function recoverstate(obj, currentNumberOfTrials, currentNumberOfTrialsParallel)
                obj.configure( ...
                    "numberOfTrials", currentNumberOfTrials, ...
                    "numberOfTrialsParallel", currentNumberOfTrialsParallel);
                obj.resetsimulation;
            end
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
            if isscalar(obj.signalsMatchFiltered{1})
                obj.setmatchfilteredsignals;
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
                            L = obj.network.pulseWidthSample(rxID, txID) - 2;
                            N = obj.network.activeReceivingNodes(rxID).numberOfSamplesPerCPI;
                    end
                    t = (-L : N)*Ts*1e6;
                    if ~isscalar(options.receivingNodeIDs)
                        plot(t, 20*log10(abs(s(:, ceil(end/2), txID))));
                    else
                        plot(t, 20*log10(abs(s(:, :, txID))));
                    end
                    hold on;
                    if strcmpi(obj.network.networkMode, "monoStatic") && obj.processingAlgorithm == 2 && obj.configurationMonostatic.CFAR
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
                if obj.processingAlgorithm ~= 2 || ~obj.configurationMonostatic.CFAR
                    yline(config.threshold_dB, 'LineStyle', '--', 'LineWidth', 2);
                end
            end
            grid off; grid on; grid minor;
            xlabel('time (\mus)'); ylabel('power (dB)');
            title('filtered signal');
            if ~isscalar(options.receivingNodeIDs) && ~(strcmpi(obj.network.networkMode, "monoStatic") && obj.processingAlgorithm == 2 && obj.configurationMonostatic.CFAR)
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
                options.variable (1, 1) string {mustBeMember(options.variable, ["time", "amplitude", "phase", "coveredZone", "indices", "weightsAmplitude", "weightsPhase"])} = "time"
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
            switch options.variable
                case "time"
                    var = squeeze(obj.hypothesizedInterface.timeDelay(options.transmittingNodeID, options.receivingNodeID, :))*1e6;
                    labelVariable = 'time delay (\mus)';
                    titleVariable = 'hypothesized time delays';
                case "amplitude"
                    var = obj.hypothesizedInterface.receivedPowerFromScatterers_dBW(options.transmittingNodeID, options.receivingNodeID, :);
                    labelVariable = 'amplitude (dB)';
                    titleVariable = 'hypothesized amplitudes';
                case "phase"
                    %%% not implemented
                    % var = 180*squeeze(obj.hypothesizedPhases(options.transmittingNodeID, options.receivingNodeID, :))/pi;
                    labelVariable = 'phase (°)';
                    titleVariable = 'hypothesized phases';
                case "coveredZone"
                    var = ~obj.blindZone;
                    labelVariable = 'isVisible';
                    titleVariable = 'visibility of cells';
                case "indices"
                    var = zeros(prod(obj.gridSize), 1);
                    var(~obj.blindZone) = squeeze(obj.integrationIndices(options.transmittingNodeID, options.receivingNodeID, :));
                    labelVariable = 'index';
                    titleVariable = 'integration indices';
                case "weightsAmplitude"
                    var = zeros(prod(obj.gridSize), 1);
                    var(~obj.blindZone) = 20*log10(abs(squeeze(obj.integrationWeights(options.transmittingNodeID, options.receivingNodeID, :))));
                    labelVariable = 'amplitude (dB)';
                    titleVariable = 'amplitude of integration weights)';
                case "weightsPhase"
                    var = zeros(prod(obj.gridSize), 1);
                    var(~obj.blindZone) = 180*angle(squeeze(obj.integrationWeights(options.transmittingNodeID, options.receivingNodeID, :)))/pi;
                    labelVariable = 'phase (°)';
                    titleVariable = 'phase of integration weights)';
            end
            switch options.dimension
                case "x-y"
                    var = reshape(var, [obj.gridSize]);
                    x = gridScan.x(:, :, 1); x = x(:)/1e3;
                    y = gridScan.y(:, :, 1); y = y(:)/1e3;
                    varPlot = var(:, :, 1);
                    plot3(x, y, varPlot(:), '.'); hold on;
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

        function visualizedictionary(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.cellInterested double {mustBeInteger, mustBePositive} = []
                options.transmittingNodeID (1, 1) double {mustBeInteger, mustBePositive} = 1
                options.plot (1, 1) string {mustBeMember(options.plot, ["magnitude", "phase"])} = "magnitude"
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            if isempty(options.cellInterested)
                switch obj.network.networkMode
                    case 'multiStatic'
                        dict = obj.dictionaryCompression(:, :, 1, options.transmittingNodeID);
                    case 'monoStatic'
                        dict = obj.dictionaryCompression{options.transmittingNodeID}(:, :, 1);
                end
            elseif obj.network.numberOfActiveReceivingNodes > 1
                switch obj.network.networkMode
                    case 'multiStatic'
                        dict = squeeze(obj.dictionaryCompression(:, options.cellInterested, :, options.transmittingNodeID));
                    case 'monoStatic'
                        dict = squeeze(obj.dictionaryCompression{options.transmittingNodeID}(:, options.cellInterested, :));
                end
            end
            if isempty(options.cellInterested) || obj.network.numberOfActiveReceivingNodes > 1
                if isempty(options.cellInterested)
                    xlabel('cells'); ylabel('sample');
                elseif obj.network.numberOfActiveReceivingNodes > 1
                    xlabel('rxID'); ylabel('sample');
                end
                switch options.plot
                    case "magnitude"
                        imagesc(abs(dict)); colorbar;
                    case "phase"
                        imagesc(180*angle(dict)/pi); colorbar;
                end
            else
                switch options.plot
                    case "magnitude"
                        plot(abs(dict), 'b');
                        xlabel('sample'); ylabel('magnitude');
                        title('magnitude of dictionary');
                    case "phase"
                        plot(180*angle(dict)/pi, '*r');
                        xlabel('sample'); ylabel('phase (degree)');
                        title('phase of dictionary');
                end
            end
            switch options.plot
                case "magnitude"
                    title('magnitude of dictionary');
                case "phase"
                    title('phase of dictionary');
            end
        end

        function visualizeintegratedsignals(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.scale (1, 1) string {mustBeMember(options.scale, ["linear", "power", "dB"])} = "dB"
                options.plotMode (1, 1) string {mustBeMember(options.plotMode, ["image", "mesh"])} = "image"
                options.trialID (1, 1) double {mustBePositive, mustBeInteger} = 1
                options.monoStaticNetworkRXID (1, 1) {mustBePositive, mustBeInteger} = 1
                options.monoStaticNetworkCHID (1, 1) {mustBePositive, mustBeInteger} = 1
                options.iterationID (1, :) {mustBeNonnegative, mustBeInteger} = 0
            end
            switch obj.network.networkMode
                case "multiStatic"
                    options.monoStaticNetworkRXID = 1;
            end
            currentDetectionAlgorithm = obj.detectionAlgorithm;
            cleanup = onCleanup(@() resetdetectionalgorithm(obj, currentDetectionAlgorithm));
            config = obj.configuration;
            if ~isempty(obj.compressionReport) && ~isempty(obj.compressionReport(options.monoStaticNetworkRXID, 1).integratedSignals)
                options.iterationID = unique(min(options.iterationID, obj.compressionReport(options.monoStaticNetworkRXID, options.trialID).numberOfIterations));
            else
                options.iterationID = 0;
            end
            if any(options.iterationID > 0) && ~isempty(obj.compressionReport(options.monoStaticNetworkRXID, options.trialID).integratedSignals)
                processedSignalType = 'compressed';
                if any(strcmpi(obj.detectionAlgorithm, ["peak", "thresholding"]))
                    switch obj.network.networkMode
                        case "multiStatic"
                            obj.detectionAlgorithm = "OMP";
                        case "monoStatic"
                            obj.detectionAlgorithm = "OMPjoint";
                    end
                end
                z = nan(prod(obj.gridSize), 1);
                if isscalar(options.iterationID)
                    z(~obj.blindZone) = obj.compressionReport(options.monoStaticNetworkRXID, options.trialID).integratedSignals(:, options.iterationID);
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
                processedSignalType = 'integrated';
                if any(strcmpi(obj.detectionAlgorithm, ["OMP", "OMPjoint", "CoSaMP"]))
                    obj.detectionAlgorithm = "peak";
                end
                if isempty(obj.signalsIntegrated)
                    fprintf('"signalsIntegrated" is empty\n');
                    return;
                end
                switch obj.network.networkMode
                    case "multiStatic"
                        z = obj.signalsIntegrated(:, :, options.trialID);
                    case "monoStatic"
                        z = obj.signalsIntegrated{options.monoStaticNetworkRXID}(:, options.monoStaticNetworkCHID, options.trialID);
                end
            end
            switch obj.network.networkMode
                case "multiStatic"
                    if options.iterationID 
                        pos = obj.estimatedPositions{options.trialID};
                        if ~isempty(pos)
                            pos = pos(:, unique(min(size(pos, 2), options.iterationID)));
                        else
                            pos = double.empty(3, 0);
                        end
                    else
                        pos = obj.estimatedPositions{options.trialID};
                    end
                case "monoStatic"
                    if options.iterationID
                        pos = obj.estimatedPositions{options.monoStaticNetworkRXID}{options.monoStaticNetworkCHID, options.trialID};
                        if size(pos, 2) >= options.iterationID
                            pos = pos(:,  options.iterationID);
                        else
                            pos = double.empty(3, 0);
                        end
                    else
                        pos = obj.estimatedPositions{options.monoStaticNetworkRXID}{options.monoStaticNetworkCHID, options.trialID};
                    end
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
                    plot(x1, 20*log10(abs(z)));
                    grid off; grid on; grid minor;
                    xlabel(xLabel); ylabel('power (dB)');
                    title([processedSignalType ' signal']);
                case 2
                    zLabel = 'z_{processed}';
                    Y = reshape(abs(z), obj.gridSize([2, 1, 3]));
                    switch options.scale
                        case "linear"
                            zLabel = [zLabel, ', linear'];
                            Y = sqrt(Y);
                        case "power"
                            zLabel = [zLabel, ', power'];
                        case "dB"
                            zLabel = [zLabel, ', power (dB)'];
                            Y = 10*log10(Y);
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
                    xlim(obj.network.boundaryListened(1, :)/1e3); ylim(obj.network.boundaryListened(2, :)/1e3);
                    grid on; grid minor; view(0, 90);
                    xlabel(xLabel); ylabel(yLabel); zlabel(zLabel);
                    switch obj.network.networkMode
                        case "multiStatic"
                            SNR = max(obj.outputSNR_dB, [], [3 4]);
                            PD = max(config.PD, [], [3 4]);
                        case "monoStatic"
                            SNR = max(obj.outputSNR_dB(:, options.monoStaticNetworkRXID, :, :), [], [3 4]);
                            PD = max(config.PD(:, options.monoStaticNetworkRXID, :, :), [], [3 4]);
                    end
                    title([processedSignalType sprintf(' signal SNR = %.3g dB, PFA = %.3e, PD = %.3g', SNR, config.PFA, PD)]); hold off;
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
            function resetdetectionalgorithm(obj, previousAlgorithm)
                obj.detectionAlgorithm = previousAlgorithm;
            end
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
            if isempty(obj.estimatedPositions)
                fprintf('"estimatedPositions" is empty\n');
                return;
            end
            switch obj.network.networkMode
                case "multiStatic"
                    options.monoStaticNetworkRXID = 1;
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
            z = abs(obj.compressionReport(options.monoStaticNetworkRXID, options.trialID).complexAmplitudes);
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
                    xlim(obj.network.boundaryListened(1, :)/1e3); ylim(obj.network.boundaryListened(2, :)/1e3);
                    grid off; grid on; grid minor;
                    xlabel(xLabel); ylabel(yLabel);
                    switch obj.network.networkMode
                        case "multiStatic"
                            SNR = max(obj.outputSNR_dB, [], [3 4]);
                            PD = max(config.PD, [], [3 4]);
                        case "monoStatic"
                            SNR = max(obj.outputSNR_dB(:, options.monoStaticNetworkRXID, :, :), [], [3 4]);
                            PD = max(config.PD(:, options.monoStaticNetworkRXID, :, :), [], [3 4]);
                    end
                    title(sprintf('compressed signal SNR = %.3g dB, PFA = %.3e, PD = %.3g', SNR, config.PFA, PD)); hold off;
                case 3
            end
            if isempty(pos)
                legend('RX', 'TX', 'targets', 'Location', 'best');
            else
                legend('RX', 'TX', 'estimations', 'targets', 'Location', 'best');
            end
            hold off; drawnow;
            figure; plot(10*log10(obj.compressionReport(options.monoStaticNetworkRXID, options.trialID).residualPowerHistory), 'b');
            hold on; yline(config.threshold_dB, 'r--');
            grid off; grid on; grid minor;
            xlabel('iteration ID'); ylabel('residual power (dB)');
            title('residual power vs iteration ID of the compressing sensing algorithm');
            legend('residual power', 'threshold', 'Location', 'best');
            hold off; drawnow;
        end

        function visualizedirectlyintegratedsignals(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.cuttingIndices double = [nan, nan, ceil(([obj.network.activeReceivingNodes(3 : end).numberOfSamplesPerCPI] + obj.network.pulseWidthSample(3 : end, 1).' - 1)/2)]
            end
            cuttingRxIDs = find(~isnan(options.cuttingIndices));
            visualizedRxIDs = setdiff(1 : obj.network.numberOfActiveReceivingNodes, cuttingRxIDs);
            Y = obj.signalsIntegratedDirectly;
            if isempty(Y)
                fprintf('"signalsIntegratedDirectly" is empty\n');
                return;
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            for rxID = cuttingRxIDs
                Y = indexer(Y, rxID, options.cuttingIndices(rxID));
            end
            function out = indexer
                dimY = numel(size(Y));
                expr = repmat(':,', 1, dimY);
                expr = expr(1 : end - 1);
                expr = [expr(1 : 2*rxID - 2), 'idx', expr(2*rxID : end)];
                out = eval(['Y(', expr, ')']);
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
                cleanup = onCleanup(@() cleanupFunction(obj, currentAlgorithm));
                for k = 1 : 6
                    obj.configure("processingAlgorithm", k);
                    PD(:, k) = obj.ROC(PFA, 10.^(.1*options.snr_dB(1)), obj.network.numberOfActiveBistaticPairs);
                end
                figure; semilogx(PFA, PD);
                grid on; grid minor;
                xlabel('PFA'); ylabel('PD');
                leg = legend(num2str((1 : 6).'), 'Location', 'best');
                title(leg, 'algorithm');
                title("Receiver operating characteristic (ROC) curve for swerling " + num2str(obj.interfaces.swerling) + " \color{blue}SNR_{out} = " + num2str(options.snr_dB(1)) + " dB");
            else
                PFA = logspace(-10, 0, nop);
                if ~isscalar(options.snr_dB)
                    PD = zeros(nop, length(options.snr_dB));
                    for k = 1 : length(options.snr_dB)
                        PD(:, k) = obj.ROC(PFA, 10.^(.1*options.snr_dB(k)), obj.network.numberOfActiveBistaticPairs);
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
                    title("Receiver operating characteristic (ROC) curve of algorithm \color{blue}L_{" + num2str(obj.processingAlgorithm) + "} for swerling " + num2str(obj.interfaces.swerling));
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
                    title("Receiver operating characteristic (ROC) curve of algorithm \color{blue}L_{" + num2str(obj.processingAlgorithm) + "} for swerling " + num2str(obj.interfaces.swerling));
                end
            end
            function cleanupFunction(obj, detectionAlgorithm)
                obj.configure("processingAlgorithm", detectionAlgorithm);
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
            title("Detection characteristic curve of algorithm \color{blue}L_{" + num2str(obj.processingAlgorithm) + "} for swerling " + num2str(obj.interfaces.swerling));
        end
    end
end