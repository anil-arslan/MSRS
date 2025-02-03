classdef interface < handle
    %interface Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = ?spu, GetAccess = public)
        network (1, :) radarNetwork = radarNetwork.empty()
        targets (1, :) target = target.empty()
        environment (1, 1) %%%%% not implemented
        configuration (1, 1) struct = struct( ...
            'noise', false, ...
            'directPath', false, ...
            'pathLoss', true)
        spatialCoherency (1, 1) string {mustBeMember(spatialCoherency, ["deterministic", "coherent", "correlatedAmplitudeFluctuations", "incoherent"])} = "deterministic"
        positions double % meters (3 x Nt x Nmcp matrix)
    end

    properties (Dependent)
        %%% target parameters
        bistaticRange double % meters (Ntx x Nrx x Nt x Nmcp matrix) major axis length
        timeDelay double % seconds (Ntx x Nrx x Nt x Nmcp matrix)
        dopplerShift double % Hz (Ntx x Nrx x Nt x Nmcp matrix)
        dopplerShiftBisector double % Hz (Ntx x Nrx x Nt x Nmcp matrix)

        transmitBackOfArrayTarget double % (Nt x Nmcp x Ntx matrix)
        transmitSteeringVector cell % (Ntx x 1 cell of M x Nt x Nmcp matrix)
        transmitBackOfArrayReceiver double % (Nrx x Ntx matrix)
        transmitSteeringVectorDirectPath cell % (Ntx x 1 cell of M x Nrx matrix)

        receiveBackOfArrayTarget double % (Nt x Nmcp x Nrx matrix)
        receiveSteeringVector cell % (Ntx x Nrx cell of M x Nt x Nmcp matrix)
        receiveBackOfArrayTransmitter double % (Ntx x Nrx matrix)
        receiveSteeringVectorDirectPath cell % (Ntx x Nrx cell of M x 1 matrix)

        transmittedBeam dobule % (Ntx x 1 x Nt x Nmcp matrix)
        transmittedBeamDirectPath dobule % (Ntx x 1 x Nrx matrix)
        receivedBeamSpaceObservations cell % (1 x Nrx cell of Ntx x 1 x Nt x Nch x Nmcp matrix)
        transmitGain_lin dobule % linear scale (Ntx x 1 x Nt x Nmcp matrix)
        transmitGainDirectPath_lin dobule % linear scale (Ntx x 1 x Nrx matrix)
        receiveGain_lin cell % linear scale % (1 x Nrx cell of Ntx x 1 x Nt x Nch x Nmcp matrix)

        powerFromScatterersMaximumGain_dBW double % dB Watt (Ntx x Nrx x Nt x Nmcp matrix)
        
        receivedPowerFromScatterers_dBW double % dB Watt (Ntx x Nrx x Nt x Nmcp matrix)
        receivedEnergyFromScatterersPerSample_dBJoule % dB Joule (Ntx x Nrx x Nt x Nmcp matrix)
        receivedEnergyFromScatterersPerPulse_dBJoule % dB Joule (Ntx x Nrx x Nt x Nmcp matrix)
        receivedPowerFromTransmittingNodes_dBW double % dB Watt (Ntx x Nrx matrix)
        inputSNR_dB double % dB (Ntx x Nrx x Nt x Nmcp matrix)
            %%% signal
            waveformReceivedFromScatterers double % complex received waveform (Ns x Nrx x Nt x Ntx x Nmcp matrix)
            waveformReceivedFromTransmitters double % complex received waveform (Ns x Nrx x Ntx matrix)
            signalReceivedFromScatterers % complex received signal (1 x Nrx cell of Ns x 1 x Nt x Ntx x M x Nmcp matrix)
            signalReceivedFromTransmittingNodes double % complex received direct path signal % (1 x Nrx cell of Ns x 1 x Ntx x M matrix)
            signalSuperposed % (1 x Nrx cell of Ns x M x Nmcp matrix)
            signalBeamformed % (1 x Nrx cell of Ns x Nch x Nmcp matrix)

        %%% ellipse
        majorAxis double % meters (Ntx x Nrx x Nt x Nmcp matrix)
        minorAxis double % meters (Ntx x Nrx x Nt x Nmcp matrix)
        eccentricity double % meters (Ntx x Nrx x Nt x Nmcp matrix)
        latusRectum double % meters (Ntx x Nrx x Nt x Nmcp matrix)

        %%% unit vectors
        unitDirectionRX (3, :, :) double % (3 x Nrx x Nt x Nmcp matrix)
        unitDirectionTX (3, :, :) double % (3 x Ntx x Nt x Nmcp matrix)
        unitDirectionBisector (3, :, :, :) double % (3 x Ntx x Nrx x Nt x Nmcp matrix)
        unitDirectionBaseline (3, :, :, :) double % (3 x Ntx x Nrx matrix)
        unitNormalBaseline (3, :, :, :) double % (3 x Ntx x Nrx x Nt x Nmcp matrix)

        %%% distances
        distanceRX (1, :, :, :) double % meters (1 x Nrx x Nt x Nmcp matrix)
        distanceTX (:, 1, :, :) double % meters (Ntx x 1 x Nt x Nmcp matrix)
        distanceBisector double % (Ntx x Nrx x Nt x Nmcp matrix)
        distanceBaseline double % (Ntx x Nrx matrix)

        %%% angles
        bistaticAngle double % (Ntx x Nrx x Nt x Nmcp matrix)
        elevationTargetRX double % (Nrx x Nt x Nmcp matrix)
        azimuthTargetRX double % (Nrx x Nt x Nmcp matrix)
        elevationTargetTX double % (Ntx x Nt x Nmcp matrix)
        azimuthTargetTX double % (Ntx x Nt x Nmcp matrix)
        angleVelocityBisector double % (Ntx x Nrx x Nt x Nmcp matrix)
        angleBaselineRX double % (Ntx x Nrx x Nt x Nmcp matrix)
        angleBaselineTX double % (Ntx x Nrx x Nt x Nmcp matrix)

        %%% factors
        triangulationFactor double {mustBeNonnegative} % (Ntx x Nrx x Nt x Nmcp matrix)
        lossFactor double {mustBeNonnegative} % (Ntx x Nrx x Nt x Nmcp matrix)
        targetRegion string % (Ntx x Nrx x Nt x Nmcp matrix)
        
        %%% sizes
        numberOfTargets (1, 1) double {mustBeInteger, mustBeNonnegative} % Nt
        numberOfReceivingNodes (1, 1) double {mustBeInteger, mustBeNonnegative} % Nrx
        numberOfTransmittingNodes (1, 1) double {mustBeInteger, mustBeNonnegative} % Ntx
        numberOfBistaticPairs (1, 1) double {mustBeInteger, mustBeNonnegative} % Nrx*Ntx
    end

    properties (Constant)
        speedOfLight (1, 1) double {mustBePositive} = physconst('LightSpeed');
    end

    properties (Access = ?spu)
        numberOfTrialsParallel (1, 1) {mustBeNonnegative, mustBeInteger} = 1
    end

    methods
        function obj = interface(options)
            arguments
                options.network (1, :) radarNetwork = radarNetwork.empty()
                options.targets (1, :) target = target.empty()
            end
            %interface Construct an instance of this class
            %   Detailed explanation goes here
            obj.network = options.network;
            obj.targets = options.targets;
            obj.settargetpositions;
        end

        %%% target parameters

        function rho = get.bistaticRange(obj) % major axis
            % (Ntx x Nrx x Nt x Nmcp matrix)
            rho = obj.distanceRX + obj.distanceTX;
        end

        function tau = get.timeDelay(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            tau = obj.bistaticRange/obj.speedOfLight;
        end

        function backOfArray = get.transmitBackOfArrayTarget(obj)
            % (Nt x Nmcp x Ntx matrix)
            transmitArrays = [obj.network.activeTransmittingNodes.array];
            backOfArray = shiftdim(pagemtimes(permute([transmitArrays.normalVector], [1 3 4 2]), 'ctranspose', permute(obj.unitDirectionTX, [1 3 4 2]), 'none') < 0, 1);
            backOfArray(:, :, [transmitArrays.backOfArray]) = false;
        end

        function a = get.transmitSteeringVector(obj)
            % (Ntx x 1 cell of M x Nt x Nmcp matrix)
            transmitArrays = [obj.network.activeTransmittingNodes.array];
            backOfArray = obj.transmitBackOfArrayTarget;
            a = cell(obj.network.numberOfActiveTransmittingNodes, 1);
            for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                % unit vector from radar to target
                a{txID} = exp(1j*2*pi*pagemtimes(transmitArrays(txID).positions, 'ctranspose', permute(obj.unitDirectionTX(:, txID, :, :), [1 3 4 2]), 'none')./obj.network.activeTransmittingNodes(txID).carrierWavelength);
                if ~transmitArrays(txID).backOfArray
                    a{txID}(:, backOfArray(:, :, txID)) = 0;
                end
            end
        end

        function backOfArray = get.transmitBackOfArrayReceiver(obj)
            % (Nrx x Ntx matrix)
            transmitArrays = [obj.network.activeTransmittingNodes.array];
            backOfArray = shiftdim(pagemtimes(permute([transmitArrays.normalVector], [1 3 2]), 'ctranspose', permute(obj.unitDirectionBaseline, [1 3 2]), 'none') < 0, 1);
            backOfArray(:, [transmitArrays.backOfArray]) = false;
        end

        function a = get.transmitSteeringVectorDirectPath(obj)
            % (Ntx x 1 cell of M x Nrx matrix)
            transmitArrays = [obj.network.activeTransmittingNodes.array];
            backOfArray = obj.transmitBackOfArrayReceiver;
            a = cell(obj.network.numberOfActiveTransmittingNodes, 1);
            for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                % unit vector from radar to target
                a{txID} = exp(1j*2*pi*transmitArrays(txID).positions'*permute(obj.unitDirectionBaseline(:, txID, :), [1 3 2])./obj.network.activeTransmittingNodes(txID).carrierWavelength);
                if ~transmitArrays(txID).backOfArray
                    a{txID}(:, backOfArray(:, txID)) = 0;
                end
            end
        end

        function backOfArray = get.receiveBackOfArrayTarget(obj)
            % (Nt x Nmcp x Nrx matrix)
            receiveArrays = [obj.network.activeReceivingNodes.array];
            backOfArray = shiftdim(pagemtimes(permute([receiveArrays.normalVector], [1 3 4 2]), 'ctranspose', permute(obj.unitDirectionRX, [1 3 4 2]), 'none') < 0, 1);
            backOfArray(:, :, [receiveArrays.backOfArray]) = false;
        end

        function a = get.receiveSteeringVector(obj)
            % (Ntx x Nrx cell of M x Nt x Nmcp matrix)
            receiveArrays = [obj.network.activeReceivingNodes.array];
            backOfArray = obj.receiveBackOfArrayTarget;
            a = cell(obj.network.numberOfActiveTransmittingNodes, obj.network.numberOfActiveReceivingNodes);
            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                    a{txID, rxID} = exp(1j*2*pi*pagemtimes(receiveArrays(rxID).positions, 'ctranspose', permute(obj.unitDirectionRX(:, rxID, :, :), [1 3 4 2]), 'none')./obj.network.activeTransmittingNodes(txID).carrierWavelength);
                    if ~receiveArrays(rxID).backOfArray
                        a{txID, rxID}(:, backOfArray(:, :, rxID)) = 0;
                    end
                end
            end
        end

        function backOfArray = get.receiveBackOfArrayTransmitter(obj)
            % (Ntx x Nrx matrix)
            receiveArrays = [obj.network.activeReceivingNodes.array];
            backOfArray = shiftdim(pagemtimes(permute([receiveArrays.normalVector], [1 3 2]), 'ctranspose', obj.unitDirectionBaseline, 'none') < 0, 1);
            backOfArray(:, [receiveArrays.backOfArray]) = false;
        end

        function a = get.receiveSteeringVectorDirectPath(obj)
            % (Ntx x Nrx cell of M x 1 matrix)
            receiveArrays = [obj.network.activeReceivingNodes.array];
            backOfArray = obj.receiveBackOfArrayTransmitter;
            a = cell(obj.network.numberOfActiveTransmittingNodes, obj.network.numberOfActiveReceivingNodes);
            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                    a{txID, rxID} = exp(1j*2*pi*receiveArrays(rxID).positions'*obj.unitDirectionBaseline(:, :, rxID)./obj.network.activeTransmittingNodes(txID).carrierWavelength);
                    if ~receiveArrays(rxID).backOfArray
                        a{txID, rxID}(:, backOfArray(:, rxID)) = 0;
                    end
                end
            end
        end

        function y = get.transmittedBeam(obj)
            % (Ntx x 1 x Nt x Nmcp matrix)
            y = zeros(obj.numberOfTransmittingNodes, 1, obj.numberOfTargets, obj.numberOfTrialsParallel);
            for txID = 1 : obj.numberOfTransmittingNodes
                y(txID, 1, :, :) = permute(pagemtimes(obj.network.activeTransmittingNodes(txID).beamformer, 'ctranspose', obj.transmitSteeringVector{txID}, 'none'), [1 4 2 3]);
            end
        end

        function y = get.transmittedBeamDirectPath(obj)
            % (Ntx x 1 x Nrx matrix)
            y = zeros(obj.numberOfTransmittingNodes, 1, obj.numberOfReceivingNodes);
            for txID = 1 : obj.numberOfTransmittingNodes
                y(txID, 1, :) = permute(obj.network.activeTransmittingNodes(txID).beamformer'*obj.transmitSteeringVectorDirectPath{txID}, [1 3 2]);
            end
        end

        function y = get.receivedBeamSpaceObservations(obj)
            % (1 x Nrx cell of Ntx x 1 x Nt x Nch x Nmcp matrix)
            y = cell(1, obj.network.numberOfActiveReceivingNodes);
            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                y{rxID} = zeros(obj.network.numberOfActiveTransmittingNodes, 1, obj.numberOfTargets, obj.network.activeReceivingNodes(rxID).numberOfTotalChannels, obj.numberOfTrialsParallel);
                for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                    beamformingCenters = permute(exp(1j*2*pi*obj.network.activeReceivingNodes(rxID).beamCenterPositions./obj.network.activeTransmittingNodes(txID).carrierWavelength), [1 3 4 2]);
                    as = exp(1j*2*pi*obj.network.activeReceivingNodes(rxID).array.steeringPositions./obj.network.activeTransmittingNodes(txID).carrierWavelength);
                    y{rxID}(txID, 1, :, :, :) = permute(pagemtimes(obj.network.activeReceivingNodes(rxID).beamformer.*beamformingCenters.*as, 'ctranspose', obj.receiveSteeringVector{txID, rxID}, 'none'), [4 5 2 1 3]);
                end
            end
        end

        function Gt = get.transmitGain_lin(obj)
            % (Ntx x 1 x Nt x Nmcp)
            arrayNodesTX = [obj.network.activeTransmittingNodes.array];
            GtMax = sqrt([obj.network.activeTransmittingNodes.peakTransmitGain].');
            Gt = obj.transmittedBeam.*GtMax./sqrt([arrayNodesTX.numberOfTotalElements]).';
        end

        function Gt = get.transmitGainDirectPath_lin(obj)
            % (Ntx x 1 x Nrx)
            arrayNodesTX = [obj.network.activeTransmittingNodes.array];
            GtMax = sqrt([obj.network.activeTransmittingNodes.peakTransmitGain].');
            Gt = obj.transmittedBeamDirectPath.*GtMax./sqrt([arrayNodesTX.numberOfTotalElements]).';
        end

        function Gr = get.receiveGain_lin(obj)
            % (1 x Nrx cell of Ntx x 1 x Nt x Nch x Nmcp matrix)
            % arrayNodesRX = [obj.network.activeReceivingNodes.array];
            % nodesTX = obj.network.activeTransmittingNodes;
            % GrMax = sqrt(4*pi*[arrayNodesRX.apertureArea]./[nodesTX.carrierWavelength].'.^2);
            % GrMax = [arrayNodesRX.numberOfTotalElements];
            Gr = cell(1, obj.network.numberOfActiveReceivingNodes);
            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                Gr{rxID} = obj.receivedBeamSpaceObservations{rxID};
            end
        end

        function Pr = get.powerFromScatterersMaximumGain_dBW(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            if isempty(obj.targets)
                Pr = -inf(obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
            else
                GtMax = [obj.network.activeTransmittingNodes.peakTransmitGain].';
                powerGains = 10*log10(GtMax) + ...
                    10*log10([obj.network.activeTransmittingNodes.inputPower_W].') + permute(obj.targets.RCS_dbms, [1 3 2]) ...
                    + 20*log10([obj.network.activeTransmittingNodes.carrierWavelength].');
                powerLosses = [obj.network.activeReceivingNodes.systemLoss_dB] + 30*log10(4*pi) + 20*log10(obj.distanceRX.*obj.distanceTX);
                Pr = powerGains - powerLosses;
            end
        end

        function Pr = get.receivedPowerFromScatterers_dBW(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            if isempty(obj.targets)
                Pr = -inf(obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
            else
                powerGains = 20*log10(abs(obj.transmitGain_lin)) + ...
                    10*log10([obj.network.activeTransmittingNodes.inputPower_W].') + permute(obj.targets.RCS_dbms, [1 3 2]) + ...
                    20*log10([obj.network.activeTransmittingNodes.carrierWavelength].');
                powerLosses = [obj.network.activeReceivingNodes.systemLoss_dB] + 30*log10(4*pi) + 20*log10(obj.distanceRX.*obj.distanceTX);
                Pr = powerGains - powerLosses;
            end
        end

        function Er = get.receivedEnergyFromScatterersPerSample_dBJoule(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            Er = obj.receivedPowerFromScatterers_dBW./[obj.network.activeReceivingNodes.samplingFrequency];
        end

        function Er = get.receivedEnergyFromScatterersPerPulse_dBJoule(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            Er = obj.receivedEnergyFromScatterersPerSample_dBJoule.*obj.network.pulseWidthSample.';
        end

        function Pr = get.receivedPowerFromTransmittingNodes_dBW(obj)
            % (Ntx x Nrx matrix)
            powerGains = permute(20*log10(abs(obj.transmitGainDirectPath_lin)), [1 3 2]) + ...
                10*log10([obj.network.activeTransmittingNodes.inputPower_W].') + ... % no RCS term
                20*log10([obj.network.activeTransmittingNodes.carrierWavelength].');
            powerLosses = [obj.network.activeReceivingNodes.systemLoss_dB] + 20*log10(4*pi*obj.distanceBaseline);
            Pr = powerGains - powerLosses;
        end

        function SNRin = get.inputSNR_dB(obj)
            % (Ntx x Nrx x Nt matrix)
            SNRin = obj.powerFromScatterersMaximumGain_dBW - [obj.network.activeReceivingNodes.noisePowerPerSample_dB];
        end

        %%% doppler

        function fd = get.dopplerShift(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            if isempty(obj.targets)
                fd = zeros(obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
            else
                vdRx = pagemtimes(permute(obj.targets.velocity, [3 1 5 2 4]), permute(obj.unitDirectionRX, [1 5 2 3 4]));
                vdTx = pagemtimes(permute(obj.targets.velocity, [3 1 5 2 4]), permute(obj.unitDirectionTX, [1 2 5 3 4]));
                fd = shiftdim((vdRx + vdTx), 1)./[obj.network.activeTransmittingNodes.carrierWavelength].';
            end
        end

        function fd = get.dopplerShiftBisector(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            if isempty(obj.targets)
                fd = zeros(obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
            else
                fd = 2*shiftdim(obj.targets.velocityMagnitude, -1).*cosd(obj.bistaticAngle/2).*cosd(obj.angleVelocityBisector)./[obj.network.activeTransmittingNodes.carrierWavelength].';
            end
        end

        %%% signal

        function waveforms = get.waveformReceivedFromTransmitters(obj)
            % (Ns x Nrx x Ntx matrix)
            tau = [obj.network.activeReceivingNodes.samplingInstants] - permute(obj.network.directPathDelay, [3 2 1]); % Ns x Nrx x Ntx matrix
            fc = shiftdim([obj.network.activeTransmittingNodes.carrierFrequency], -1); % 1 x 1 x Ntx matrix
            waveforms = zeros(size(tau)); % Ns x Nrx x Ntx matrix
            for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                waveforms(:, :, txID) = obj.network.activeTransmittingNodes(txID).waveform(tau(:, :, txID));
            end
            waveforms = waveforms.*exp(-1j*2*pi*fc.*tau); % Ns x Nrx x Ntx matrix
        end

        function waveforms = get.waveformReceivedFromScatterers(obj)
            % (Ns x Nrx x Nt x Ntx x Nmcp matrix)
            tau = [obj.network.activeReceivingNodes.samplingInstants] - permute(obj.timeDelay, [5 2 3 1 4]); % Ns x Nrx x Nt x Ntx x Nmcp matrix
            fd = permute(obj.dopplerShift, [5 2 3 1 4]); % 1 x Nrx x Nt x Ntx matrix
            fc = shiftdim([obj.network.activeTransmittingNodes.carrierFrequency], -2); % 1 x 1 x 1 x Ntx matrix
            waveforms = zeros(size(tau)); % Ns x Nrx x Nt x Ntx matrix
            for txID = 1 : obj.numberOfTransmittingNodes
                waveforms(:, :, :, txID, :) = obj.network.activeTransmittingNodes(txID).waveform(-tau(:, :, :, txID, :));
            end
            waveforms = waveforms.*exp(-1j*2*pi*(fc + fd).*tau);
            waveforms(isinf(tau)) = 0;
        end

        function s = get.signalReceivedFromScatterers(obj)
            % (1 x Nrx cell of Ns x 1 x Nt x Ntx x M x Nmcp matrix)
            s = cell(1, obj.network.numberOfActiveReceivingNodes);
            a = obj.complexgaussian([1, obj.numberOfTargets, obj.numberOfTrialsParallel]);
            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                if obj.configuration.pathLoss
                    A = permute(10.^(.05*obj.receivedPowerFromScatterers_dBW(:, rxID, :, :)), [1 3 4 2]); % Ntx x Nt x Nmcp matrix
                else
                    A = 1;
                end
                switch obj.network.networkMode
                    case "multiStatic"
                        switch obj.spatialCoherency
                            case "deterministic"
                                A = A.*1;
                            case "coherent"
                                % For each RX-TX pair phases and amplitudes
                                % are known relatively
                                A = A.*a;
                            case "correlatedAmplitudeFluctuations" % Similar to swerling 1
                                % For each RX-TX pair phases are independent
                                % uniformly distributed, amplitude are known relatively
                                A = A.*a.*exp(1j*2*pi*rand([obj.numberOfTransmittingNodes, obj.numberOfTargets, obj.numberOfTrialsParallel]));
                            case "incoherent" % Similar to swerling 2
                                % For each target, and RX-TX pair
                                % independent complex random variable
                                A = A.*obj.complexgaussian([obj.numberOfTransmittingNodes, obj.numberOfTargets, obj.numberOfTrialsParallel]);
                        end
                    case "monoStatic"
                        monoStaticTXIDs = obj.network.monoStaticTransmitterIDs(rxID);
                        nonMonoStaticTXIDs = setdiff(1 : obj.numberOfTransmittingNodes, monoStaticTXIDs);
                        switch obj.spatialCoherency
                            case "deterministic"
                                A(nonMonoStaticTXIDs, :, :) = A(nonMonoStaticTXIDs, :, :).*obj.complexgaussian([length(nonMonoStaticTXIDs), obj.numberOfTargets, obj.numberOfTrialsParallel]);
                            case {"coherent", "correlatedAmplitudeFluctuations"}
                                A(nonMonoStaticTXIDs, :, :) = A(nonMonoStaticTXIDs, :, :).*obj.complexgaussian([length(nonMonoStaticTXIDs), obj.numberOfTargets, obj.numberOfTrialsParallel]);
                                A(monoStaticTXIDs, :, :) = A(monoStaticTXIDs, :, :).*exp(1j*2*pi*rand([length(monoStaticTXIDs), obj.numberOfTargets, obj.numberOfTrialsParallel]));
                            case "incoherent"
                                A = A.*obj.complexgaussian([obj.numberOfTransmittingNodes, obj.numberOfTargets, obj.numberOfTrialsParallel]);
                        end
                end
                Aw = permute(A, [4 5 2 1 3]).*obj.waveformReceivedFromScatterers(:, rxID, :, :, :); % Ns x 1 x Nt x Ntx x Nmcp matrix
                M = obj.network.activeReceivingNodes(rxID).array.numberOfTotalElements;
                Ns = obj.network.activeReceivingNodes(rxID).numberOfSamplesPerCPI;
                if obj.configuration.noise
                    n = 10.^(.05*obj.network.activeReceivingNodes(rxID).noisePowerPerSample_dB).*obj.complexgaussian([Ns, 1, 1, 1, M, obj.numberOfTrialsParallel]);
                else
                    n = 0;
                end
                s{rxID} = zeros(Ns, 1, obj.numberOfTargets, obj.network.numberOfActiveTransmittingNodes, M, obj.numberOfTrialsParallel);
                for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                    s{rxID}(:, :, :, txID, :, :) = permute(Aw(:, :, :, txID, :), [1 2 3 4 6 5]).*permute(obj.receiveSteeringVector{txID, rxID}, [6 4 2 5 1 3]) + n;
                end
            end
        end

        function s = get.signalReceivedFromTransmittingNodes(obj)
            % (1 x Nrx cell of Ns x 1 x Ntx x M matrix)
            s = cell(1, obj.network.numberOfActiveReceivingNodes);
            for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                if obj.configuration.pathLoss
                    A = permute(10.^(.05*obj.receivedPowerFromTransmittingNodes_dBW(:, rxID, :)), [3, 2, 1]); % 1 x 1 x Ntx matrix
                else
                    A = 1;
                end
                Aw = A.*obj.waveformReceivedFromTransmitters(:, rxID, :); % Ns x 1 x Ntx matrix
                Aw(:, :, isinf(A)) = 0;
                M = obj.network.activeReceivingNodes(rxID).array.numberOfTotalElements;
                Ns = obj.network.activeReceivingNodes(rxID).numberOfSamplesPerCPI;
                s{rxID} = zeros(Ns, 1, obj.network.numberOfActiveTransmittingNodes, M);
                for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                    s{rxID}(:, :, txID, :) = Aw(:, :, txID).*permute(obj.receiveSteeringVectorDirectPath{txID, rxID}, [3 4 2 1]);
                end
            end
        end

        function s = get.signalSuperposed(obj)
            % (1 x Nrx cell of Ns x M x Nmcp matrix)
            s = cell(1, obj.network.numberOfActiveReceivingNodes);
            for rxID = 1 : obj.numberOfReceivingNodes
                switch obj.network.networkMode
                    case "multiStatic"
                        if obj.configuration.directPath
                            s{rxID} = permute(sum(obj.signalReceivedFromScatterers{rxID} + permute(obj.signalReceivedFromTransmittingNodes{rxID}, [1 2 5 3 4]), [3, 4]), [1 5 6 3 4 2]);
                        else
                            s{rxID} = permute(sum(obj.signalReceivedFromScatterers{rxID}, [3, 4]), [1 5 6 3 4 2]);
                        end
                    case "monoStatic"
                        x = obj.signalReceivedFromScatterers{rxID};
                        s{rxID} = zeros(size(x, 1), size(x, 5), size(x, 6));
                        %%% other TX signals not implemented
                        s{rxID} = permute(sum(x(:, :, :, obj.network.monoStaticTransmitterIDs(rxID), :, :), 3), [1 5 6 3 4 2]);
                end
            end
        end

        function s = get.signalBeamformed(obj)
            % (1 x Nrx cell of Ns x Nch x Nmcp matrix)
            switch obj.network.beamformingMode
                case 'conventional'
                    s = cell(1, obj.network.numberOfActiveReceivingNodes);
                    for rxID = 1 : obj.numberOfReceivingNodes
                        x = obj.signalSuperposed{rxID}; % Ns x M
                        switch obj.network.networkMode
                            case "multiStatic"
                                %%% Different frequency is not implemented
                                lambda = obj.network.activeTransmittingNodes(1).carrierWavelength;
                            case "monoStatic"
                                lambda = obj.network.activeTransmittingNodes(obj.network.monoStaticTransmitterIDs(rxID)).carrierWavelength;
                        end
                        beamformingCenters = exp(1j*2*pi*obj.network.activeReceivingNodes(rxID).beamCenterPositions./lambda);
                        as = exp(1j*2*pi*obj.network.activeReceivingNodes(rxID).array.steeringPositions./lambda);
                        s{rxID} = pagemtimes(x, conj(obj.network.activeReceivingNodes(rxID).beamformer.*beamformingCenters.*as));
                    end
                case 'bypass'
                    s = obj.signalSuperposed;
            end
        end

        %%% ellipse

        function a = get.majorAxis(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            a = obj.bistaticRange;
        end

        function b = get.minorAxis(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            b = sqrt(obj.majorAxis.^2 - obj.distanceBaseline.^2);
        end

        function e = get.eccentricity(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            e = obj.distanceBaseline./obj.bistaticRange;
        end

        function l = get.latusRectum(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            l = obj.majorAxis.*(1 - obj.eccentricity.^2);
        end

        %%% unit vectors

        function u = get.unitDirectionRX(obj) % unit vector from receiving nodes to target
            % (3 x Nrx x Nt x Nmcp matrix)
            if isempty(obj.targets)
                u = zeros(3, obj.numberOfReceivingNodes);
            else
                u = permute(obj.positions, [1 4 2 3]) - [obj.network.activeReceivingNodes.position];
                u = u./sqrt(sum(abs(u).^2));
                u(isnan(u)) = 0;
            end
        end

        function u = get.unitDirectionTX(obj) % unit vector from transmitting nodes to target
            % (3 x Ntx x Nt x Nmcp matrix)
            if isempty(obj.targets)
                u = zeros(3, obj.numberOfTransmittingNodes);
            else
                u = permute(obj.positions, [1 4 2 3]) - [obj.network.activeTransmittingNodes.position];
                u = u./sqrt(sum(abs(u).^2));
                u(isnan(u)) = 0;
            end
        end

        function g = get.unitDirectionBisector(obj)
            % (3 x Ntx x Nrx x Nt x Nmcp matrix)
            g = permute(obj.unitDirectionRX, [1 5 2 3 4]) + permute(obj.unitDirectionTX, [1 2 5 3 4]);
            g = g./sqrt(sum(abs(g).^2, 1));
            g(isnan(g)) = 0;
        end

        function t = get.unitDirectionBaseline(obj) % unit vector from transmitting nodes to receiving nodes
            % (3 x Ntx x Nrx matrix)
            t = obj.network.positionReceivingNode - obj.network.positionTransmittingNode; % baseline vector
            t = t./sqrt(sum(abs(t).^2, 1));
            t(isnan(t)) = 0;
        end

        function n = get.unitNormalBaseline(obj)
            % (3 x Ntx x Nrx x Nt matrix)
            t = repmat(obj.unitDirectionBaseline, [1 1 1 obj.numberOfTargets, obj.numberOfTrialsParallel]); % baseline vector
            u = repmat(permute(obj.unitDirectionRX, [1 5 2 3 4]), 1, obj.numberOfTransmittingNodes); % rx-target vector
            n = u - dot(t, u).*t;
            n = n./sqrt(sum(abs(n).^2, 1));
            n(isnan(n)) = 0;
        end

        %%% distances

        function R = get.distanceRX(obj)
            % (1 x Nrx x Nt x Nmcp matrix)
            if isempty(obj.targets)
                R = -inf(1, obj.numberOfReceivingNodes, obj.numberOfTransmittingNodes);
            else
                R = sqrt(sum((permute(obj.positions, [1 4 2 3]) - [obj.network.activeReceivingNodes.position]).^2));
            end
        end

        function R = get.distanceTX(obj)
            % (Ntx x 1 x Nt x Nmcp matrix)
            if isempty(obj.targets)
                R = -inf(1, obj.numberOfReceivingNodes, obj.numberOfTransmittingNodes);
            else
                R = permute(sqrt(sum((permute(obj.positions, [1 4 2 3]) - [obj.network.activeTransmittingNodes.position]).^2)), [2 1 3 4]);
            end
        end

        function d = get.distanceBisector(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            if isempty(obj.targets)
                d = -inf(obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
            else
                m = repmat(obj.network.positionReceivingNode + obj.network.positionTransmittingNode, [1, 1, 1, obj.numberOfTargets])/2;
                d = permute(sqrt(sum(abs(m - permute(obj.positions, [1 4 5 2 3])).^2)), [2 3 4 5 1]);
            end
        end

        function d = get.distanceBaseline(obj)
            % (Ntx x Nrx x Nt matrix)
            d = obj.network.distanceBaseline;
        end

        %%% angles

        function ksi = get.bistaticAngle(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            ksi = permute(real(acosd(pagemtimes(permute(obj.unitDirectionTX, [2 1 5 3 4]), permute(obj.unitDirectionRX, [1 5 2 3 4])))), [1 3 4 5 2]);
            %ksi = obj.localAngleBaselineRX + obj.localAngleBaselineTX;
            %ksi - 90 needed sometimes
        end

        function theta = get.elevationTargetRX(obj)
            % (Nrx x Nt x Nmcp matrix)
            theta = asind(shiftdim(obj.unitDirectionRX(3, :, :, :), 1));
        end

        function phi = get.azimuthTargetRX(obj)
            % (Nrx x Nt x Nmcp matrix)
            phiSign = sign(asind(shiftdim(obj.unitDirectionRX(2, :, :, :), 1)./cosd(obj.elevationTargetRX)));
            phiSign(isnan(phiSign)) = 1;
            phiSign(~phiSign) = 1;
            phi = phiSign.*real(acosd(shiftdim(obj.unitDirectionRX(1, :, :, :), 1)./cosd(obj.elevationTargetRX)));
        end

        function theta = get.elevationTargetTX(obj)
            % (Ntx x Nt x Nmcp matrix)
            theta = asind(shiftdim(obj.unitDirectionTX(3, :, :, :), 1));
        end

        function phi = get.azimuthTargetTX(obj)
            % (Ntx x Nt x Nmcp matrix)
            phiSign = sign(asind(shiftdim(obj.unitDirectionTX(2, :, :, :), 1)./cosd(obj.elevationTargetTX)));
            phiSign(isnan(phiSign)) = 1;
            phiSign(~phiSign) = 1;
            phi = phiSign.*real(acosd(shiftdim(obj.unitDirectionTX(1, :, :, :), 1)./cosd(obj.elevationTargetTX)));
        end

        function ksi = get.angleVelocityBisector(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            if isempty(obj.targets)
                ksi = nan(obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
            else
                ksi = permute(real(acosd(pagemtimes(permute(obj.unitDirectionBisector, [2 1 3 4 5]), permute(obj.targets.velocityUnitDirection, [1 3 4 2])))), [1 3 4 5 2]);
            end
        end

        function ksi = get.angleBaselineRX(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            ksi = permute(real(acosd(pagemtimes(permute(obj.unitNormalBaseline, [2 1 4 3 5]), permute(obj.unitDirectionRX, [1 5 3 2 4])))), [1 4 3 5 2]);
        end

        function ksi = get.angleBaselineTX(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            ksi = permute(real(acosd(pagemtimes(permute(obj.unitNormalBaseline, [3 1 4 2 5]), permute(obj.unitDirectionTX, [1 5 3 2 4])))), [4 1 3 5 2]);
        end

        %%% factors

        function F = get.triangulationFactor(obj)
            % (Ntx x Nrx x Nt matrix)
            F = sind(obj.bistaticAngle/2);
        end

        function D = get.lossFactor(obj)
            % (Ntx x Nrx x Nt matrix)
            D = sqrt(1 - obj.triangulationFactor.^2);
        end

        function region = get.targetRegion(obj)
            % (Ntx x Nrx x Nt matrix)
            D = obj.triangulationFactor;
            region = strings(size(D));
            region(D <= 0.1) = "baseline";
            region(D >= 0.9) = "lateral";
            region(D > 0.1 & D < 0.9) = "broadside";
        end

        %%% sizes

        function N = get.numberOfTargets(obj) % Nt
            if isempty(obj.targets)
                N = 1;
            else
                N = obj.targets.numberOfTargets;
            end
        end

        function N = get.numberOfReceivingNodes(obj) % Nrx
            N = obj.network.numberOfActiveReceivingNodes;
        end

        function N = get.numberOfTransmittingNodes(obj) % Ntx
            N = obj.network.numberOfActiveTransmittingNodes;
        end

        function N = get.numberOfBistaticPairs(obj) % Nrx*Ntx
            N = obj.network.numberOfActiveBistaticPairs;
        end

        %%% set methods

        function settargets(obj, targets)
            arguments
                obj
                targets (1, 1) target = obj.targets
            end
            obj.targets = targets;
            obj.settargetpositions;
        end

        function settargetpositions(obj, options)
            arguments
                obj
                options.width (:, 1) {mustBeNonnegative} = zeros(3, 1)
            end
            if any(options.width)
                obj.positions = obj.targets.position + options.width.*(randn(3, obj.numberOfTargets, obj.numberOfTrialsParallel) - 0.5);
            else
                obj.positions = obj.targets.position;
            end
        end

        %%% configuration

        function configure(obj, options)
            arguments
                obj
                options.noise (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.noise, [0, 1])} = obj.configuration.noise
                options.directPath (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.directPath, [0, 1])} = obj.configuration.directPath
                options.pathLoss (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.pathLoss, [0, 1])} = obj.configuration.pathLoss
                options.spatialCoherency (1, 1) string {mustBeMember(options.spatialCoherency, ["deterministic", "coherent", "correlatedAmplitudeFluctuations", "incoherent"])} = obj.spatialCoherency
            end
            obj.configuration.noise = options.noise;
            obj.configuration.directPath = options.directPath;
            obj.configuration.pathLoss = options.pathLoss;
            obj.spatialCoherency = options.spatialCoherency;
        end

        %%% visualization methods

        function visualizescenario(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.showVelocity (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.showVelocity, [0, 1])} = false
                options.showUnitNormalBaseline (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.showUnitNormalBaseline, [0, 1])} = false
                options.showUnitDirectionBisector (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.showUnitDirectionBisector, [0, 1])} = false
                options.showUnitNormalArray (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.showUnitNormalArray, [0, 1])} = false
                options.showSteeringDirection (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.showSteeringDirection, [0, 1])} = false
                options.showPattern (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.showPattern, [0, 1])} = false
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            posRX = [obj.network.activeReceivingNodes.position]/1e3;
            posTX = [obj.network.activeTransmittingNodes.position]/1e3;
            plot3(posRX(1, :), posRX(2, :), posRX(3, :), 'vb', 'LineWidth', 2, 'MarkerSize', 10);
            hold on;
            plot3(posTX(1, :), posTX(2, :), posTX(3, :), 'vr', 'LineWidth', 2, 'MarkerSize', 10);
            if ~isempty(obj.targets)
                x = squeeze(obj.positions(1, :, :))/1e3;
                y = squeeze(obj.positions(2, :, :))/1e3;
                z = squeeze(obj.positions(3, :, :))/1e3;
                vx = obj.targets.velocity(1, :)/1e3;
                vy = obj.targets.velocity(2, :)/1e3;
                vz = obj.targets.velocity(3, :)/1e3;
                if obj.numberOfTargets < 11
                    plot3(x, y, z, '+k', 'LineWidth', 2, 'MarkerSize', 10);
                elseif obj.numberOfTargets < 21
                    plot3(x, y, z, '*k', 'LineWidth', 1, 'MarkerSize', 10);
                else
                    plot3(x, y, z, '.k');
                end
                if options.showVelocity
                    quiver3(x, y, z, vx, vy, vz, 'Color', 'g', 'LineWidth', 2);
                end
                if options.showUnitNormalBaseline
                    n = reshape(obj.unitNormalBaseline, [3, obj.numberOfBistaticPairs, obj.numberOfTargets]);
                    xt = repmat(permute(x, [1 3 2]), [1, obj.numberOfBistaticPairs, 1]);
                    yt = repmat(permute(y, [1 3 2]), [1, obj.numberOfBistaticPairs, 1]);
                    zt = repmat(permute(z, [1 3 2]), [1, obj.numberOfBistaticPairs, 1]);
                    quiver3(xt, yt, zt, n(1, :, :), n(2, :, :), n(3, :, :), 'Color', 'm', 'LineWidth', 2);     
                end
                if options.showUnitDirectionBisector
                    n = reshape(obj.unitDirectionBisector, [3, obj.numberOfBistaticPairs, obj.numberOfTargets]);
                    xt = repmat(permute(x, [1 3 2]), [1, obj.numberOfBistaticPairs, 1]);
                    yt = repmat(permute(y, [1 3 2]), [1, obj.numberOfBistaticPairs, 1]);
                    zt = repmat(permute(z, [1 3 2]), [1, obj.numberOfBistaticPairs, 1]);
                    quiver3(xt, yt, zt, n(1, :, :), n(2, :, :), n(3, :, :), 'Color', 'c', 'LineWidth', 2); 
                end
            end
            text(posRX(1, :), posRX(2, :), posRX(3, :), num2str(find(obj.network.receivingNodeActivity).'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            text(posTX(1, :), posTX(2, :), posTX(3, :), num2str(find(obj.network.transmittingNodeActivity).'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            if obj.numberOfTargets < 11 && ~isempty(obj.targets)
                text(x, y, z, num2str((1 : obj.numberOfTargets).'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            end
            if options.showPattern
                if obj.numberOfTargets > 20
                    for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                        visualizationNormalization = max(obj.distanceTX(txID, :, :), [], 'all')/1e3;
                        G = visualizationNormalization + permute(abs(obj.transmittedBeam(txID, :, :, 1).^2./obj.network.activeTransmittingNodes(txID).array.numberOfTotalElements), [1 3 2]);
                        G(G < 0) = nan;
                        u = permute(obj.unitDirectionTX, [1 3 2 4]);
                        xu = u(1, :, txID, 1);
                        yu = u(2, :, txID, 1);
                        zu = u(3, :, txID, 1);
                        plot3(G.*xu + posTX(1, txID), G.*yu + posTX(2, txID), G.*zu + posTX(3, txID), '.r');
                        for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                            visualizationNormalization = max(obj.distanceRX(:, rxID, :, :), [], 'all')/1e3;
                            G = visualizationNormalization + permute(abs(obj.receivedBeamSpaceObservations{rxID}(txID, :, :, :, 1).^2./obj.network.activeReceivingNodes(rxID).array.numberOfTotalElements), [3 4 1 2]);
                            G(G < 0) = nan;
                            u = permute(obj.unitDirectionRX, [3 1 2 4]);
                            xu = u(:, 1, rxID, 1);
                            yu = u(:, 2, rxID, 1);
                            zu = u(:, 3, rxID, 1);
                            plot3(G.*xu + posRX(1, rxID), G.*yu + posRX(2, rxID), G.*zu + posRX(3, rxID), '.b');
                        end
                    end
                else
                    phi = linspace(-180, 180, 1000);
                    d = max(sqrt(sum(obj.network.boundary.^2)));
                    for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                        targetsTemp = target( ...
                            'position', [ ...
                            d*cosd(phi) + posTX(1, txID)*1e3; ...
                            d*sind(phi) + posTX(2, txID)*1e3; ...
                            zeros(1, 1000) + posTX(3, txID)*1e3]);
                        intTemp = interface( ...
                            'network', obj.network, ...
                            'targets', targetsTemp);
                        visualizationNormalization = max(intTemp.distanceTX(txID, :, :, :), [], 'all')/1e3;
                        G = visualizationNormalization.*permute(abs(intTemp.transmittedBeam(txID, :, :).^2./intTemp.network.activeTransmittingNodes(txID).array.numberOfTotalElements), [1 3 2]);
                        u = permute(intTemp.unitDirectionTX, [1 3 2 4]);
                        plot3(G.*u(1, :, txID, 1) + posTX(1, txID), G.*u(2, :, txID, 1) + posTX(2, txID), u(3, :, txID, 1) + posTX(3, txID), 'r');
                        for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                            targetsTemp = target( ...
                                'position', [ ...
                                d*cosd(phi) + posRX(1, rxID)*1e3; ...
                                d*sind(phi) + posRX(2, rxID)*1e3; ...
                                zeros(1, 1000) + posRX(3, rxID)*1e3]);
                            intTemp = interface( ...
                                'network', obj.network, ...
                                'targets', targetsTemp);
                            visualizationNormalization = max(intTemp.distanceRX(:, rxID, :, :), [], 'all')/1e3;
                            G = visualizationNormalization.*permute(abs(intTemp.receivedBeamSpaceObservations{rxID}(txID, :, :, :).^2./intTemp.network.activeReceivingNodes(rxID).array.numberOfTotalElements), [3 4 1 2]);
                            u = permute(intTemp.unitDirectionRX, [3 1 2 4]);
                            plot3(G.*u(:, 1, rxID, 1) + posRX(1, rxID), G.*u(:, 2, rxID, 1) + posRX(2, rxID), u(:, 3, rxID, 1) + posRX(3, rxID), 'b');
                        end
                    end
                end
            end
            if options.showUnitNormalArray
                for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                    n = obj.network.activeReceivingNodes(rxID).array.normalVector;
                    quiver3(posRX(1, rxID), posRX(2, rxID), posRX(3, rxID), n(1), n(2), n(3), 'Color', 'm', 'LineWidth', 2);
                end
                for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                    n = obj.network.activeTransmittingNodes(txID).array.normalVector;
                    quiver3(posTX(1, txID), posTX(2, txID), posTX(3, txID), n(1), n(2), n(3), 'Color', 'm', 'LineWidth', 2);
                end
            end
            if options.showSteeringDirection
                for rxID = 1 : obj.network.numberOfActiveReceivingNodes
                    s = obj.network.activeReceivingNodes(rxID).array.steeringUnitDirection;
                    quiver3(posRX(1, rxID), posRX(2, rxID), posRX(3, rxID), s(1), s(2), s(3), 'Color', 'm', 'LineWidth', 2, 'LineStyle', '--');
                end
                for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                    s = obj.network.activeTransmittingNodes(txID).array.steeringUnitDirection;
                    quiver3(posTX(1, txID), posTX(2, txID), posTX(3, txID), s(1), s(2), s(3), 'Color', 'm', 'LineWidth', 2, 'LineStyle', '--');
                end
            end
            posRX = repelem(posRX, 1, obj.numberOfTransmittingNodes);
            posTX = repmat(posTX, 1, obj.numberOfReceivingNodes);
            line([posRX(1, :); posTX(1, :)], [posRX(2, :); posTX(2, :)], [posRX(3, :); posTX(3, :)], 'lineStyle', '--', 'Color', 'k');
            grid off; grid on; grid minor; title('scenario');
            xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
            if ~isempty(obj.targets)
                legend('RX', 'TX', 'targets', 'Location', 'best');
            else
                legend('RX', 'TX', 'Location', 'best');
            end
            view(0, 90); hold off; drawnow;
        end

        function visualizereceivedsignals(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.receivingNodeIDs (1, :) double {mustBeInteger, mustBeNonnegative} = 1
                options.trialID (1, 1) double {mustBePositive, mustBeInteger} = 1
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            mustBeInRange(options.receivingNodeIDs, 1, obj.numberOfReceivingNodes);
            s = obj.signalSuperposed(options.receivingNodeIDs);
            for rxID = options.receivingNodeIDs
                t = obj.network.activeReceivingNodes(rxID).samplingInstants;
                if isscalar(options.receivingNodeIDs)
                    y = 20*log10(abs(s{rxID}(:, :, options.trialID)));
                else
                    y = 20*log10(abs(s{rxID}(:, ceil(end/2), options.trialID)));
                end
                y = y - obj.network.activeReceivingNodes(rxID).noisePowerPerSample_dB;
                plot(t, y); hold on;
            end
            grid off; grid on; grid minor;
            xlabel('time (s)'); ylabel('SNR (dB)');
            title('received signal');
            if ~isscalar(options.receivingNodeIDs)
                leg = legend(num2str(options.receivingNodeIDs.'), 'Location', 'best');
                title(leg, 'RX ID');
            end
            hold off; drawnow;
        end

        function visualizebeamformedsignals(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.receivingNodeIDs (1, :) double {mustBeInteger, mustBeNonnegative} = 1
                options.trialID (1, 1) double {mustBePositive, mustBeInteger} = 1
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            mustBeInRange(options.receivingNodeIDs, 1, obj.numberOfReceivingNodes);
            s = obj.signalBeamformed(options.receivingNodeIDs);
            for rxID = options.receivingNodeIDs
                t = obj.network.activeReceivingNodes(rxID).samplingInstants;
                if isscalar(options.receivingNodeIDs)
                    y = 20*log10(abs(s{rxID}(:, :, options.trialID)));
                else
                    y = 20*log10(abs(s{rxID}(:, ceil(end/2), options.trialID)));
                end
                y = y - obj.network.activeReceivingNodes(rxID).noisePowerPerSample_dB;
                plot(t, y); hold on;
            end
            grid off; grid on; grid minor;
            xlabel('time (s)'); ylabel('SNR (dB)');
            title('received signal');
            if ~isscalar(options.receivingNodeIDs)
                leg = legend(num2str(options.receivingNodeIDs.'), 'Location', 'best');
                title(leg, 'RX ID');
            end
            hold off; drawnow;
        end

        function visualizeisocurves(obj, options)
            arguments
                obj
                options.curves (1, :) string {mustBeMember(options.curves, ["delay", "power", "doppler", "bistaticAngle", "triangulation", "lossFactor", "all"])} = "delay"
                options.numberOfCurves (1, 1) double {mustBeInteger, mustBePositive} = 10
                options.transmittingNodeID (1, 1) double {mustBeInteger, mustBePositive} = 1
                options.receivingNodeID (1, 1) double {mustBeInteger, mustBePositive} = 1
            end
            mustBeInRange(options.transmittingNodeID, 1, obj.numberOfTransmittingNodes);
            mustBeInRange(options.receivingNodeID, 1, obj.numberOfReceivingNodes);
            if isempty(obj.targets)
                fprintf('There are no targets in the scene\n');
                return;
            end
            axID1 = obj.targets.firstAxisID;
            axID2 = obj.targets.secondAxisID;
            L1 = obj.targets.firstGridLength;
            L2 = obj.targets.secondGridLength;
            if L1 == 1 || L2 == 1
                warning('targets must be scanned as a grid');
                return;
            end
            switch axID1
                case 1
                    firstAxisLabel = 'x-pos (km)';
                case 2
                    firstAxisLabel = 'y-pos (km)';
                case 3
                    firstAxisLabel = 'z-pos (km)';
            end
            switch axID2
                case 1
                    secondAxisLabel = 'x-pos (km)';
                case 2
                    secondAxisLabel = 'y-pos (km)';
                case 3
                    secondAxisLabel = 'z-pos (km)';
            end
            posRX = obj.network.activeReceivingNodes(options.receivingNodeID).position/1e3;
            posTX = obj.network.activeTransmittingNodes(options.transmittingNodeID).position/1e3;
            x = reshape(obj.targets.position(axID1, :), L1, L2)/1e3;
            y = reshape(obj.targets.position(axID2, :), L1, L2)/1e3;
            if any(strcmpi(options.curves, "delay")) || any(strcmpi(options.curves, "all"))
                figure; contourf(x, y, reshape(-obj.timeDelay(options.transmittingNodeID, options.receivingNodeID, :, :)*1e6, L1, L2), options.numberOfCurves, "ShowText", true, "LabelFormat", "%3.3g usec");
                hold on; plot(posRX(axID1, :), posRX(axID2, :), 'ob', 'LineWidth', 2, 'MarkerSize', 10);
                plot(posTX(axID1, :), posTX(axID2, :), 'or', 'LineWidth', 2, 'MarkerSize', 10);
                xlabel(firstAxisLabel); ylabel(secondAxisLabel);
                xlim tight; ylim tight; zlim tight; view(2);
                title('iso delay curves'); drawnow;
            end
            if any(strcmpi(options.curves, "power")) || any(strcmpi(options.curves, "all"))
                figure; contourf(x, y, reshape(obj.inputSNR_dB(options.transmittingNodeID, options.receivingNodeID, :), L1, L2), options.numberOfCurves, "ShowText", true, "LabelFormat", "%3.3g dB");
                hold on; plot(posRX(axID1, :), posRX(axID2, :), 'ob', 'LineWidth', 2, 'MarkerSize', 10);
                plot(posTX(axID1, :), posTX(axID2, :), 'or', 'LineWidth', 2, 'MarkerSize', 10);
                xlabel(firstAxisLabel); ylabel(secondAxisLabel);
                xlim tight; ylim tight; zlim tight; view(2);
                title('iso power curves'); drawnow;
            end
            if any(strcmpi(options.curves, "doppler")) || any(strcmpi(options.curves, "all"))
                figure; contourf(x, y, reshape(obj.dopplerShift(options.transmittingNodeID, options.receivingNodeID, :), L1, L2), options.numberOfCurves, "ShowText", true, "LabelFormat", "%3.3g Hz");
                hold on; plot(posRX(axID1, :), posRX(axID2, :), 'ob', 'LineWidth', 2, 'MarkerSize', 10);
                plot(posTX(axID1, :), posTX(axID2, :), 'or', 'LineWidth', 2, 'MarkerSize', 10);
                xlabel(firstAxisLabel); ylabel(secondAxisLabel);
                xlim tight; ylim tight; zlim tight; view(2);
                title('iso doppler curves'); drawnow;
            end
            if any(strcmpi(options.curves, "bistaticAngle")) || any(strcmpi(options.curves, "all"))
                figure; contourf(x, y, reshape(obj.bistaticAngle(options.transmittingNodeID, options.receivingNodeID, :), L1, L2), linspace(0, 180, options.numberOfCurves + 3), "ShowText", true, "LabelFormat", "%3.3g");
                hold on; plot(posRX(axID1, :), posRX(axID2, :), 'ob', 'LineWidth', 2, 'MarkerSize', 10);
                plot(posTX(axID1, :), posTX(axID2, :), 'or', 'LineWidth', 2, 'MarkerSize', 10);
                xlabel(firstAxisLabel); ylabel(secondAxisLabel);
                xlim tight; ylim tight; zlim tight; view(2);
                title('iso bistatic angle curves'); drawnow;
            end
            if any(strcmpi(options.curves, "triangulation")) || any(strcmpi(options.curves, "all"))
                figure; contourf(x, y, reshape(obj.triangulationFactor(options.transmittingNodeID, options.receivingNodeID, :), L1, L2), linspace(0, 1, options.numberOfCurves + 1), "ShowText", true, "LabelFormat", "%3.3g");
                hold on; plot(posRX(axID1, :), posRX(axID2, :), 'ob', 'LineWidth', 2, 'MarkerSize', 10);
                plot(posTX(axID1, :), posTX(axID2, :), 'or', 'LineWidth', 2, 'MarkerSize', 10);
                xlabel(firstAxisLabel); ylabel(secondAxisLabel);
                xlim tight; ylim tight; zlim tight; view(2);
                title('iso triangulation curves'); drawnow;
            end
            if any(strcmpi(options.curves, "lossFactor")) || any(strcmpi(options.curves, "all"))
                figure; contourf(x, y, reshape(obj.lossFactor(options.transmittingNodeID, options.receivingNodeID, :), L1, L2), linspace(0, 1, options.numberOfCurves + 1), "ShowText", true, "LabelFormat", "%3.3g");
                hold on; plot(posRX(axID1, :), posRX(axID2, :), 'ob', 'LineWidth', 2, 'MarkerSize', 10);
                plot(posTX(axID1, :), posTX(axID2, :), 'or', 'LineWidth', 2, 'MarkerSize', 10);
                xlabel(firstAxisLabel); ylabel(secondAxisLabel);
                xlim tight; ylim tight; zlim tight; view(2);
                title('iso lossFactor curves'); drawnow;
            end
        end
    end

    methods(Static)
        function n = complexgaussian(sizeRV)
            if ~nargin
                sizeRV = 1;
            end
            n = (randn(sizeRV) + 1j*randn(sizeRV))/sqrt(2);
        end
    end
end