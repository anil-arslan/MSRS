classdef interface < handle
    %interface Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = ?fusionCenter, GetAccess = public)
        network (1, :) radarNetwork = radarNetwork.empty()
        targets (1, :) target = target.empty()
        environment (1, 1) %%%%% not implemented
        configuration (1, 1) struct = struct( ...
            'noise', false, ...
            'directPath', false, ...
            'pathLoss', true, ...
            'fractionalDelay', "off", ...
            'timeSynchronization', true, ...
            'frequencySynchronization', true ...
            )
        spatialCoherency (1, 1) string {mustBeMember(spatialCoherency, ["deterministic", "coherent", "correlatedAmplitudeFluctuations", "noncoherent"])} = "deterministic"
        swerling (1, 1) double {mustBeInteger, mustBeMember(swerling, 0 : 4)} = 0
        positions double % meters (3 x Nt x Nmcp matrix)
    end

    properties (Dependent)
        %%% target parameters
        carrierWavelength double % meters (Ntx x 1) beamforming
        carrierWavelengthRealized double % meters (Ntx x 1 x 1 x Nmcp) signal modeling
        bistaticRange double % meters (Ntx x Nrx x Nt x Nmcp matrix) major axis length
        timeDelay double % seconds (Ntx x Nrx x Nt x Nmcp matrix)
        dopplerShift double % Hz (Ntx x Nrx x Nt x Nmcp matrix)
        dopplerShiftBisector double % Hz (Ntx x Nrx x Nt x Nmcp matrix)

        transmitBackOfArrayTarget double % (Nt x Nmcp x Ntx matrix)
        transmitSteeringVector cell % (Ntx x 1 cell of M x Nt x Nmcp matrix)
        transmitBackOfArrayReceiver double % (Nrx x Ntx matrix)
        transmitSteeringVectorDirectPath cell % (Ntx x 1 cell of M x Nrx x Nmcp matrix)

        receiveBackOfArrayTarget double % (Nt x Nmcp x Nrx matrix)
        receiveSteeringVector cell % (Ntx x Nrx cell of M x Nt x Nmcp matrix)
        receiveBackOfArrayTransmitter double % (Ntx x Nrx matrix)
        receiveSteeringVectorDirectPath cell % (Ntx x Nrx cell of M x 1 x Nmcp matrix)

        transmittedBeam dobule % (Ntx x 1 x Nt x Nmcp matrix)
        transmittedBeamDirectPath dobule % (Ntx x 1 x Nrx x Nmcp matrix)
        receivedBeamSpaceObservations cell % (1 x Nrx cell of Ntx x 1 x Nt x Nrxch x Nmcp matrix)
        transmitGain_lin dobule % linear scale (Ntx x 1 x Nt x Nmcp matrix)
        transmitGainDirectPath_lin dobule % linear scale (Ntx x 1 x Nrx x Nmcp matrix)
        receiveGain_lin cell % linear scale % (1 x Nrx cell of Ntx x 1 x Nt x Nrxch x Nmcp matrix)

        receivedPowerFromScatterersMaximumGain_dBW double % dB Watt (Ntx x Nrx x Nt x Nmcp matrix)
        
        receivedPowerFromScatterers_dBW double % dB Watt (Ntx x Nrx x Nt x Nmcp matrix)
        receivedPowerFromTransmittingNodes_dBW double % dB Watt (Ntx x Nrx matrix)
        averageSNR_dB double % dB (Ntx x Nrx x Nt x Nmcp matrix)
        
        receivedEnergyFromScatterersPerPulse_dBJoule % dB Joule (Ntx x Nrx x Nt x Nmcp matrix)
        receivedEnergyFromScatterersPerSample_dBJoule % dB Joule (Ntx x Nrx x Nt x Nmcp matrix)
            %%% signal modeling
            waveformReceivedFromScatterers double % complex received waveform (Ns x Nrx x Nt x Ntx x Nmcp matrix)
            waveformReceivedFromTransmitters double % complex received waveform (Ns x Nrx x Ntx matrix)
            signalReceivedFromScatterers % complex received signal (1 x Nrx cell of Ns x 1 x Nt x Ntx x M x Nmcp matrix)
            fluctuation % complex fluctutation (Ntx x Nrx x Nt x Nmcp matrix)
            signalReceivedFromTransmittingNodes double % complex received direct path signal % (1 x Nrx cell of Ns x 1 x Ntx x M matrix)
            signalSuperposed % (1 x Nrx cell of Ns x M x Ntxch x Nmcp matrix)
            signalBeamformed % (1 x Nrx cell of Ns x Nrxch x Ntxch x Nmcp matrix)

        %%% ellipse
        majorAxis double % meters (Ntx x Nrx x Nt x Nmcp matrix)
        minorAxis double % meters (Ntx x Nrx x Nt x Nmcp matrix)
        eccentricity double % meters (Ntx x Nrx x Nt x Nmcp matrix)
        latusRectum double % meters (Ntx x Nrx x Nt x Nmcp matrix)
        ellipseRotationMatrix double % meters (3 x 3 x Ntx x Nrx matrix)
        ellipse double % (3 x Nscan x Ntx x Nrx x Nt x Nmcp matrix)
        ellipseResolution double % (3 x Nscan x Ntx x Nrx x Nt x Nmcp x 2 matrix)

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
        centerBaseline double % (3 x Ntx x Nrx matrix)
        azimuthBaseline double % degrees (Ntx x Nrx matrix)

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

    properties (Access = ?fusionCenter)
        numberOfTrialsParallel (1, 1) double {mustBeNonnegative, mustBeInteger} = 1
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

        %%% signal modeling important variables

        function lambda = get.carrierWavelength(obj) % beamforming
            % meters (Ntx x 1) beamforming
            lambda = obj.speedOfLight./[obj.network.activeTransmittingNodes.carrierFrequency].';
        end

        function lambda = get.carrierWavelengthRealized(obj) % signal modeling
            % meters (Ntx x 1 x 1 x Nmcp) signal modeling
            if ~obj.configuration.frequencySynchronization
                frequencyDeviation = zeros(obj.numberOfTransmittingNodes, 1, 1, obj.numberOfTrialsParallel);
                deviatedFrequencies = logical(obj.network.frequencyDeviationActiveBistaticPairs);
                frequencyDeviation(deviatedFrequencies, :, :, :) = exp(log(obj.network.frequencyDeviationActiveBistaticPairs(deviatedFrequencies)).*randn(nnz(deviatedFrequencies), 1, 1, obj.numberOfTrialsParallel));
                frequencyOffset = obj.network.frequencyOffsetActiveBistaticPairs;
                lambda = obj.speedOfLight./([obj.network.activeTransmittingNodes.carrierFrequency].' + frequencyOffset + frequencyDeviation);
            else
                lambda = obj.speedOfLight./[obj.network.activeTransmittingNodes.carrierFrequency].';
            end
        end

        function rho = get.bistaticRange(obj) % major axis
            % (Ntx x Nrx x Nt x Nmcp matrix)
            rho = obj.distanceRX + obj.distanceTX;
        end

        function tau = get.timeDelay(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            tau = obj.bistaticRange/obj.speedOfLight;
        end

        %%% steering vectors

        function backOfArray = get.transmitBackOfArrayTarget(obj)
            % (Nt x Nmcp x Ntx matrix)
            transmitArrays = [obj.network.activeTransmittingNodes.array];
            if isempty(transmitArrays)
                backOfArray = false(obj.numberOfTargets, 1);
            else
                backOfArray = shiftdim(pagemtimes(permute([transmitArrays.normalVector], [1 3 4 2]), 'ctranspose', permute(obj.unitDirectionTX, [1 3 4 2]), 'none') < shiftdim(cosd([transmitArrays.backOfArrayRegion]/2), -2), 1);
                backOfArray(:, :, [transmitArrays.backOfArray]) = false;
            end
        end

        function a = get.transmitSteeringVector(obj)
            % (Ntx x 1 cell of M x Nt x Nmcp matrix)
            transmitArrays = [obj.network.activeTransmittingNodes.array];
            backOfArray = obj.transmitBackOfArrayTarget;
            a = cell(obj.numberOfTransmittingNodes, 1);
            for txID = 1 : obj.numberOfTransmittingNodes
                % unit vector from radar to target
                a{txID} = exp(1j*2*pi*pagemtimes(transmitArrays(txID).positions, 'ctranspose', permute(obj.unitDirectionTX(:, txID, :, :), [1 3 4 2]), 'none')./permute(obj.carrierWavelengthRealized(txID, :, :, :), [1 2 4 3]));
                if ~transmitArrays(txID).backOfArray
                    a{txID}(:, backOfArray(:, :, txID)) = 0;
                end
            end
        end

        function backOfArray = get.transmitBackOfArrayReceiver(obj)
            % (Nrx x Ntx matrix)
            transmitArrays = [obj.network.activeTransmittingNodes.array];
            if isempty(transmitArrays)
                backOfArray = false(obj.numberOfReceivingNodes, 1);
            else
                backOfArray = shiftdim(pagemtimes(permute([transmitArrays.normalVector], [1 3 2]), 'ctranspose', permute(obj.unitDirectionBaseline, [1 3 2]), 'none') < cosd([transmitArrays.backOfArrayRegion]/2), 1);
                backOfArray(:, [transmitArrays.backOfArray]) = false;
            end
        end

        function a = get.transmitSteeringVectorDirectPath(obj)
            % (Ntx x 1 cell of M x Nrx x Nmcp matrix)
            transmitArrays = [obj.network.activeTransmittingNodes.array];
            backOfArray = obj.transmitBackOfArrayReceiver;
            a = cell(obj.numberOfTransmittingNodes, 1);
            for txID = 1 : obj.numberOfTransmittingNodes
                % unit vector from transmitting node to receiving node
                a{txID} = exp(1j*2*pi*transmitArrays(txID).positions'*permute(obj.unitDirectionBaseline(:, txID, :), [1 3 2])./permute(obj.carrierWavelengthRealized(txID, :, :, :), [1 2 4 3]));
                if ~transmitArrays(txID).backOfArray
                    a{txID}(:, backOfArray(:, txID)) = 0;
                end
            end
        end

        function backOfArray = get.receiveBackOfArrayTarget(obj)
            % (Nt x Nmcp x Nrx matrix)
            receiveArrays = [obj.network.activeReceivingNodes.array];
            if isempty(receiveArrays)
                backOfArray = false(obj.numberOfTargets, 1);
            else
                backOfArray = shiftdim(pagemtimes(permute([receiveArrays.normalVector], [1 3 4 2]), 'ctranspose', permute(obj.unitDirectionRX, [1 3 4 2]), 'none') < shiftdim(cosd([receiveArrays.backOfArrayRegion]/2), -2), 1);
                backOfArray(:, :, [receiveArrays.backOfArray]) = false;
            end
        end

        function a = get.receiveSteeringVector(obj)
            % (Ntx x Nrx cell of M x Nt x Nmcp matrix)
            receiveArrays = [obj.network.activeReceivingNodes.array];
            backOfArray = obj.receiveBackOfArrayTarget;
            a = cell(obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
            for rxID = 1 : obj.numberOfReceivingNodes
                for txID = 1 : obj.numberOfTransmittingNodes
                    a{txID, rxID} = exp(1j*2*pi*pagemtimes(receiveArrays(rxID).positions, 'ctranspose', permute(obj.unitDirectionRX(:, rxID, :, :), [1 3 4 2]), 'none')./permute(obj.carrierWavelengthRealized(txID, :, :, :), [1 2 4 3]));
                    if ~receiveArrays(rxID).backOfArray
                        a{txID, rxID}(:, backOfArray(:, :, rxID)) = 0;
                    end
                end
            end
        end

        function backOfArray = get.receiveBackOfArrayTransmitter(obj)
            % (Ntx x Nrx matrix)
            receiveArrays = [obj.network.activeReceivingNodes.array];
            if isempty(receiveArrays)
                backOfArray = false(obj.numberOfTransmittingNodes, 1);
            else
                backOfArray = shiftdim(pagemtimes(permute([receiveArrays.normalVector], [1 3 2]), 'ctranspose', obj.unitDirectionBaseline, 'none') < cosd([receiveArrays.backOfArrayRegion]/2), 1);
                backOfArray(:, [receiveArrays.backOfArray]) = false;
            end
        end

        function a = get.receiveSteeringVectorDirectPath(obj)
            % (Ntx x Nrx cell of M x 1 x Nmcp matrix)
            receiveArrays = [obj.network.activeReceivingNodes.array];
            backOfArray = obj.receiveBackOfArrayTransmitter;
            a = cell(obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
            for rxID = 1 : obj.numberOfReceivingNodes
                for txID = 1 : obj.numberOfTransmittingNodes
                    a{txID, rxID} = exp(1j*2*pi*receiveArrays(rxID).positions'*obj.unitDirectionBaseline(:, txID, rxID)./permute(obj.carrierWavelengthRealized(txID, :, :, :), [1 2 4 3]));
                    if ~receiveArrays(rxID).backOfArray
                        a{txID, rxID}(:, backOfArray(txID, rxID)) = 0;
                    end
                end
            end
        end

        %%% beams with signal gains close to real directivity

        function y = get.transmittedBeam(obj)
            % (Ntx x 1 x Nt x Nmcp matrix)
            y = zeros(obj.numberOfTransmittingNodes, 1, obj.numberOfTargets, obj.numberOfTrialsParallel);
            for txID = 1 : obj.numberOfTransmittingNodes
                beam = obj.network.activeTransmittingNodes(txID).beamform(permute(obj.transmitSteeringVector{txID}, [4 1 2 3]));
                if size(beam, 4) == obj.numberOfTrialsParallel
                    y(txID, 1, :, :) = beam;
                else
                    y(txID, 1, :, :) = repmat(beam, [1 1 1 obj.numberOfTrialsParallel]);
                end
            end
        end

        function y = get.transmittedBeamDirectPath(obj)
            % (Ntx x 1 x Nrx x Nmcp matrix)
            y = zeros(obj.numberOfTransmittingNodes, 1, obj.numberOfReceivingNodes, obj.numberOfTrialsParallel);
            for txID = 1 : obj.numberOfTransmittingNodes
                beam = obj.network.activeTransmittingNodes(txID).beamform(permute(obj.transmitSteeringVectorDirectPath{txID}, [4 1 2 3]));
                if size(beam, 4) == obj.numberOfTrialsParallel
                    y(txID, 1, :, :) = beam;
                else
                    y(txID, 1, :, :) = repmat(beam, [1 1 1 obj.numberOfTrialsParallel]);
                end
            end
        end

        function y = get.receivedBeamSpaceObservations(obj)
            % (1 x Nrx cell of Ntx x 1 x Nt x Nrxch x Nmcp matrix)
            y = cell(1, obj.numberOfReceivingNodes);
            for rxID = 1 : obj.numberOfReceivingNodes
                y{rxID} = zeros(obj.numberOfTransmittingNodes, 1, obj.numberOfTargets, obj.network.activeReceivingNodes(rxID).numberOfTotalChannels, obj.numberOfTrialsParallel);
                for txID = 1 : obj.numberOfTransmittingNodes
                    % M --> Nrxch (bypass: Nrxch = M)
                    beam = obj.network.activeReceivingNodes(rxID).beamform(permute(obj.receiveSteeringVector{txID, rxID}, [4 1 2 3]), obj.carrierWavelength(txID));
                    if size(beam, 5) == obj.numberOfTrialsParallel
                        y{rxID}(txID, 1, :, :, :) = beam;
                    else
                        y{rxID}(txID, 1, :, :, :) = repmat(beam, [1 1 1 obj.numberOfTrialsParallel]);
                    end
                end
            end
        end

        %%% to simulate directivity

        function Gt = get.transmitGain_lin(obj)
            % (Ntx x 1 x Nt x Nmcp)
            arrayNodesTX = [obj.network.activeTransmittingNodes.array];
            if isempty(arrayNodesTX)
                Gt = zeros(1, 1, obj.numberOfTargets);
            else
                % GtMax = sqrt([obj.network.activeTransmittingNodes.peakTransmitGain].');
                GtMax = sqrt(4*pi*[arrayNodesTX.apertureArea]./obj.carrierWavelengthRealized.^2);
                Gt = obj.transmittedBeam.*GtMax./sqrt([arrayNodesTX.numberOfTotalElements]).';
            end
        end

        function Gt = get.transmitGainDirectPath_lin(obj)
            % (Ntx x 1 x Nrx x Nmcp)
            arrayNodesTX = [obj.network.activeTransmittingNodes.array];
            if isempty(arrayNodesTX)
                Gt = zeros(1, 1, obj.numberOfReceivingNodes);
            else
                % GtMax = sqrt([obj.network.activeTransmittingNodes.peakTransmitGain].');
                GtMax = sqrt(4*pi*[arrayNodesTX.apertureArea]./obj.carrierWavelengthRealized.^2);
                Gt = obj.transmittedBeamDirectPath.*GtMax./sqrt([arrayNodesTX.numberOfTotalElements]).';
            end
        end

        function Gr = get.receiveGain_lin(obj)
            % (1 x Nrx cell of Ntx x 1 x Nt x Nrxch x Nmcp matrix)
            arrayNodesRX = [obj.network.activeReceivingNodes.array];
            if isempty(arrayNodesRX)
                Gr = {zeros(obj.numberOfTransmittingNodes, 1, obj.numberOfTargets)};
            else
                GrMax = sqrt(4*pi*[arrayNodesRX.apertureArea]./permute(obj.carrierWavelengthRealized, [1 2 3 5 4]).^2);
                Gr = cell(1, obj.numberOfReceivingNodes);
                for rxID = 1 : obj.numberOfReceivingNodes
                    Gr{rxID} = obj.receivedBeamSpaceObservations{rxID}.*GrMax./sqrt(arrayNodesRX(rxID).numberOfTotalElements);
                end
            end
        end

        %%% power and energy

        function Pr = get.receivedPowerFromScatterersMaximumGain_dBW(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            if isempty(obj.targets)
                Pr = -inf(obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
            else
                powerGains = 10*log10([obj.network.activeTransmittingNodes.peakTransmitGain].') + ...
                    10*log10([obj.network.activeTransmittingNodes.peakPower_W]./[obj.network.activeTransmittingNodes.pulseWidth]).' + permute(obj.targets.RCS_dbms, [1 3 2]) + ...
                    20*log10(obj.carrierWavelengthRealized);
                powerLosses = [obj.network.activeReceivingNodes.systemLoss_dB] + 30*log10(4*pi) + 20*log10(obj.distanceRX.*obj.distanceTX);
                Pr = powerGains - powerLosses;
                Pr(isinf(Pr) & Pr > 0) = 20*log10(1/eps);
            end
        end

        function Pr = get.receivedPowerFromScatterers_dBW(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            if isempty(obj.targets)
                Pr = -inf(obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
            else
                powerGains = 20*log10(abs(obj.transmittedBeam)) + ...
                    10*log10([obj.network.activeTransmittingNodes.peakPower_W]).' + permute(obj.targets.RCS_dbms, [1 3 2]) + ...
                    20*log10(obj.carrierWavelengthRealized);
                powerLosses = [obj.network.activeReceivingNodes.systemLoss_dB] + 30*log10(4*pi) + 20*log10(obj.distanceRX.*obj.distanceTX);
                Pr = powerGains - powerLosses;
                Pr(isinf(Pr) & Pr > 0) = 20*log10(1/eps);
            end
            CPIcontrol = obj.timeDelay - [obj.network.activeReceivingNodes.CPI] - [obj.network.activeReceivingNodes.samplingPeriod] > 0;
            if size(CPIcontrol, 4) == 1
                CPIcontrol = repmat(CPIcontrol, [1 1 1 size(Pr, 4)]);
            end
            Pr(CPIcontrol) = -inf;
            if ~strcmpi(obj.network.surveillanceMode, "rotating")
                backOfArrayControl = permute(obj.transmitBackOfArrayTarget, [3 4 1 2]) | permute(obj.receiveBackOfArrayTarget, [4 3 1 2]);
                if size(backOfArrayControl, 4) == 1
                    backOfArrayControl = repmat(backOfArrayControl, [1 1 1 size(Pr, 4)]);
                end
            else
                backOfArrayControl = false(size(Pr));
            end
            Pr(backOfArrayControl) = -inf;
        end

        function Pr = get.receivedPowerFromTransmittingNodes_dBW(obj)
            % (Ntx x Nrx x 1 x Nmcp matrix)
            powerGains = permute(20*log10(abs(obj.transmittedBeamDirectPath)), [1 3 2 4]) + ...
                10*log10([obj.network.activeTransmittingNodes.peakPower_W].') + ... % no RCS term
                20*log10(obj.carrierWavelengthRealized);
            powerLosses = [obj.network.activeReceivingNodes.systemLoss_dB] + 20*log10(4*pi*obj.distanceBaseline);
            Pr = powerGains - powerLosses;
        end

        function SNRin = get.averageSNR_dB(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            SNRin = obj.receivedPowerFromScatterers_dBW - [obj.network.activeReceivingNodes.noisePowerPerSample_dB];
        end

        function Er = get.receivedEnergyFromScatterersPerSample_dBJoule(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            Er = obj.receivedPowerFromScatterers_dBW./[obj.network.activeReceivingNodes.samplingFrequency];
        end

        function Er = get.receivedEnergyFromScatterersPerPulse_dBJoule(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            Er = obj.receivedEnergyFromScatterersPerSample_dBJoule.*obj.network.pulseWidthSample.';
        end

        %%% doppler shift

        function fd = get.dopplerShift(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            if isempty(obj.targets)
                fd = zeros(obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
            else
                vdRx = pagemtimes(permute(obj.targets.velocity, [3 1 5 2 4]), permute(obj.unitDirectionRX, [1 5 2 3 4]));
                vdTx = pagemtimes(permute(obj.targets.velocity, [3 1 5 2 4]), permute(obj.unitDirectionTX, [1 2 5 3 4]));
                fd = -shiftdim((vdRx + vdTx), 1)./obj.carrierWavelengthRealized;
            end
        end

        function fd = get.dopplerShiftBisector(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            if isempty(obj.targets)
                fd = zeros(obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
            else
                fd = 2*shiftdim(obj.targets.velocityMagnitude, -1).*cosd(obj.bistaticAngle/2).*cosd(obj.angleVelocityBisector)./obj.carrierWavelengthRealized;
            end
        end

        %%% signal modeling

        function waveforms = get.waveformReceivedFromScatterers(obj)
            % (Ns x Nrx x Nt x Ntx x Nmcp matrix)
            Ts = [obj.network.activeReceivingNodes.samplingPeriod]; % 1 x Nrx matrix
            tau = [obj.network.activeReceivingNodes.samplingInstants] - permute(obj.timeDelay, [5 2 3 1 4]); % Ns x Nrx x Nt x Ntx x Nmcp matrix
            if ~obj.configuration.timeSynchronization
                tau = tau + permute(obj.network.timeOffsetActiveBistaticPairs, [3 2 4 1]) + permute(obj.network.timeDeviationActiveBistaticPairs, [3 2 4 1]).*randn(1, obj.numberOfReceivingNodes, 1, obj.numberOfTransmittingNodes, obj.numberOfTrialsParallel);
            end
            fd = permute(obj.dopplerShift, [5 2 3 1 4]); % 1 x Nrx x Nt x Ntx matrix
            fc = shiftdim([obj.network.activeTransmittingNodes.carrierFrequency], -2); % 1 x 1 x 1 x Ntx vector
            modulatorSpatial = exp(-1j*2*pi*fc.*tau(1, :, :, :, :)); % Ns x Nrx x Nt x Ntx x Nmcp matrix
            modulatorDoppler = exp(-1j*2*pi*fd.*tau);
            % across time it is very small.
            waveforms = zeros(size(tau)); % Ns x Nrx x Nt x Ntx matrix
            for txID = 1 : obj.numberOfTransmittingNodes
                waveforms(:, :, :, txID, :) = obj.network.activeTransmittingNodes(txID).waveform(-tau(:, :, :, txID, :), unique(Ts));
            end
            waveforms = waveforms.*modulatorDoppler.*modulatorSpatial;
            waveforms(isinf(tau)) = 0;
            waveforms = waveforms.*sqrt([obj.network.activeTransmittingNodes.pulseWidth]./[obj.network.activeReceivingNodes.samplingPeriod]);

            switch obj.configuration.fractionalDelay
                case 'sinc-based'
                    notAlignedSampling = tau(1, :, :, :, :).*[obj.network.activeReceivingNodes.samplingFrequency];
                    notAlignedSampling = notAlignedSampling - round(notAlignedSampling);
        
                    % Sinc-based fractional delay filter
                    N = unique([obj.network.pulseWidthSample]);  % Number of neighboring samples to consider
                    n = (0 : N - 1).';  % Filter taps
                    h = sinc(n - notAlignedSampling - ceil(N/2));  % Fractional delay filter
        
                    % Create sampled signal and apply delay filter
                    for waveformID = 1 : prod(size(waveforms, [2 3 4 5]))
                        waveforms(:, waveformID) = conv(waveforms(:, waveformID), h(:, waveformID), 'same');  % Apply fractional delay filter
                    end
                case 'lagrange-based'
                case 'off'
            end
        end

        function waveforms = get.waveformReceivedFromTransmitters(obj)
            % (Ns x Nrx x Ntx matrix)
            Ts = [obj.network.activeReceivingNodes.samplingPeriod]; % 1 x Nrx matrix
            tau = [obj.network.activeReceivingNodes.samplingInstants] - permute(obj.network.directPathDelay, [3 2 1]); % Ns x Nrx x Ntx matrix
            if ~obj.configuration.timeSynchronization
                %%% not implemented
            end
            fc = shiftdim([obj.network.activeTransmittingNodes.carrierFrequency], -1); % 1 x 1 x Ntx vector
            modulatorSpatial = exp(-1j*2*pi*fc.*tau(1, :, :)); % Ns x Nrx x Ntx matrix
            waveforms = zeros(size(tau)); % Ns x Nrx x Ntx matrix
            for txID = 1 : obj.numberOfTransmittingNodes
                waveforms(:, :, txID) = obj.network.activeTransmittingNodes(txID).waveform(tau(:, :, txID), unique(Ts));
            end
            waveforms = waveforms.*modulatorSpatial;
            waveforms(isinf(tau)) = 0;

            switch obj.configuration.fractionalDelay
                case 'sinc-based'
                    notAlignedSampling = tau(1, :, :).*[obj.network.activeReceivingNodes.samplingFrequency];
                    notAlignedSampling = notAlignedSampling - round(notAlignedSampling);
        
                    % Sinc-based fractional delay filter
                    N = 10;  % Number of neighboring samples to consider
                    n = (-N : N).';  % Filter taps
                    h = sinc(n - notAlignedSampling);  % Fractional delay filter
        
                    % Create sampled signal and apply delay filter
                    for waveformID = 1 : prod(size(waveforms, [2 3 4 5]))
                        waveforms(:, waveformID) = conv(waveforms(:, waveformID), h(:, waveformID), 'same');  % Apply fractional delay filter
                    end
                case 'lagrange-based'
                case 'off'
            end
        end

        function y = get.signalReceivedFromScatterers(obj)
            % (1 x Nrx cell of Ns x 1 x Nt x Ntx x M x Nmcp matrix)
            y = cell(1, obj.numberOfReceivingNodes);
            targetComplexFluctuation = obj.fluctuation;
            if obj.configuration.pathLoss
                Pr = permute(10.^(.05*obj.receivedPowerFromScatterers_dBW), [1 3 4 2]); % Ntx x Nt x Nmcp x Nrx matrix
            else
                Pr = ones(obj.numberOfTransmittingNodes, obj.numberOfTargets, obj.numberOfTrialsParallel, obj.numberOfReceivingNodes);
            end
            % a = obj.receiveSteeringVector;
            s = obj.waveformReceivedFromScatterers;
            for rxID = 1 : obj.numberOfReceivingNodes
                alpha = permute(targetComplexFluctuation(:, rxID, :, :), [1 3 4 2]);
                Aw = permute(alpha.*Pr(:, :, :, rxID), [4 5 2 1 3]).*s(:, rxID, :, :, :); % Ns x 1 x Nt x Ntx x Nmcp matrix
                M = obj.network.activeReceivingNodes(rxID).array.numberOfTotalElements;
                Ns = obj.network.activeReceivingNodes(rxID).numberOfSamplesPerCPI;
                if obj.configuration.noise
                    n = 10.^(.05*obj.network.activeReceivingNodes(rxID).noisePowerPerSample_dB).*obj.complexgaussian([Ns, 1, 1, 1, M, obj.numberOfTrialsParallel]);
                else
                    n = 0;
                end
                y{rxID} = zeros(Ns, 1, obj.numberOfTargets, obj.numberOfTransmittingNodes, M, obj.numberOfTrialsParallel);
                % for txID = 1 : obj.numberOfTransmittingNodes
                %     y{rxID}(:, :, :, txID, :, :) = permute(Aw(:, :, :, txID, :), [1 2 3 4 6 5]).*permute(a{txID, rxID}, [6 4 2 5 1 3]) + n;
                % end
                y{rxID} = permute(Aw, [1 2 3 4 6 5]) + n;
            end
        end

        function f = get.fluctuation(obj)
            % (Ntx x Nrx x Nt x Nmcp matrix)
            switch obj.network.networkMode
                case "multiStatic"
                    switch obj.spatialCoherency
                        case "deterministic" % both phase and amplitude is fixed (non random)
                            f = 1;
                        case "coherent"
                            % For each RX-TX pair phases and amplitudes
                            % are known relatively
                            switch obj.swerling
                                case 0 % only absolute phase fluctuates
                                    f = exp(1j*2*pi*rand([1, 1, obj.numberOfTargets, obj.numberOfTrialsParallel]));
                                case 1 % both absolute phase and absolute amplitude fluctuates
                                    f = obj.complexgaussian([1, 1, obj.numberOfTargets, obj.numberOfTrialsParallel]);
                                case 3
                                    error('chi square with 4 DOF not implemented');
                                case {2, 4}
                                    error("multi pulse not implemented");
                            end
                        case "correlatedAmplitudeFluctuations" % Similar to swerling 1 when different receivers are thought as different pulses
                            % For each RX-TX pair phases are independent
                            % uniformly distributed, amplitude are known relatively
                            switch obj.swerling
                                case 0 % whole phase samples fluctuates
                                    f = exp(1j*2*pi*rand([obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes, obj.numberOfTargets, obj.numberOfTrialsParallel]));
                                case 1 % whole phase samples and absolute amplitude fluctuates
                                    f = abs(obj.complexgaussian([1, 1, obj.numberOfTargets, obj.numberOfTrialsParallel])).*exp(1j*2*pi*rand([obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes, obj.numberOfTargets, obj.numberOfTrialsParallel]));
                                case 3
                                    error('chi square with 4 DOF not implemented');
                                case {2, 4}
                                    error("multi pulse not implemented");
                            end
                        case "noncoherent" % Similar to swerling 2 when different receivers are thought as different pulses
                            % For each target, and RX-TX pair
                            % independent complex random variable
                            switch obj.swerling
                                case 0 % whole phase samples fluctuates, MEANINGLESS
                                    f = exp(1j*2*pi*rand([obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes, obj.numberOfTargets, obj.numberOfTrialsParallel]));
                                case 1 % whole complex samples fluctuates
                                    f = obj.complexgaussian([obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes, obj.numberOfTargets, obj.numberOfTrialsParallel]);
                                case 3
                                    error('chi square with 4 DOF not implemented');
                                case {2, 4}
                                    error("multi pulse not implemented");
                            end
                    end
                    if size(f, 2) ~= obj.numberOfReceivingNodes
                        f = repmat(f, [1 obj.numberOfReceivingNodes]);
                    end
                case "monoStatic"
                    f = ones(obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes, obj.numberOfTargets, obj.numberOfTrialsParallel);
                    for rxID = 1 : obj.numberOfReceivingNodes
                        monoStaticTXID = obj.network.monoStaticTransmitterIDs(rxID);
                        nonMonoStaticTXIDs = setdiff(1 : obj.numberOfTransmittingNodes, monoStaticTXID);
                        % Signals from other transmitters will be random
                        f(nonMonoStaticTXIDs, rxID, :, :) = obj.complexgaussian([length(nonMonoStaticTXIDs), 1, obj.numberOfTargets, obj.numberOfTrialsParallel]);
                        switch obj.spatialCoherency
                            case "deterministic"
                                switch obj.swerling
                                    case 0 % both phase and amplitude is fixed (non random)
                                    case 1 % only absolute amplitude fluctuates
                                        f(monoStaticTXID, rxID, :, :) = abs(obj.complexgaussian([1, 1, obj.numberOfTargets, obj.numberOfTrialsParallel]));
                                    case 3
                                        error('chi square with 4 DOF not implemented');
                                    case {2, 4}
                                        error("multi pulse not implemented");
                                end
                            case {"coherent", "correlatedAmplitudeFluctuations", "noncoherent"}
                                switch obj.swerling
                                    case 0 % only absolute phase fluctuates
                                        f(monoStaticTXID, rxID, :, :) = exp(1j*2*pi*rand([1, 1, obj.numberOfTargets, obj.numberOfTrialsParallel]));
                                    case 1 % both absolute phase and absolute amplitude fluctuates
                                        f(monoStaticTXID, rxID, :, :) = obj.complexgaussian([1, 1, obj.numberOfTargets, obj.numberOfTrialsParallel]);
                                    case 3
                                        error('chi square with 4 DOF not implemented');
                                    case {2, 4}
                                        error("multi pulse not implemented");
                                end
                            otherwise
                                error('signal spatial coherency must be coherent for monostatic networks, reconfigure the interface');
                        end
                    end
            end
        end

        function y = get.signalReceivedFromTransmittingNodes(obj)
            % (1 x Nrx cell of Ns x 1 x Ntx x M x Nmcp matrix)
            y = cell(1, obj.numberOfReceivingNodes);
            s = obj.waveformReceivedFromTransmitters;
            for rxID = 1 : obj.numberOfReceivingNodes
                if obj.configuration.pathLoss
                    Pr = permute(10.^(.05*obj.receivedPowerFromTransmittingNodes_dBW(:, rxID, :, :)), [3 2 1 5 4]); % 1 x 1 x Ntx x 1 x Nmcp matrix
                else
                    Pr = 1;
                end
                Aw = Pr.*s(:, rxID, :); % Ns x 1 x Ntx x 1 x Nmcp matrix
                Aw(isinf(Pr)) = 0;
                M = obj.network.activeReceivingNodes(rxID).array.numberOfTotalElements;
                Ns = obj.network.activeReceivingNodes(rxID).numberOfSamplesPerCPI;
                y{rxID} = zeros(Ns, 1, obj.numberOfTransmittingNodes, M, obj.numberOfTrialsParallel);
                for txID = 1 : obj.numberOfTransmittingNodes
                    y{rxID}(:, :, txID, :, :) = Aw(:, :, txID, :, :).*permute(obj.receiveSteeringVectorDirectPath{txID, rxID}, [5 4 2 1 3]);
                end
            end
        end

        function y = get.signalSuperposed(obj)
            % (1 x Nrx cell of Ns x M x Ntxch x Nmcp matrix)
            y = cell(1, obj.numberOfReceivingNodes);
            % transmitted signals perfectly resolved in carrier frequency
            sumDimension = 3; % superpose only signals reflected from targets
            % sumDimension = [3 4]; % superpose all signals transmitted and reflected from targets
            dimensionPermutation = [1 5 4 6 3 2]; % [Ns x 1 x Nt x Ntxch x M x Nmcp] --> [Ns x M x Ntxch x Nmcp x Nt x 1]

            scatteredSignals = obj.signalReceivedFromScatterers;
            if obj.configuration.directPath && strcmp(obj.network.networkMode, "multiStatic")
                signalsDirectPath = obj.signalReceivedFromTransmittingNodes;
            else
                signalsDirectPath = 0;
            end

            for rxID = 1 : obj.numberOfReceivingNodes
                switch obj.network.networkMode
                    case "multiStatic"
                        if obj.configuration.directPath
                            y{rxID} = permute(sum(scatteredSignals{rxID} + permute(signalsDirectPath{rxID}, [1 2 5 3 4]), sumDimension), dimensionPermutation);
                        else
                            y{rxID} = permute(sum(scatteredSignals{rxID}, sumDimension), dimensionPermutation);
                        end
                    case "monoStatic"
                        %%% other TX signals not implemented
                        y{rxID} = permute(sum(scatteredSignals{rxID}(:, :, :, obj.network.monoStaticTransmitterIDs(rxID), :, :), sumDimension), dimensionPermutation);
                end
            end
        end

        function y = get.signalBeamformed(obj)
            % (1 x Nrx cell of Ns x Nrxch x Ntxch x Nmcp matrix)
            y = cell(1, obj.numberOfReceivingNodes);
            s = obj.signalSuperposed;
            for rxID = 1 : obj.numberOfReceivingNodes
                switch obj.network.networkMode
                    case "multiStatic"
                        wavelengths = shiftdim(obj.carrierWavelength, -2);
                    case "monoStatic"
                        wavelengths = obj.carrierWavelength(obj.network.monoStaticTransmitterIDs(rxID));
                end
                % M --> Nrxch (bypass: Nrxch = M)
                y{rxID} = obj.network.activeReceivingNodes(rxID).beamform(s{rxID}, wavelengths);
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

        function R = get.ellipseRotationMatrix(obj)
            % meters (3 x 3 x Ntx x Nrx matrix)
            phi = shiftdim(obj.azimuthBaseline, -2);
            R = [cosd(phi) -sind(phi); sind(phi) cosd(phi)];
        end

        function v = get.ellipse(obj)
            % (3 x Nscan x Ntx x Nrx x Nt x Nmcp matrix)
            t = -180 : .1 : 180;
            v = [shiftdim(obj.majorAxis, -2).*cosd(t)/2; shiftdim(obj.minorAxis, -2).*sind(t)/2];
            v = pagemtimes(obj.ellipseRotationMatrix, v) + permute(obj.centerBaseline(1 : 2, :, :), [1 4 2 3]);
        end

        function v = get.ellipseResolution(obj)
            % (3 x Nscan x Ntx x Nrx x Nt x Nmcp matrix)
            t = -180 : .1 : 180;
            % aForward = obj.bistaticRange + [obj.network.activeReceivingNodes.samplingPeriod]*obj.speedOfLight;
            aForward = obj.bistaticRange + obj.speedOfLight./[obj.network.transmittingNodes.bandWidth]/2;
            bForward = sqrt(aForward.^2 - obj.distanceBaseline.^2);
            vForward = [shiftdim(aForward, -2).*cosd(t)/2; shiftdim(bForward, -2).*sind(t)/2];
            vForward = pagemtimes(obj.ellipseRotationMatrix, vForward) + permute(obj.centerBaseline(1 : 2, :, :), [1 4 2 3]);
            % aBackward = obj.bistaticRange - [obj.network.activeReceivingNodes.samplingPeriod]*obj.speedOfLight;
            aBackward = obj.bistaticRange - obj.speedOfLight./[obj.network.transmittingNodes.bandWidth]/2;
            bBackward = sqrt(aBackward.^2 - obj.distanceBaseline.^2);
            vBackward = [shiftdim(aBackward, -2).*cosd(t)/2; shiftdim(bBackward, -2).*sind(t)/2];
            vBackward = pagemtimes(obj.ellipseRotationMatrix, vBackward) + permute(obj.centerBaseline(1 : 2, :, :), [1 4 2 3]);
            v = cat(7, real(vBackward), vForward);
        end

        %%% unit vectors

        function u = get.unitDirectionRX(obj) % unit vector from receiving nodes to target
            % (3 x Nrx x Nt x Nmcp matrix)
            if isempty(obj.targets) || ~logical(obj.numberOfReceivingNodes)
                if logical(obj.numberOfReceivingNodes)
                    u = zeros(3, obj.numberOfReceivingNodes);
                elseif ~isempty(obj.targets)
                    u = zeros(3, 1, obj.numberOfTargets);
                else
                    u = zeros(3, 1);
                end
            else
                u = permute(obj.positions, [1 4 2 3]) - [obj.network.activeReceivingNodes.position];
                u = u./sqrt(sum(abs(u).^2));
                u(isnan(u)) = 0;
            end
        end

        function u = get.unitDirectionTX(obj) % unit vector from transmitting nodes to target
            % (3 x Ntx x Nt x Nmcp matrix)
            if isempty(obj.targets) || ~logical(obj.numberOfTransmittingNodes)
                if logical(obj.numberOfTransmittingNodes)
                    u = zeros(3, obj.numberOfTransmittingNodes);
                elseif ~isempty(obj.targets)
                    u = zeros(3, 1, obj.numberOfTargets);
                else
                    u = zeros(3, 1);
                end
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
            if ~logical(obj.numberOfReceivingNodes) || ~logical(obj.numberOfTransmittingNodes)
                if logical(obj.numberOfReceivingNodes)
                    t = zeros(3, 1, obj.numberOfReceivingNodes);
                elseif logical(obj.numberOfTransmittingNodes)
                    t = zeros(3, obj.numberOfTransmittingNodes);
                end
            else
                t = obj.network.positionReceivingNode - obj.network.positionTransmittingNode; % baseline vector
                t = t./sqrt(sum(abs(t).^2, 1));
                t(isnan(t)) = 0;
            end
        end

        function n = get.unitNormalBaseline(obj)
            % (3 x Ntx x Nrx x Nt matrix)
            u = repmat(permute(obj.unitDirectionRX, [1 5 2 3 4]), 1, obj.numberOfTransmittingNodes); % rx-target vector
            t = repmat(obj.unitDirectionBaseline, [1 1 1 obj.numberOfTargets, size(u, 5)]); % baseline vector
            n = u - dot(t, u).*t;
            n = n./sqrt(sum(abs(n).^2, 1));
            n(isnan(n)) = 0;
        end

        %%% distances

        function R = get.distanceRX(obj)
            % (1 x Nrx x Nt x Nmcp matrix)
            if isempty(obj.targets) || ~logical(obj.numberOfReceivingNodes)
                if logical(obj.numberOfReceivingNodes)
                    R = inf(1, obj.numberOfReceivingNodes);
                elseif ~isempty(obj.targets)
                    R = inf(1, 1, obj.numberOfTargets);
                else
                    R = inf;
                end
            else
                R = sqrt(sum((permute(obj.positions, [1 4 2 3]) - [obj.network.activeReceivingNodes.position]).^2));
            end
        end

        function R = get.distanceTX(obj)
            % (Ntx x 1 x Nt x Nmcp matrix)
            if isempty(obj.targets) || ~logical(obj.numberOfTransmittingNodes)
                if logical(obj.numberOfTransmittingNodes)
                    R = inf(obj.numberOfTransmittingNodes, 1);
                elseif ~isempty(obj.targets)
                    R = inf(1, 1, obj.numberOfTargets);
                else
                    R = inf;
                end
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
            % (Ntx x Nrx matrix)
            d = obj.network.distanceBaseline;
        end

        function d = get.centerBaseline(obj)
            % (3 x Ntx x Nrx matrix)
            d = obj.network.centerBaseline;
        end

        function phi = get.azimuthBaseline(obj)
            % degrees (Ntx x Nrx matrix)
            phi = shiftdim(atan2d(obj.unitDirectionBaseline(2, :, :), obj.unitDirectionBaseline(1, :, :)), 1);
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
                obj.positions = obj.targets.position + options.width.*(rand(3, obj.numberOfTargets, obj.numberOfTrialsParallel) - 0.5);
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
                options.fractionalDelay (1, 1) string {mustBeMember(options.fractionalDelay, ["sinc-based", "lagrange-based", "off"])} = obj.configuration.fractionalDelay
                options.timeSynchronization (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.timeSynchronization, [0, 1])} = obj.configuration.timeSynchronization
                options.frequencySynchronization (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.frequencySynchronization, [0, 1])} = obj.configuration.frequencySynchronization
                options.spatialCoherency (1, 1) string {mustBeMember(options.spatialCoherency, ["deterministic", "coherent", "correlatedAmplitudeFluctuations", "noncoherent"])} = obj.spatialCoherency
                options.swerling (1, 1) double {mustBeInteger, mustBeMember(options.swerling, 0 : 4)} = obj.swerling
            end
            obj.configuration.noise = options.noise;
            obj.configuration.directPath = options.directPath;
            obj.configuration.pathLoss = options.pathLoss;
            obj.configuration.fractionalDelay = options.fractionalDelay;
            obj.configuration.timeSynchronization = options.timeSynchronization;
            obj.configuration.frequencySynchronization = options.frequencySynchronization;
            obj.spatialCoherency = options.spatialCoherency;
            obj.swerling = options.swerling;
        end

        %%% visualization methods

        function fig = visualizescenario(obj, options)
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
                fig = figure;
            else
                fig = figure(options.figureID);
            end
            posRX = [obj.network.activeReceivingNodes.position]/1e3;
            posTX = [obj.network.activeTransmittingNodes.position]/1e3;
            plot3(posRX(1, :), posRX(2, :), posRX(3, :), 'xb', 'LineWidth', 2, 'MarkerSize', 10);
            hold on;
            plot3(posTX(1, :), posTX(2, :), posTX(3, :), '+r', 'LineWidth', 2, 'MarkerSize', 10);
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
            % text(posTX(1, :), posTX(2, :), posTX(3, :), num2str(find(obj.network.transmittingNodeActivity).'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            % if obj.numberOfTargets < 11 && ~isempty(obj.targets)
            %     text(x, y, z, num2str((1 : obj.numberOfTargets).'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            % end
            if options.showPattern
                if obj.numberOfTargets > 20
                    for txID = 1 : obj.numberOfTransmittingNodes
                        visualizationNormalization = max(obj.distanceTX(txID, :, :), [], 'all')/1e3;
                        G = visualizationNormalization + permute(abs(obj.transmittedBeam(txID, :, :, 1).^2./obj.network.activeTransmittingNodes(txID).array.numberOfTotalElements), [1 3 2]);
                        G(G < 0) = nan;
                        u = permute(obj.unitDirectionTX, [1 3 2 4]);
                        xu = u(1, :, txID, 1);
                        yu = u(2, :, txID, 1);
                        zu = u(3, :, txID, 1);
                        plot3(G.*xu + posTX(1, txID), G.*yu + posTX(2, txID), G.*zu + posTX(3, txID), '.r');
                        for rxID = 1 : obj.numberOfReceivingNodes
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
                    np = 1441;
                    phi = linspace(-180, 180, np);
                    d = max(sqrt(sum(obj.network.boundaryListened.^2)))/10;
                    for txID = 1 : obj.numberOfTransmittingNodes
                        targetsTemp = target( ...
                            'position', [ ...
                            d*cosd(phi) + posTX(1, txID)*1e3; ...
                            d*sind(phi) + posTX(2, txID)*1e3; ...
                            zeros(1, np) + posTX(3, txID)*1e3]);
                        intTemp = interface( ...
                            'network', obj.network, ...
                            'targets', targetsTemp);
                        visualizationNormalization = max(intTemp.distanceTX(txID, :, :, :), [], 'all')/1e3;
                        G = visualizationNormalization.*permute(abs(intTemp.transmittedBeam(txID, :, :).^2./intTemp.network.activeTransmittingNodes(txID).array.numberOfTotalElements), [1 3 2]);
                        u = permute(intTemp.unitDirectionTX, [1 3 2 4]);
                        plot3(G.*u(1, :, txID, 1) + posTX(1, txID), G.*u(2, :, txID, 1) + posTX(2, txID), u(3, :, txID, 1) + posTX(3, txID), 'r');
                        for rxID = 1 : obj.numberOfReceivingNodes
                            targetsTemp = target( ...
                                'position', [ ...
                                d*cosd(phi) + posRX(1, rxID)*1e3; ...
                                d*sind(phi) + posRX(2, rxID)*1e3; ...
                                zeros(1, np) + posRX(3, rxID)*1e3]);
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
                for rxID = 1 : obj.numberOfReceivingNodes
                    n = obj.network.activeReceivingNodes(rxID).array.normalVector;
                    quiver3(posRX(1, rxID), posRX(2, rxID), posRX(3, rxID), n(1), n(2), n(3), 'Color', 'm', 'LineWidth', 2);
                end
                for txID = 1 : obj.numberOfTransmittingNodes
                    n = obj.network.activeTransmittingNodes(txID).array.normalVector;
                    quiver3(posTX(1, txID), posTX(2, txID), posTX(3, txID), n(1), n(2), n(3), 'Color', 'm', 'LineWidth', 2);
                end
            end
            if options.showSteeringDirection
                for rxID = 1 : obj.numberOfReceivingNodes
                    s = obj.network.activeReceivingNodes(rxID).array.steeringUnitDirection;
                    quiver3(posRX(1, rxID), posRX(2, rxID), posRX(3, rxID), s(1), s(2), s(3), 'Color', 'm', 'LineWidth', 2, 'LineStyle', '--');
                end
                for txID = 1 : obj.numberOfTransmittingNodes
                    s = obj.network.activeTransmittingNodes(txID).array.steeringUnitDirection;
                    quiver3(posTX(1, txID), posTX(2, txID), posTX(3, txID), s(1), s(2), s(3), 'Color', 'm', 'LineWidth', 2, 'LineStyle', '--');
                end
            end
            posRX = repelem(posRX, 1, obj.numberOfTransmittingNodes);
            posTX = repmat(posTX, 1, obj.numberOfReceivingNodes);
            line([posRX(1, :); posTX(1, :)], [posRX(2, :); posTX(2, :)], [posRX(3, :); posTX(3, :)], 'lineStyle', '--', 'Color', 'k');
            xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
            if ~isempty(obj.targets)
                legend('RX', 'TX', 'targets', 'Location', 'best');
            else
                legend('RX', 'TX', 'Location', 'best');
            end
            grid off; grid on; grid minor; axis equal; % title('scenario');
            view(0, 90); hold off; drawnow;
        end

        function fig = visualizeellipses(obj, options)
            arguments
                obj
                options.ellipseType (1, 1) string {mustBeMember(options.ellipseType, ["target", "resolution"])} = "target"
                options.targetID (1, 1) double {mustBePositive, mustBeInteger} = 1;
                options.trialID (1, 1) double {mustBePositive, mustBeInteger} = 1;
            end
            switch options.ellipseType
                case "target"
                    x = obj.ellipse(1, :, :, :, options.targetID, options.trialID)/1e3; % (3 x Nscan x Ntx x Nrx x Nt x Nmcp matrix)
                    y = obj.ellipse(2, :, :, :, options.targetID, options.trialID)/1e3; % (3 x Nscan x Ntx x Nrx x Nt x Nmcp matrix)
                case "resolution"
                    x = obj.ellipseResolution(1, :, :, :, options.targetID, options.trialID, :)/1e3; % (3 x Nscan x Ntx x Nrx x Nt x Nmcp x 2 matrix)
                    y = obj.ellipseResolution(2, :, :, :, options.targetID, options.trialID, :)/1e3; % (3 x Nscan x Ntx x Nrx x Nt x Nmcp x 2 matrix)
            end
            xTarget = obj.positions(1, options.targetID, options.trialID)/1e3;
            yTarget = obj.positions(2, options.targetID, options.trialID)/1e3;
            fig = figure;
            hold on;
            grid off; grid on; grid minor;
            for rxID = 1 : obj.numberOfReceivingNodes
                for txID = 1 : obj.numberOfTransmittingNodes
                    switch options.ellipseType
                        case "target"
                            plot(x(:, :, txID, rxID, :, :), y(:, :, txID, rxID, :, :));
                        case "resolution"
                            xPlot = [x(:, :, txID, rxID, :, :, 1), fliplr(x(:, :, txID, rxID, :, :, 2))];
                            yPlot = [y(:, :, txID, rxID, :, :, 1), fliplr(y(:, :, txID, rxID, :, :, 2))];
                            % patch(xPlot, yPlot, 'r', 'FaceAlpha', 0.2);
                            plot(x(:, :, txID, rxID, :, 1), y(:, :, txID, rxID, :, 1), 'k', 'LineWidth', 2); hold on;
                            plot(x(:, :, txID, rxID, :, 2), y(:, :, txID, rxID, :, 2), 'k', 'LineWidth', 2);
                            % patch(xPlot, yPlot, 'y');
                    end
                end
            end
            switch options.ellipseType
                case "target"
                    plot(xTarget, yTarget, 'ok', 'LineWidth', 2, 'MarkerSize', 20);
                    % title('iso bistatic range ellipse of the target');
                case "resolution"
                    % plot(xTarget, yTarget, 'k+', 'LineWidth', 2, 'MarkerSize', 10);
                    % title('resolution ellipses of the target');
            end
            xlabel('x-pos (km)'); ylabel('y-pos (km)'); zlabel('z-pos (km)');
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
                t = obj.network.activeReceivingNodes(rxID).samplingInstants*1e6;
                if isscalar(options.receivingNodeIDs)
                    y = 20*log10(abs(s{rxID}(:, :, options.trialID)));
                else
                    y = 20*log10(abs(s{rxID}(:, ceil(end/2), options.trialID)));
                end
                y = y - obj.network.activeReceivingNodes(rxID).noisePowerPerSample_dB;
                plot(t, y, 'LineWidth', 2); hold on;
            end
            grid off; grid on; grid minor;
            xlabel('time delay (\mus)'); ylabel('Per Sample SNR (dB)');
            xline(20, 'LineStyle', '--'); 
            xline(29.5, 'LineStyle', '--');
            % title('received signal');
            % if ~isscalar(options.receivingNodeIDs)
            %     leg = legend(num2str(options.receivingNodeIDs.'), 'Location', 'best');
            %     title(leg, 'RX ID');
            % end
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
                t = obj.network.activeReceivingNodes(rxID).samplingInstants*1e6;
                if isscalar(options.receivingNodeIDs)
                    y = 20*log10(abs(s{rxID}(:, :, options.trialID)));
                else
                    y = 20*log10(abs(s{rxID}(:, ceil(end/2), options.trialID)));
                end
                y = y - obj.network.activeReceivingNodes(rxID).noisePowerPerSample_dB;
                plot(t, y); hold on;
            end
            grid off; grid on; grid minor;
            xlabel('time delay (\mus)'); ylabel('Per Sample SNR (dB)');
            xline(120, 'LineStyle', '--'); 
            xline(129.5, 'LineStyle', '--');
            % title('received and beamformed signal');
            % if ~isscalar(options.receivingNodeIDs)
            %     leg = legend(num2str(options.receivingNodeIDs.'), 'Location', 'best');
            %     title(leg, 'RX ID');
            % end
            hold off; drawnow;
        end

        function visualizeisocurves(obj, options)
            arguments
                obj
                options.curves (1, :) string {mustBeMember(options.curves, ["delay", "power", "doppler", "bistaticAngle", "triangulation", "lossFactor", "all"])} = "delay"
                options.numberOfCurves (1, 1) double {mustBeInteger, mustBePositive} = 10
                options.transmittingNodeID (1, 1) double {mustBeInteger, mustBePositive} = 1
                options.receivingNodeID (1, 1) double {mustBeInteger, mustBePositive} = 1
                options.saveFigure (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.saveFigure, [0, 1])} = false
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
                fig = figure; contourf(x, y, reshape(-obj.timeDelay(options.transmittingNodeID, options.receivingNodeID, :, :)*1e6, L1, L2), options.numberOfCurves, "ShowText", true, "LabelFormat", "%3.3g us");
                hold on; plot(posRX(axID1, :), posRX(axID2, :), '.b', 'LineWidth', 2, 'MarkerSize', 5);
                plot(posTX(axID1, :), posTX(axID2, :), '.r', 'LineWidth', 2, 'MarkerSize', 5);
                text(posRX(axID1, :), posRX(axID2, :), 'RX', 'Color', 'b', 'VerticalAlignment', 'bottom');
                text(posTX(axID1, :), posTX(axID2, :), 0, 'TX', 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
                xlabel(firstAxisLabel); ylabel(secondAxisLabel);
                xlim tight; ylim tight; zlim tight; view(2);
                if options.saveFigure
                    figureName = 'fund_iso_delay';
                    savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
                    saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
                else
                    title('iso delay curves'); drawnow;
                end
            end
            if any(strcmpi(options.curves, "power")) || any(strcmpi(options.curves, "all"))
                fig = figure; contourf(x, y, obj.network.beamformingGain_dB + reshape(obj.averageSNR_dB(options.transmittingNodeID, options.receivingNodeID, :), L1, L2), options.numberOfCurves, "ShowText", true, "LabelFormat", "%3.3g dB");
                % hold on; plot(posRX(axID1, :), posRX(axID2, :), 'ob', 'LineWidth', 2, 'MarkerSize', 10);
                % plot(posTX(axID1, :), posTX(axID2, :), 'or', 'LineWidth', 2, 'MarkerSize', 10);
                xlabel(firstAxisLabel); ylabel(secondAxisLabel);
                xlim tight; ylim tight; zlim tight; view(2);
                if options.saveFigure
                    figureName = 'fund_iso_power';
                    savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
                    saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
                else
                    title('iso power curves'); drawnow;
                end
            end
            if any(strcmpi(options.curves, "doppler")) || any(strcmpi(options.curves, "all"))
                fig = figure; contourf(x, y, reshape(obj.dopplerShift(options.transmittingNodeID, options.receivingNodeID, :), L1, L2), options.numberOfCurves, "ShowText", true, "LabelFormat", "%3.3g Hz");
                % hold on; plot(posRX(axID1, :), posRX(axID2, :), 'ob', 'LineWidth', 2, 'MarkerSize', 10);
                % plot(posTX(axID1, :), posTX(axID2, :), 'or', 'LineWidth', 2, 'MarkerSize', 10);
                xlabel(firstAxisLabel); ylabel(secondAxisLabel);
                xlim tight; ylim tight; zlim tight; view(2);
                if options.saveFigure
                    figureName = 'fund_iso_doppler';
                    savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
                    saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
                else
                    title('iso doppler curves'); drawnow;
                end
            end
            if any(strcmpi(options.curves, "bistaticAngle")) || any(strcmpi(options.curves, "all"))
                fig = figure; contourf(x, y, reshape(obj.bistaticAngle(options.transmittingNodeID, options.receivingNodeID, :), L1, L2), linspace(0, 180, options.numberOfCurves + 3), "ShowText", true, "LabelFormat", "%3.3g°");
                % hold on; plot(posRX(axID1, :), posRX(axID2, :), 'ob', 'LineWidth', 2, 'MarkerSize', 10);
                % plot(posTX(axID1, :), posTX(axID2, :), 'or', 'LineWidth', 2, 'MarkerSize', 10);
                xlabel(firstAxisLabel); ylabel(secondAxisLabel);
                xlim tight; ylim tight; zlim tight; view(2);
                if options.saveFigure
                    figureName = 'fund_iso_bistatic_angle';
                    savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
                    saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
                else
                    title('iso bistatic angle curves'); drawnow;
                end
            end
            if any(strcmpi(options.curves, "triangulation")) || any(strcmpi(options.curves, "all"))
                fig = figure; contourf(x, y, reshape(obj.triangulationFactor(options.transmittingNodeID, options.receivingNodeID, :), L1, L2), linspace(0, 1, options.numberOfCurves + 1), "ShowText", true, "LabelFormat", "%3.3g");
                hold on; plot(posRX(axID1, :), posRX(axID2, :), 'ob', 'LineWidth', 2, 'MarkerSize', 10);
                plot(posTX(axID1, :), posTX(axID2, :), 'or', 'LineWidth', 2, 'MarkerSize', 10);
                xlabel(firstAxisLabel); ylabel(secondAxisLabel);
                xlim tight; ylim tight; zlim tight; view(2);
                if options.saveFigure
                    figureName = 'fund_iso_triangulation';
                    savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
                    saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
                else
                    title('iso triangulation curves'); drawnow;
                end
            end
            if any(strcmpi(options.curves, "lossFactor")) || any(strcmpi(options.curves, "all"))
                fig = figure; contourf(x, y, reshape(obj.lossFactor(options.transmittingNodeID, options.receivingNodeID, :), L1, L2), linspace(0, 1, options.numberOfCurves + 1), "ShowText", true, "LabelFormat", "%3.3g");
                hold on; plot(posRX(axID1, :), posRX(axID2, :), 'ob', 'LineWidth', 2, 'MarkerSize', 10);
                plot(posTX(axID1, :), posTX(axID2, :), 'or', 'LineWidth', 2, 'MarkerSize', 10);
                xlabel(firstAxisLabel); ylabel(secondAxisLabel);
                xlim tight; ylim tight; zlim tight; view(2);
                if options.saveFigure
                    figureName = 'fund_iso_lossFactor';
                    savefig(fig, ['C:\GitRepo\MSRS\figures\' figureName '.fig']);
                    saveas(fig, ['C:\GitRepo\MSRS\figures\' figureName '.eps'], 'epsc');
                else
                    title('iso lossFactor curves'); drawnow;
                end
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