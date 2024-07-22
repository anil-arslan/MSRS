classdef interface < handle
    %interface Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = private, GetAccess = public)
        network (1, :) radarNetwork = radarNetwork.empty()
        targets (1, 1) target = target;
        environment (1, 1) %%%%% not implemented
        configuration (1, 1) struct = struct('noise', false, 'directPath', false)
    end

    properties (Dependent)
        %%% target parameters
        bistaticRange double % meters (Ntx x Nrx x Nt matrix) major axis length
        timeDelay double % seconds (Ntx x Nrx x Nt matrix)
        dopplerShift double % Hz (Ntx x Nrx x Nt matrix)
        dopplerShiftBisector double % Hz (Ntx x Nrx x Nt matrix)
        receivedPowerFromScatterers_dBW double % dB watt (Ntx x Nrx x Nt matrix)
        receivedPowerFromTransmittingNodes_dBW double % dB watt (Ntx x Nrx matrix)
        inputSNR_dB double % dB (Ntx x Nrx x Nt matrix)
            %%% signal
            signalReceivedFromScatterers double % complex received signal (Ns x Nrx x Nt x Ntx matrix)
            signalReceivedFromTransmittingNodes double % complex received direct path signal (Ns x Nrx x Ntx matrix)
            signalSuperposed double % (Ns x Nrx matrix)

        %%% ellipse
        majorAxis double % meters (Ntx x Nrx x Nt matrix)
        minorAxis double % meters (Ntx x Nrx x Nt matrix)
        eccentricity double % meters (Ntx x Nrx x Nt matrix)
        latusRectum double % meters (Ntx x Nrx x Nt matrix)

        %%% unit vectors
        unitDirectionRX (3, :, :) double % (3 x Nrx x Nt matrix)
        unitDirectionTX (3, :, :) double % (3 x Ntx x Nt matrix)
        unitDirectionBisector (3, :, :, :) double % (3 x Ntx x Nrx x Nt matrix)
        unitDirectionBaseline (3, :, :, :) double % (3 x Ntx x Nrx x Nt matrix)
        unitNormalBaseline (3, :, :, :) double % (3 x Ntx x Nrx x Nt matrix)

        %%% distances
        distanceRX (1, :, :) double % meters (1 x Nrx x Nt matrix)
        distanceTX (:, 1, :) double % meters (Ntx x 1 x Nt matrix)
        distanceBisector double % (Ntx x Nrx x Nt matrix)
        distanceBaseline double % (Ntx x Nrx matrix)

        %%% angles
        bistaticAngle double % (Ntx x Nrx x Nt matrix)
        elevationTargetRX double % (Nrx x Nt matrix)
        azimuthTargetRX double % (Nrx x Nt matrix)
        elevationTargetTX double % (Ntx x Nt matrix)
        azimuthTargetTX double % (Ntx x Nt matrix)
        angleVelocityBisector double % (Ntx x Nrx x Nt matrix)
        angleBaselineRX double % (Ntx x Nrx x Nt matrix)
        angleBaselineTX double % (Ntx x Nrx x Nt matrix)

        %%% factors
        triangulationFactor double {mustBeNonnegative} % (Ntx x Nrx x Nt matrix)
        lossFactor double {mustBeNonnegative} % (Ntx x Nrx x Nt matrix)
        targetRegion string % (Ntx x Nrx x Nt matrix)
        
        %%% sizes
        numberOfTargets (1, 1) double {mustBeInteger, mustBeNonnegative} % Nt
        numberOfReceivingNodes (1, 1) double {mustBeInteger, mustBeNonnegative} % Nrx
        numberOfTransmittingNodes (1, 1) double {mustBeInteger, mustBeNonnegative} % Ntx
        numberOfBistaticPairs (1, 1) double {mustBeInteger, mustBeNonnegative} % Nrx*Ntx
    end

    properties (Constant)
        speedOfLight (1, 1) double {mustBePositive} = physconst('LightSpeed');
    end

    methods
        function obj = interface(options)
            arguments
                options.network (1, 1) radarNetwork = radarNetwork.empty()
                options.targets (1, 1) target = target;
            end
            %interface Construct an instance of this class
            %   Detailed explanation goes here
            obj.network = options.network;
            obj.targets = options.targets;
        end

        %%% target parameters

        function rho = get.bistaticRange(obj) % major axis
            % (Ntx x Nrx x Nt matrix)
            rho = reshape(obj.distanceRX + obj.distanceTX, obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes, obj.numberOfTargets);
        end

        function tau = get.timeDelay(obj)
            % (Ntx x Nrx x Nt matrix)
            tau = obj.bistaticRange/obj.speedOfLight;
        end

        function Pr = get.receivedPowerFromScatterers_dBW(obj)
            % (Ntx x Nrx x Nt matrix)
            powerGains = 10*log10([obj.network.activeTransmittingNodes.inputPower_W].') + permute(obj.targets.RCS_dbms, [1 3 2]) ...
            + 20*log10([obj.network.activeTransmittingNodes.carrierWavelength].');
            powerLosses = [obj.network.activeReceivingNodes.systemLoss_dB] + 30*log10(4*pi) + 20*log10(obj.distanceRX.*obj.distanceTX);
            Pr = powerGains - powerLosses;
        end

        function Pr = get.receivedPowerFromTransmittingNodes_dBW(obj)
            % (Ntx x Nrx matrix)
            %%%%% not implemented %%%%%
            powerGains = 10*log10([obj.network.activeTransmittingNodes.inputPower_W].') ...
            + 20*log10([obj.network.activeTransmittingNodes.carrierWavelength].');
            powerLosses = [obj.network.activeReceivingNodes.systemLoss_dB] + 30*log10(4*pi) + 20*log10(obj.distanceBaseline);
            Pr = powerGains - powerLosses;
        end

        function SNRin = get.inputSNR_dB(obj)
            % (Ntx x Nrx x Nt matrix)
            SNRin = obj.receivedPowerFromScatterers_dBW - [obj.network.activeReceivingNodes.noisePowerPerSample_dB];
        end

        %%% doppler

        function fd = get.dopplerShift(obj)
            % (Ntx x Nrx x Nt matrix)
            vdRx = pagemtimes(permute(obj.targets.velocity, [3 1 4 2]), permute(obj.unitDirectionRX, [1 4 2 3]));
            vdTx = pagemtimes(permute(obj.targets.velocity, [3 1 4 2]), permute(obj.unitDirectionTX, [1 2 4 3]));
            fd = shiftdim((vdRx + vdTx), 1)./[obj.network.activeTransmittingNodes.carrierWavelength].';
        end

        function fd = get.dopplerShiftBisector(obj)
            % (Ntx x Nrx x Nt matrix)
            fd = 2*shiftdim(obj.targets.velocityMagnitude, -1).*cosd(obj.bistaticAngle/2).*cosd(obj.angleVelocityBisector)./[obj.network.activeTransmittingNodes.carrierWavelength].';
        end

        %%% signal

        function s = get.signalReceivedFromScatterers(obj)
            % (Ns x Nrx x Nt x Ntx matrix)
            tau = [obj.network.activeReceivingNodes.samplingInstants] - permute(obj.timeDelay, [4, 2, 3, 1]); % Ns x Nrx x Nt x Ntx matrix
            fd = permute(obj.dopplerShift, [4, 2, 3, 1]); % 1 x Nrx x Nt x Ntx matrix
            fc = shiftdim([obj.network.activeTransmittingNodes.carrierFrequency], -2); % 1 x 1 x 1 x Ntx matrix
            A = permute(10.^(.05*obj.receivedPowerFromScatterers_dBW), [4, 2, 3, 1]); % 1 x Nrx x Nt x Ntx matrix
            waveforms = zeros(size(tau)); % Ns x Nrx x Nt x Ntx matrix
            for txID = 1 : obj.numberOfTransmittingNodes
                waveforms(:, :, :, txID) = obj.network.activeTransmittingNodes(txID).waveform(-tau(:, :, :, txID));
            end
            s = A.*waveforms.*exp(-1j*2*pi*(fc + fd).*tau); % Ns x Nrx x Nt x Ntx matrix
            if obj.configuration.noise
                Ns = [obj.network.activeReceivingNodes.numberOfSamplesPerCPI];
                NsUnique = unique(Ns);
                switch length(NsUnique)
                    case 1
                        s = s + 10.^(.05*[obj.network.activeReceivingNodes.noisePowerPerSample_dB])/sqrt(2).*(randn(NsUnique, obj.numberOfReceivingNodes) + 1j*randn(NsUnique, obj.numberOfReceivingNodes));
                    otherwise
                        error('not implemented');
                end
            end
        end

        function s = get.signalReceivedFromTransmittingNodes(obj)
            % (Ns x Nrx x Ntx matrix)
            %%%%% not implemented %%%%%
            tau = [obj.network.activeReceivingNodes.samplingInstants] - permute(obj.network.directPathDelay, [3, 2, 1]); % Ns x Nrx x Ntx matrix
            fc = shiftdim([obj.network.activeTransmittingNodes.carrierFrequency], -1); % 1 x 1 x Ntx matrix
            A = permute(10.^(.05*obj.receivedPowerFromTransmittingNodes_dBW), [3, 2, 1]); % 1 x Nrx x Ntx matrix
            waveforms = zeros(size(tau)); % Ns x Nrx x Ntx matrix
            for txID = 1 : obj.network.numberOfActiveTransmittingNodes
                waveforms(:, :, txID) = obj.network.activeTransmittingNodes(txID).waveform(tau(:, :, txID));
            end
            s = A.*waveforms.*exp(1j*2*pi*fc.*tau); % Ns x Nrx x Ntx matrix
        end

        function s = get.signalSuperposed(obj)
            % (Ns x Nrx matrix)
            if obj.configuration.directPath
                % s = sum(obj.signalReceivedFromScatterers + permute(obj.signalReceivedFromTransmittingNodes, [1 2 4 3]), [3, 4]);
                error('not implemented');
            else
                s = sum(obj.signalReceivedFromScatterers, [3, 4]);
            end
        end

        %%% ellipse

        function a = get.majorAxis(obj)
            % (Ntx x Nrx x Nt matrix)
            a = obj.bistaticRange;
        end

        function b = get.minorAxis(obj)
            % (Ntx x Nrx x Nt matrix)
            b = sqrt(obj.majorAxis.^2 - obj.distanceBaseline.^2);
        end

        function e = get.eccentricity(obj)
            % (Ntx x Nrx x Nt matrix)
            e = obj.distanceBaseline./obj.bistaticRange;
        end

        function l = get.latusRectum(obj)
            % (Ntx x Nrx x Nt matrix)
            l = obj.majorAxis.*(1 - obj.eccentricity.^2);
        end

        %%% unit vectors

        function u = get.unitDirectionRX(obj) % unit vector from radar to target
            % (3 x Nrx x Nt matrix)
            u = permute(obj.targets.position, [1, 3, 2]) - [obj.network.activeReceivingNodes.position];
            u = u./sqrt(sum(abs(u).^2));
            u(isnan(u)) = 0;
        end

        function u = get.unitDirectionTX(obj) % unit vector from radar to target
            % (3 x Ntx x Nt matrix)
            u = permute(obj.targets.position, [1, 3, 2]) - [obj.network.activeTransmittingNodes.position];
            u = u./sqrt(sum(abs(u).^2));
            u(isnan(u)) = 0;
        end

        function g = get.unitDirectionBisector(obj)
            % (3 x Ntx x Nrx x Nt matrix)
            g = permute(obj.unitDirectionRX, [1 4 2 3]) + permute(obj.unitDirectionTX, [1 2 4 3]);
            g = g./sqrt(sum(abs(g).^2, 1));
            g(isnan(g)) = 0;
        end

        function t = get.unitDirectionBaseline(obj)
            % (3 x Ntx x Nrx x Nt matrix)
            t = repmat(obj.network.positionReceivingNode - obj.network.positionTransmittingNode, [1 1 1 obj.numberOfTargets]); % baseline vector
            t = t./sqrt(sum(abs(t).^2, 1));
            t(isnan(t)) = 0;
        end

        function n = get.unitNormalBaseline(obj)
            % (3 x Ntx x Nrx x Nt matrix)
            t = obj.unitDirectionBaseline; % baseline vector
            u = repmat(permute(obj.unitDirectionRX, [1, 4, 2, 3]), 1, obj.numberOfTransmittingNodes); % rx-target vector
            n = u - dot(t, u).*t;
            n = n./sqrt(sum(abs(n).^2, 1));
            n(isnan(n)) = 0;
        end

        %%% distances

        function R = get.distanceRX(obj)
            % (1 x Nrx x Nt matrix)
            R = sqrt(sum((permute(obj.targets.position, [1, 3, 2]) - [obj.network.activeReceivingNodes.position]).^2));
        end

        function R = get.distanceTX(obj)
            % (Ntx x 1 x Nt matrix)
            R = permute(sqrt(sum((permute(obj.targets.position, [1, 3, 2]) - [obj.network.activeTransmittingNodes.position]).^2)), [2, 1, 3]);
        end

        function d = get.distanceBisector(obj)
            % (Ntx x Nrx x Nt matrix)
            m = repmat(obj.network.positionReceivingNode + obj.network.positionTransmittingNode, [1, 1, 1, obj.numberOfTargets])/2;
            d = permute(sqrt(sum(abs(m - permute(obj.targets.position, [1 3 4 2])).^2)), [2 3 4 1]);
        end

        function d = get.distanceBaseline(obj)
            % (Ntx x Nrx x Nt matrix)
            d = obj.network.distanceBaseline;
        end

        %%% angles

        function ksi = get.bistaticAngle(obj)
            % (Ntx x Nrx x Nt matrix)
            ksi = permute(real(acosd(pagemtimes(permute(obj.unitDirectionTX, [2 1 4 3]), permute(obj.unitDirectionRX, [1 4 2 3])))), [1 3 4 2]);
            %ksi = obj.localAngleBaselineRX + obj.localAngleBaselineTX;
            %ksi - 90 needed sometimes
        end

        function theta = get.elevationTargetRX(obj)
            % (Nrx x Nt matrix)
            theta = asind(shiftdim(obj.unitDirectionRX(3, :, :), 1));
        end

        function phi = get.azimuthTargetRX(obj)
            % (Nrx x Nt matrix)
            phiSign = sign(asind(shiftdim(obj.unitDirectionRX(2, :, :), 1)./cosd(obj.elevationTargetRX)));
            phiSign(isnan(phiSign)) = 1;
            phiSign(~phiSign) = 1;
            phi = phiSign.*real(acosd(shiftdim(obj.unitDirectionRX(1, :, :), 1)./cosd(obj.elevationTargetRX)));
        end

        function theta = get.elevationTargetTX(obj)
            % (Ntx x Nt matrix)
            theta = asind(shiftdim(obj.unitDirectionTX(3, :, :), 1));
        end

        function phi = get.azimuthTargetTX(obj)
            % (Ntx x Nt matrix)
            phiSign = sign(asind(shiftdim(obj.unitDirectionTX(2, :, :), 1)./cosd(obj.elevationTargetTX)));
            phiSign(isnan(phiSign)) = 1;
            phiSign(~phiSign) = 1;
            phi = phiSign.*real(acosd(shiftdim(obj.unitDirectionTX(1, :, :), 1)./cosd(obj.elevationTargetTX)));
        end

        function ksi = get.angleVelocityBisector(obj)
            % (Ntx x Nrx x Nt matrix)
            ksi = permute(real(acosd(pagemtimes(permute(obj.unitDirectionBisector, [2 1 3 4]), permute(obj.targets.velocityUnitDirection, [1 3 4 2])))), [1 3 4 2]);
        end

        function ksi = get.angleBaselineRX(obj)
            % (Ntx x Nrx x Nt matrix)
            ksi = permute(real(acosd(pagemtimes(permute(obj.unitNormalBaseline, [2 1 4 3]), permute(obj.unitDirectionRX, [1 4 3 2])))), [1 4 3 2]);
        end

        function ksi = get.angleBaselineTX(obj)
            % (Ntx x Nrx x Nt matrix)
            ksi = permute(real(acosd(pagemtimes(permute(obj.unitNormalBaseline, [3 1 4 2]), permute(obj.unitDirectionTX, [1 4 3 2])))), [4 1 3 2]);
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
            N = obj.targets.numberOfTargets;
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

        %%% configuration

        function configure(obj, options)
            arguments
                obj
                options.noise (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.noise, [0, 1])} = obj.configuration.noise
                options.directPath (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.directPath, [0, 1])} = obj.configuration.directPath
            end
            obj.configuration.noise = options.noise;
            obj.configuration.directPath = options.directPath;
        end

        %%% visualization methods

        function visualizescenario(obj, options)
            arguments
                obj
                options.showVelocity (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.showVelocity, [0, 1])} = true
                options.showUnitNormalBaseline (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.showUnitNormalBaseline, [0, 1])} = false
                options.showUnitDirectionBisector (1, 1) logical {mustBeNumericOrLogical, mustBeMember(options.showUnitDirectionBisector, [0, 1])} = false
            end
            posRX = [obj.network.activeReceivingNodes.position]/1e3;
            posTX = [obj.network.activeTransmittingNodes.position]/1e3;
            figure; plot3(posRX(1, :), posRX(2, :), posRX(3, :), 'vb', 'LineWidth', 2, 'MarkerSize', 10);
            hold on; plot3(posTX(1, :), posTX(2, :), posTX(3, :), 'vr', 'LineWidth', 2, 'MarkerSize', 10);
            x = obj.targets.position(1, :)/1e3;
            y = obj.targets.position(2, :)/1e3;
            z = obj.targets.position(3, :)/1e3;
            vx = obj.targets.velocity(1, :)/1e3;
            vy = obj.targets.velocity(2, :)/1e3;
            vz = obj.targets.velocity(3, :)/1e3;
            if obj.numberOfTargets < 11
                plot3(x, y, z, '+k', 'LineWidth', 2, 'MarkerSize', 10);
            else
                plot3(x, y, z, '*k', 'LineWidth', 1, 'MarkerSize', 10);
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
            text(posRX(1, :), posRX(2, :), posRX(3, :), num2str(find(obj.network.receivingNodeActivity).'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            text(posTX(1, :), posTX(2, :), posTX(3, :), num2str(find(obj.network.transmittingNodeActivity).'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            if obj.numberOfTargets < 11
                text(x, y, z, num2str((1 : obj.numberOfTargets).'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            end
            posRX = repelem(posRX, 1, obj.numberOfTransmittingNodes);
            posTX = repmat(posTX, 1, obj.numberOfReceivingNodes);
            line([posRX(1, :); posTX(1, :)], [posRX(2, :); posTX(2, :)], [posRX(3, :); posTX(3, :)], 'lineStyle', '--', 'Color', 'k');
            grid on; grid minor; title('scenario');
            xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
            legend('RX', 'TX', 'targets', 'velocity', 'Location', 'best');
            if options.showVelocity
                if options.showUnitNormalBaseline
                    if options.showUnitDirectionBisector
                        legend('rx', 'tx', 'target', 'velocity', 'n_{baseline}', 'n_{bisector}', 'Location', 'best');
                    else
                        legend('rx', 'tx', 'target', 'velocity', 'n_{baseline}', 'Location', 'best');
                    end
                elseif options.showUnitDirectionBisector
                    legend('rx', 'tx', 'target', 'velocity', 'n_{bisector}', 'Location', 'best');
                end
            else
                if options.showUnitNormalBaseline
                    if options.showUnitDirectionBisector
                        legend('rx', 'tx', 'target', 'n_{baseline}', 'n_{bisector}', 'Location', 'best');
                    else
                        legend('rx', 'tx', 'target', 'n_{baseline}', 'Location', 'best');
                    end
                elseif options.showUnitDirectionBisector
                    legend('rx', 'tx', 'target', 'n_{bisector}', 'Location', 'best');
                end
            end
        end

        function visualizereceivedsignals(obj, options)
            arguments
                obj
                options.receivingNodeIDs (1, :) double {mustBeInteger, mustBeNonnegative} = 1
            end
            mustBeInRange(options.receivingNodeIDs, 1, obj.numberOfReceivingNodes);
            s = obj.signalSuperposed(:, options.receivingNodeIDs);
            figure; hold on;
            for rxID = options.receivingNodeIDs
                t = obj.network.activeReceivingNodes(rxID).samplingInstants;
                plot(t, 20*log10(abs(s)));
            end
            grid on; grid minor;
            xlabel('time (s)'); ylabel('power (dB)');
            title('received signal');
            leg = legend(num2str(options.receivingNodeIDs.'), 'Location', 'best');
            title(leg, 'RX ID'); hold off;
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
            axID1 = obj.targets.firstAxisID;
            axID2 = obj.targets.secondAxisID;
            L1 = obj.targets.firstGridLength;
            L2 = obj.targets.secondGridLength;
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
                figure; contourf(x, y, reshape(-obj.timeDelay(options.transmittingNodeID, options.receivingNodeID, :)*1e6, L1, L2), options.numberOfCurves, "ShowText", true, "LabelFormat", "%3.3g usec");
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
end