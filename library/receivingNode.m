classdef receivingNode < handle
    %receivingNode Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = private, GetAccess = public)
        beamformingMode (1, 1) string {mustBeMember(beamformingMode, ["conventional", "bypass"])} = "conventional"
        array planarArray {mustBeScalarOrEmpty} = planarArray.empty()
        position (3, 1) double = zeros(3, 1)
        systemLoss_dB (1, 1) double {mustBeNonnegative} = 0 % dB
        noiseFigure_dB (1, 1) double {mustBeNonnegative} = 0 % dB
        temperature (1, 1) double {mustBeNonnegative} = 300 % Kelvin
    end

    properties (SetAccess = private, GetAccess = public)
        beamCentersElevation (1, :) double {mustBeInRange(beamCentersElevation, -90, 90)} = 0
        beamCentersAzimuth (1, :) double {mustBeInRange(beamCentersAzimuth, -180, 180)} = 0
        taperTypeSpatial (1, 1) string {mustBeMember(taperTypeSpatial, ["rectwin", "hann", "hamming", "taylorwin"])} = "taylorwin"
    end

    properties (Dependent)
        taperSpatial (:, 1) double
        beamformer double % M x Nch matrix
        beamformingGain_dB (1, 1) double {mustBeNonnegative} % dB
        beamCenterUnitVector (2, :) double
        beamCenterUnitDirection (3, :) double
        beamCenterPositions double
        numberOfTotalChannels {mustBeInteger, mustBePositive}
        noisePowerPerSample_dB (1, 1) double % dB watt
    end

    properties (SetAccess = private, GetAccess = public)
        CPI (1, 1) double {mustBeNonnegative} = 20e-5 % sec
    end

    properties (Dependent)
        listenedRadius (1, 1) double {mustBeNonnegative} % meter
    end

    properties (SetAccess = private, GetAccess = public)
        samplingFrequency (1, 1) double {mustBeNonnegative} = 4e7 % Hz
    end

    properties (Dependent)
        samplingPeriod (1, 1) double {mustBeNonnegative} % second
        samplingInstants (:, 1) double % sec (Ns x 1) vector
        numberOfSamplesPerCPI (1, 1) double {mustBeInteger, mustBeNonnegative} % during CPI
    end

    properties (GetAccess = public, SetAccess = ?radarNetwork)
        directionFinder struct = struct( ...
            'table', [], ...
            'elevation', [], ...
            'azimuth', [], ...
            'azimuthResolution', [], ...
            'elevationResolution', [])
    end

    properties (Constant)
        speedOfLight (1, 1) double {mustBePositive} = physconst('LightSpeed');
        constantBoltzmann (1, 1) double {mustBePositive} = physconst('Boltzmann');
    end

    methods
        function obj = receivingNode(options)
            arguments
                options.position (3, :) double = zeros(3, 1)
                options.array (1, :) planarArray = planarArray
                options.systemLoss (1, :) double {mustBeNonnegative} = 0 % dB
                options.samplingFrequency (1, :) double {mustBeNonnegative} = 4e7 % Hz
                options.temperature (1, :) double {mustBeNonnegative} = 300 % Kelvin
                options.noiseFigure (1, :) double {mustBeNonnegative} = 0 % dB
                options.CPIsecond (1, :) double {mustBeNonnegative} = 20e-5 % sec
            end
            %receivingNode Construct an instance of this class
            %   Detailed explanation goes here
            numberOfNodes = max([size(options.position, 2), numel(options.array), numel(options.systemLoss), numel(options.samplingFrequency), ...
                numel(options.temperature), numel(options.noiseFigure), numel(options.CPIsecond)]);
            if numberOfNodes == 1
                obj.position = options.position;
                obj.systemLoss_dB = options.systemLoss;
                obj.samplingFrequency = options.samplingFrequency;
                obj.temperature = options.temperature;
                obj.noiseFigure_dB = options.noiseFigure;
                obj.CPI = options.CPIsecond;
                if isempty(options.array.node)
                    obj.array = options.array;
                else
                    obj.array = planarArray( ...
                        "numberOfElements", options.array.numberOfElements, ...
                        "spacing", options.array.spacing, ...
                        "rpm", options.array.rpm, ...
                        "scanDirection", options.array.scanDirection, ...
                        "backOfArray", options.array.backOfArray);
                end
                obj.array.node = obj;
            else
                if size(options.position, 2) == 1
                    options.position = repmat(options.position, 1, numberOfNodes);
                elseif size(options.position, 2) ~= numberOfNodes
                    error('number of receiving nodes is %d', numberOfNodes);
                end
                for field = ["array", "systemLoss", "samplingFrequency", "temperature", "noiseFigure", "CPIsecond"]
                    if isscalar(options.(field))
                        options.(field) = repmat(options.(field), 1, numberOfNodes);
                    elseif size(options.(field), 2) ~= numberOfNodes
                        error('number of receiving nodes is %d', numberOfNodes);
                    end
                end
                obj = receivingNode.empty(0, numberOfNodes);
                for nodeID = 1 : numberOfNodes
                    obj(nodeID) = receivingNode( ...
                        'position', options.position(:, nodeID), ...
                        'array', options.array(nodeID), ...
                        'systemLoss', options.systemLoss(nodeID), ...
                        'samplingFrequency', options.samplingFrequency(nodeID), ...
                        'temperature', options.temperature(nodeID), ...
                        'noiseFigure', options.noiseFigure(nodeID), ...
                        'CPIsecond', options.CPIsecond(nodeID));
                end
            end
        end

        function T = get.taperSpatial(obj)
            funcTaper = eval("@" + obj.taperTypeSpatial);
            T = funcTaper(obj.array.numberOfTotalElements);
            T = T./norm(T);
        end

        function H = get.beamformer(obj)
            H = obj.taperSpatial;
        end

        function G = get.beamformingGain_dB(obj)
            switch obj.beamformingMode
                case 'conventional'
                    G = 20*log10(sum(obj.taperSpatial, 1));
                case 'bypass'
                    G = 10*log10([obj.array.numberOfTotalElements].');
            end
        end

        function g = get.beamCenterUnitVector(obj)
            g = [cosd(obj.beamCentersElevation).*sind(obj.beamCentersAzimuth); sind(obj.beamCentersElevation)];
        end

        function g = get.beamCenterUnitDirection(obj)
            g = obj.array.rotationMatrix*[cosd(obj.beamCentersElevation).*cosd(obj.beamCentersAzimuth); cosd(obj.beamCentersElevation).*sind(obj.beamCentersAzimuth); sind(obj.beamCentersElevation)];
        end

        function p = get.beamCenterPositions(obj)
            p = obj.array.positionsYZ.'*obj.beamCenterUnitVector;
        end

        function Nch = get.numberOfTotalChannels(obj)
            switch obj.beamformingMode
                case 'conventional'
                    Nch = numel(obj.beamCentersElevation);
                case 'bypass'
                    Nch = obj.array.numberOfTotalElements;
            end
        end

        function Pn = get.noisePowerPerSample_dB(obj)
            Pn = 10*log10(obj.samplingFrequency*obj.temperature*obj.constantBoltzmann) + obj.noiseFigure_dB;
        end

        function R = get.listenedRadius(obj)
            R = obj.CPI.*obj.speedOfLight/2;
        end

        function Ts = get.samplingPeriod(obj)
            Ts = 1/obj.samplingFrequency;
        end

        function t = get.samplingInstants(obj)
            t = (1/obj.samplingFrequency : 1/obj.samplingFrequency : obj.CPI).';
        end

        function N = get.numberOfSamplesPerCPI(obj)
            N = round(obj.CPI.*obj.samplingFrequency);
        end

        %%% set methods

        function setbeamcenters(obj,options)
            arguments
                obj
                options.beamCentersElevation (1, :) = {obj.beamCentersElevation}
                options.beamCentersAzimuth (1, :) = {obj.beamCentersAzimuth}
            end
            numberOfReceivingNodes = numel(obj);
            for rxID = 1 : numberOfReceivingNodes
                if iscell(options.beamCentersElevation)
                    obj(rxID).beamCentersElevation = options.beamCentersElevation{rxID};
                else
                    obj(rxID).beamCentersElevation = options.beamCentersElevation;
                end
                if iscell(options.beamCentersAzimuth)
                    obj(rxID).beamCentersAzimuth = options.beamCentersAzimuth{rxID};
                else
                    obj(rxID).beamCentersAzimuth = options.beamCentersAzimuth;
                end
                if isscalar(obj(rxID).beamCentersElevation) && ~isscalar(obj(rxID).beamCentersAzimuth)
                    obj(rxID).beamCentersElevation = repmat(obj(rxID).beamCentersElevation, 1, length(obj(rxID).beamCentersAzimuth));
                elseif isscalar(obj(rxID).beamCentersAzimuth) && ~isscalar(obj(rxID).beamCentersElevation)
                    obj(rxID).beamCentersAzimuth = repmat(obj(rxID).beamCentersAzimuth, 1, length(obj(rxID).beamCentersElevation));
                elseif length(obj(rxID).beamCentersAzimuth) ~= length(obj(rxID).beamCentersElevation)
                    error('length of beam elevation and azimuth centers must be the same');
                end
            end
        end

        function settaper(obj, options)
            arguments
                obj
                options.taperTypeSpatial (1, :) string {mustBeMember(options.taperTypeSpatial, ["rectwin", "hann", "hamming"])} = "rectwin"
            end
            numberOfReceivingNodes = numel(obj);
            if isscalar(options.taperTypeSpatial)
                options.taperTypeSpatial = repmat(options.taperTypeSpatial, 1, numberOfReceivingNodes);
            elseif numel(options.taperTypeSpatial) ~= numberOfReceivingNodes
                error('taperTypeSpatial vector must have a size of either %d or %d', numberOfReceivingNodes, 1);
            end
            for rxID = 1 : numberOfReceivingNodes
                obj(rxID).taperTypeSpatial = options.taperTypeSpatial(rxID);
            end
        end

        function settingsreceiver(obj, options)
            arguments
                obj
                options.beamformingMode (1, 1) string {mustBeMember(options.beamformingMode, ["conventional", "bypass"])} = obj.beamformingMode
            end
            numberOfReceivingNodes = numel(obj);
            if isscalar(options.beamformingMode)
                options.beamformingMode = repmat(options.beamformingMode, 1, numberOfReceivingNodes);
            elseif numel(options.beamformingMode) ~= numberOfReceivingNodes
                error('beamformingMode vector must have a size of either %d or %d', numberOfReceivingNodes, 1);
            end
            for rxID = 1 : numberOfReceivingNodes
                obj(rxID).beamformingMode = options.beamformingMode(rxID);
            end
        end

        %%% utility methods

        function beamSpaceSignal = beamform(obj, elementSpaceSignal, wavelength)
            switch obj.beamformingMode
                case 'conventional'
                    beamformingCenters = exp(1j*2*pi*obj.beamCenterPositions./wavelength);
                    steeringVector = exp(1j*2*pi*obj.array.steeringPositions./wavelength);
                    beamSpaceSignal = pagemtimes(elementSpaceSignal, conj(obj.beamformer.*beamformingCenters.*steeringVector));
                case 'bypass'
                    beamSpaceSignal = elementSpaceSignal;
            end
        end
    end
end