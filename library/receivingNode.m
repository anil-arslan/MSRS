classdef receivingNode < handle
    %receivingNode Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = private, GetAccess = public)
        position (3, 1) double = zeros(3, 1)
        systemLoss_dB (1, 1) double {mustBeNonnegative} = 0 % dB
        noiseFigure_dB (1, 1) double {mustBeNonnegative} = 0 % dB
        temperature (1, 1) double {mustBeNonnegative} = 300 % Kelvin
    end

    properties (Dependent)
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

    properties (Constant)
        speedOfLight (1, 1) double {mustBePositive} = physconst('LightSpeed');
        constantBoltzmann (1, 1) double {mustBePositive} = physconst('Boltzmann');
    end

    methods
        function obj = receivingNode(options)
            arguments
                options.position (3, :) double = zeros(3, 1)
                options.systemLoss (1, :) double {mustBeNonnegative} = 0 % dB
                options.samplingFrequency (1, :) double {mustBeNonnegative} = 4e7 % Hz
                options.temperature (1, :) double {mustBeNonnegative} = 300 % Kelvin
                options.noiseFigure (1, :) double {mustBeNonnegative} = 0 % dB
                options.CPIsecond (1, :) double {mustBeNonnegative} = 20e-5 % sec
            end
            %receivingNode Construct an instance of this class
            %   Detailed explanation goes here
            numberOfNodes = max([size(options.position, 2), numel(options.systemLoss), numel(options.samplingFrequency), ...
                numel(options.temperature), numel(options.noiseFigure), numel(options.CPIsecond)]);
            if numberOfNodes == 1
                obj.position = options.position;
                obj.systemLoss_dB = options.systemLoss;
                obj.samplingFrequency = options.samplingFrequency;
                obj.temperature = options.temperature;
                obj.noiseFigure_dB = options.noiseFigure;
                obj.CPI = options.CPIsecond;
            else
                if size(options.position, 2) == 1
                    options.position = repmat(options.position, 1, numberOfNodes);
                elseif size(options.position, 2) ~= numberOfNodes
                    error('number of receiving nodes is %d', numberOfNodes);
                end
                for field = ["systemLoss", "samplingFrequency", "temperature", "noiseFigure", "CPIsecond"]
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
                        'systemLoss', options.systemLoss(nodeID), ...
                        'samplingFrequency', options.samplingFrequency(nodeID), ...
                        'temperature', options.temperature(nodeID), ...
                        'noiseFigure', options.noiseFigure(nodeID), ...
                        'CPIsecond', options.CPIsecond(nodeID));
                end
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
    end
end