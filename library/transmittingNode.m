classdef transmittingNode < handle & dynamicprops
    %transmittingNode Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = private, GetAccess = public)
        array planarArray {mustBeScalarOrEmpty} = planarArray.empty()
        position (3, 1) double = zeros(3, 1)
        inputPower_W (1, 1) double {mustBePositive} = 1e3 % Watt
    end

    properties (Dependent)
        transmittedEnergy (1, 1) double {mustBePositive} % Joule
        transmissionType (1, 1) string {mustBeMember(transmissionType, ["pulsed", "continuous"])}
    end

    properties (SetAccess = private, GetAccess = public)
        carrierFrequency (1, 1) double {mustBeNonnegative} = 1e9 % Hz
    end

    properties (Dependent)
        carrierWavelength (1, 1) double {mustBeNonnegative} % m
        peakTransmitGain (1, 1) double {mustBeNonnegative} % power scale
    end

    properties (Dependent)
        steeringVector (:, 1) double
        beamformer (:, 1) double
        beamformingGain_dB (1, 1) double {mustBeNonnegative} % dB
    end

    properties (SetAccess = private, GetAccess = public)
        numberOfPulses (1, 1) double {mustBeInteger, mustBePositive} = 1
        dutyCycle (1, 1) double {mustBeInRange(dutyCycle, 0, 1)} = 0.1
        pulseWidth (1, 1) double {mustBeNonnegative} = 1e-6 % sec
    end

    properties (Dependent)
        PRI (1, 1) double {mustBeNonnegative} % sec
        PRF (1, 1) double {mustBeNonnegative} % Hertz
        beamTime (1, 1) double {mustBeNonnegative} % sec
    end

    properties (SetAccess = private, GetAccess = public)
        modulationType (1, 1) string {mustBeMember(modulationType, ["unmodulated", "linearFrequencyModulated"])} = "unmodulated"
        taperTypeFastTime (1, 1) string {mustBeMember(taperTypeFastTime, ["rectwin", "hann", "hamming"])} = "rectwin"
        taperTypeSpatial (1, 1) string {mustBeMember(taperTypeSpatial, ["rectwin", "hann", "hamming", "taylorwin"])} = "taylorwin"
    end

    properties (Dependent)
        taperSpatial (:, 1) double
        waveform (1, 1) function_handle
    end

    properties (Constant)
        speedOfLight (1, 1) double {mustBePositive} = physconst('LightSpeed');
    end

    properties (Access = private, Hidden)
        dynamicProperties (1, 1) struct
    end

    methods
        function obj = transmittingNode(options)
            arguments
                options.position (3, :) double = zeros(3, 1)
                options.array (1, 1) planarArray = planarArray
                options.inputPower_W (1, :) double {mustBePositive} = 1e3 % Watt
                options.carrierFrequency (1, :) double {mustBePositive} = 1e9 % Hz
                options.pulseWidth (1, :) double {mustBeNonnegative} = 1e-6 % sec
            end
            %transmittingNode Construct an instance of this class
            %   Detailed explanation goes here
            numberOfNodes = max([size(options.position, 2), numel(options.array), numel(options.inputPower_W), numel(options.carrierFrequency), ...
                numel(options.pulseWidth)]);
            if numberOfNodes == 1
                obj.position = options.position;
                obj.inputPower_W = options.inputPower_W;
                obj.carrierFrequency = options.carrierFrequency;
                obj.pulseWidth = options.pulseWidth;
                obj.setunmodulation;
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
                    error('number of transmitting nodes is %d', numberOfNodes);
                end
                for fieldName = ["array", "inputPower_W", "carrierFrequency", "pulseWidth"]
                    if isscalar(options.(fieldName))
                        options.(fieldName) = repmat(options.(fieldName), 1, numberOfNodes);
                    elseif numel(options.(fieldName)) ~= numberOfNodes
                        error('number of transmitting nodes is %d', numberOfNodes);
                    end
                end
                obj = transmittingNode.empty(0, numberOfNodes);
                for nodeID = 1 : numberOfNodes
                    obj(nodeID) = transmittingNode( ...
                        'position', options.position(:, nodeID), ...
                        'array', options.array(nodeID), ...
                        'inputPower_W', options.inputPower_W(nodeID), ...
                        'carrierFrequency', options.carrierFrequency(nodeID), ...
                        'pulseWidth', options.pulseWidth(nodeID));
                end
            end
        end

        %%% get methods

        function E = get.transmittedEnergy(obj)
            E = obj.inputPower_W.*obj.pulseWidth.*obj.numberOfPulses;
        end

        function out = get.transmissionType(obj)
            if isinf(obj.pulseWidth)
                out = "continuous";
            else
                out = "pulsed";
            end
        end

        function lambda = get.carrierWavelength(obj)
            lambda = obj.speedOfLight/obj.carrierFrequency;
        end

        function Gt = get.peakTransmitGain(obj)
            Gt = 4*pi*obj.array.apertureArea./obj.carrierWavelength.^2;
        end

        function a = get.steeringVector(obj)
            a = exp(1j*2*pi*obj.array.steeringPositions./obj.carrierWavelength);
        end

        function h = get.beamformer(obj)
            h = obj.taperSpatial.*obj.steeringVector;
        end

        function G = get.beamformingGain_dB(obj)
            G = 20*log10(sum(obj.taperSpatial));
        end

        function pri = get.PRI(obj)
            pri = obj.pulseWidth./obj.dutyCycle;
        end

        function prf = get.PRF(obj)
            prf = 1./obj.PRI;
        end

        function t = get.beamTime(obj)
            t = obj.PRI.*obj.numberOfPulses;
        end

        function T = get.taperSpatial(obj)
            funcTaper = eval("@" + obj.taperTypeSpatial);
            T = funcTaper(obj.array.numberOfTotalElements);
            T = T./norm(T);
        end

        function s = get.waveform(obj)
            u = @(t) double(t > 0); % unit step function
            switch obj.transmissionType
                case "pulsed"
                    r = @(t) u(t + obj.pulseWidth) - u(t);
                case "continuous"
                    r = @(t) u(-t);
            end
            switch obj.modulationType % complex modulation waveforms
                case "unmodulated"
                    s = @(t) obj.taperFastTimeSamples(t).*r(t);
                case "linearFrequencyModulated"
                    s = @(t) LFM(t);
            end
            function sig = LFM(t)
                sig = zeros(size(t));
                finiteInstants = ~isinf(t);
                t = t(finiteInstants);
                switch obj.frequencyDirection
                    case "increasing"
                        sig(finiteInstants) = obj.taperFastTimeSamples(t).*r(t).*exp(1j*pi*(obj.bandWidth/obj.pulseWidth)*(t + obj.pulseWidth).^2 + 1j*2*pi*obj.frequencyOffset*(t + obj.pulseWidth));
                    case "decreasing"
                        sig(finiteInstants) = obj.taperFastTimeSamples(t).*r(t).*exp(-1j*pi*(obj.bandWidth/obj.pulseWidth)*t.^2 + 1j*2*pi*obj.frequencyOffset.*t);
                end
            end
        end

        function T = taperFastTimeSamples(obj, samplingInstants)
            pulseInstants = samplingInstants <= 0 & samplingInstants > -permute([obj.pulseWidth], [1 4 3 2]);
            numberOfReceivingNodes = size(samplingInstants, 2);
            numberOfTransmittingNodes = numel(obj);
            numberOfTargets = size(samplingInstants, 3);
            T = zeros(size(samplingInstants, 1), numberOfTargets, numberOfReceivingNodes, numberOfTransmittingNodes);
            for rxID = 1 : numberOfReceivingNodes
                for txID = 1: numberOfTransmittingNodes
                    switch obj(txID).taperTypeFastTime
                        case "rectwin"
                            T = 1;
                        otherwise
                            targetInstants = any(pulseInstants(:, rxID, :, txID), 1);
                            t = eval(obj(txID).taperTypeFastTime + "(max(sum(pulseInstants(:, rxID, :, txID))))");
                            K = size(samplingInstants, 1)*numberOfTargets;
                            idx = find(pulseInstants(:, rxID, targetInstants, txID)) + (rxID - 1)*K + (txID - 1)*(rxID - 1)*K;
                            T(idx) = repmat(t, sum(targetInstants), 1);
                    end
                end
            end
            T = permute(T, [1 3 2 4]);
        end

        %%% set methods

        function settaper(obj, options)
            arguments
                obj
                options.taperTypeFastTime (1, :) string {mustBeMember(options.taperTypeFastTime, ["rectwin", "hann", "hamming"])} = "rectwin"
                options.taperTypeSpatial (1, :) string {mustBeMember(options.taperTypeSpatial, ["rectwin", "hann", "hamming"])} = "rectwin"
            end
            numberOfTransmittingNodes = numel(obj);
            if isscalar(options.taperTypeFastTime)
                options.taperTypeFastTime = repmat(options.taperTypeFastTime, 1, numberOfTransmittingNodes);
            elseif numel(options.taperTypeFastTime) ~= numberOfTransmittingNodes
                error('taperTypeFastTime vector must have a size of either %d or %d', numberOfTransmittingNodes, 1);
            end
            if isscalar(options.taperTypeSpatial)
                options.taperTypeSpatial = repmat(options.taperTypeSpatial, 1, numberOfTransmittingNodes);
            elseif numel(options.taperTypeSpatial) ~= numberOfTransmittingNodes
                error('taperTypeSpatial vector must have a size of either %d or %d', numberOfTransmittingNodes, 1);
            end
            for txID = 1 : numberOfTransmittingNodes
                obj(txID).taperTypeFastTime = options.taperTypeFastTime(txID);
                obj(txID).taperTypeSpatial = options.taperTypeSpatial(txID);
            end
        end

        function settransmission(obj, options)
            arguments
                obj
                options.inputPower_W (1, :) double {mustBePositive} = [] % W
                options.carrierFrequency (1, :) double {mustBePositive} = [] % Hz
                options.pulseWidth (1, :) double {mustBeNonnegative} = [] % sec
                options.dutyCycle (1, :) double {mustBeNonnegative} = [] % sec
            end
            numberOfTransmittingNodes = numel(obj);
            fieldNames = ["inputPower_W", "carrierFrequency", "pulseWidth", "dutyCycle"];
            for fieldName = fieldNames
                if ~isempty(options.(fieldName))
                    if isscalar(options.(fieldName))
                        options.(fieldName) = repmat(options.(fieldName), 1, numberOfTransmittingNodes);
                    elseif size(options.(fieldName), 2) ~= numberOfTransmittingNodes
                        error('length of %s transmitting nodes is %d', fieldName, numberOfTransmittingNodes);
                    end
                else
                    options.(fieldName) = [obj.(fieldName)];
                end
            end
            for txID = 1 : numberOfTransmittingNodes
                for fieldName = fieldNames
                    obj(txID).(fieldName) = options.(fieldName)(txID);
                end
            end
        end

        function setunmodulation(obj)
            numberOfTransmittingNodes = numel(obj);
            fieldNames = ["frequencyOffset", "bandWidth", "frequencyDirection", "timeBandwidthProduct"];
            for txID = 1 : numberOfTransmittingNodes
                for fieldName = fieldNames
                    if isprop(obj(txID), fieldName)
                        delete(obj(txID).dynamicProperties.(fieldName));
                        obj(txID).dynamicProperties = rmfield(obj(txID).dynamicProperties, fieldName);
                    end
                end
                obj(txID).modulationType = "unmodulated";
            end
        end

        function setLFM(obj, options)
            arguments
                obj
                options.frequencyOffset (1, :) double {mustBeNonnegative} = 0 % Hz
                options.bandWidth (1, :) double {mustBeNonnegative} = 5e6 % Hz
                options.frequencyDirection (1, :) string {mustBeMember(options.frequencyDirection, ["increasing", "decreasing"])} = "increasing"
            end
            numberOfTransmittingNodes = numel(obj);
            fieldNames = ["frequencyOffset", "bandWidth", "frequencyDirection"];
            dependentPropertyNames = "timeBandwidthProduct";
            for fieldName = fieldNames
                if isscalar(options.(fieldName))
                    options.(fieldName) = repmat(options.(fieldName), 1, numberOfTransmittingNodes);
                elseif numel(options.(fieldName)) ~= numberOfTransmittingNodes
                    error('length of %s vector must be 1 or %d', fieldname, numberOfTransmittingNodes);
                end
            end
            for txID = 1 : numberOfTransmittingNodes
                if ~strcmpi(obj(txID).modulationType, "linearFrequencyModulated")
                    for fieldName = fieldNames
                        if ~isprop(obj(txID), fieldName)
                            obj(txID).dynamicProperties.(fieldName) = obj(txID).addprop(fieldName);
                            obj(txID).dynamicProperties.(fieldName).SetAccess = 'private';
                        end
                        obj(txID).(fieldName) = options.(fieldName)(txID);
                    end
                    for fieldName = dependentPropertyNames
                        if ~isprop(obj(txID), fieldName)
                            obj(txID).dynamicProperties.(fieldName) = obj(txID).addprop(fieldName);
                            obj(txID).dynamicProperties.(fieldName).Dependent = true;
                            switch fieldName
                                case "timeBandwidthProduct"
                                    obj(txID).dynamicProperties.timeBandwidthProduct.GetMethod = @(obj) obj.pulseWidth.*obj.bandWidth;
                            end
                        end
                    end
                end
                obj(txID).modulationType = "linearFrequencyModulated";
            end
        end

        %%% utility methods

        function beamSpaceSignal = beamform(obj, elementSpaceSignal)
            beamSpaceSignal = pagemtimes(elementSpaceSignal, conj(obj.beamformer));
        end

        %%% visualization methods

        function visualizewaveform(obj, options)
            arguments
                obj
                options.domain (1, 1) string {mustBeMember(options.domain, ["time", "frequency"])} = "time"
                options.plot (1, 1) string {mustBeMember(options.plot, ["magnitude", "phase", "real", "imaginary"])} = "magnitude"
            end
            figure; hold on;
            numberOfTransmittingNodes = numel(obj);
            for txID = 1 : numberOfTransmittingNodes
                switch options.domain
                    case "time"
                        Ns = 1024;
                        if isinf(obj(txID).pulseWidth)
                            ts = 1/obj(txID).frequencyOffset;
                            t = -Ns*ts : ts : 0;
                        else
                            ts = obj(txID).pulseWidth/Ns;
                            t = -obj(txID).pulseWidth + ts : ts : 0;
                        end
                        switch options.plot
                            case "magnitude"
                                plot(t*1e6, abs(obj(txID).waveform(t)));
                            case "phase"
                                plot(t*1e6, 180*unwrap(angle(obj(txID).waveform(t)))/pi);
                            case "real"
                                plot(t*1e6, real(obj(txID).waveform(t)));
                            case "imaginary"
                                plot(t*1e6, imag(obj(txID).waveform(t)));
                        end
                    case "frequency"
                        Nsover = 4;
                        if isprop(obj(txID), 'timeBandwidthProduct')
                            Ns = 2.^(nextpow2(obj(txID).timeBandwidthProduct) + 1);
                        else
                            Ns = 16;
                        end
                        ts = obj(txID).pulseWidth/Ns;
                        t = -obj(txID).pulseWidth + ts : ts : 0;
                        nfft = Ns*Nsover;
                        f = (-nfft/2 : nfft/2 - 1)/ts/nfft/1e6;
                        switch options.plot
                            case "magnitude"
                                plot(f, 20*log10(fftshift(abs(fft(obj(txID).waveform(t), nfft)))));
                            case "phase"
                                plot(f, 180*unwrap(fftshift(angle(fft(obj(txID).waveform(t), nfft))))/pi);
                            case "real"
                                plot(f, 20*log10(fftshift(real(fft(obj(txID).waveform(t), nfft)))));
                            case "imaginary"
                                plot(f, 20*log10(fftshift(imag(fft(obj(txID).waveform(t), nfft)))));
                        end
                end
            end
            grid on; grid minor;
            switch options.domain
                case "time"
                    xlabel('time (us)');
                    strTitle = 'waveform in time domain';
                case "frequency"
                    xlabel('frequency (MHz)');
                    strTitle = 'waveform in frequency domain';
            end
            switch options.plot
                case "magnitude"
                    ylabel('power (dB)');
                    strTitle = ['magnitude of ' strTitle];
                case "phase"
                    ylabel('phase (°)');
                    strTitle = ['phase of ' strTitle];
                case "real"
                    ylabel('real part');
                    strTitle = ['real part of ' strTitle];
                case "imaginary"
                    ylabel('imaginary part');
                    strTitle = ['imaginary part of ' strTitle];
            end
            title(strTitle);
            leg = legend(num2str((1 : numberOfTransmittingNodes).'), 'Location', 'best');
            title(leg, 'TX ID');
            hold off; drawnow;
        end
    end
end