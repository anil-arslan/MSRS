classdef uniformPlanarArray < handle
    %uniformPlanarArray Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = ?radarNetwork, GetAccess = public)
        positionsYZ (2, :) double = zeros(2, 1) % meters y-z plane
        roll (1, 1) double {mustBeInRange(roll, -90, 90)} = 0 % degrees
        pitch (1, 1) double {mustBeInRange(pitch, -90, 90)} = 0 % degrees
        yaw (1, 1) double {mustBeInRange(yaw, -180, 180)} = 0 % degrees
        rpm (1, 1) double {mustBeNonnegative} = 0 % rev per min
        scanDirection (1, 1) string {mustBeMember(scanDirection, ["cw", "ccw"])} = "ccw"
        backOfArray (1, 1) logical {mustBeNumericOrLogical} = false
        backOfArrayRegion (1, 1) double {mustBeInRange(backOfArrayRegion, 0, 360)} = 180
    end

    properties (Dependent)
        positions (3, :) double % meters y-z plane
    end

    properties (SetAccess = ?radarNetwork, GetAccess = public)
        spacing (2, 1) double {mustBeNonnegative} = zeros(2, 1)
    end

    properties (SetAccess = private, GetAccess = public)
        numberOfElements (2, 1) double {mustBePositive, mustBeInteger} = ones(2, 1)
    end

    properties (Dependent)
        numberOfTotalElements (1, 1) double {mustBePositive, mustBeInteger}
        apertureArea (1, 1) double
        rotationMatrix (3, 3) double
        normalVector (3, 1) double
    end

    properties (SetAccess = {?receivingNode, ?transmittingNode}, GetAccess = public)
        steeringElevation (1, 1) double {mustBeInRange(steeringElevation, -90, 90)} = 0
        steeringAzimuth (1, 1) double {mustBeInRange(steeringAzimuth, -180, 180)} = 0
        steeringAzimuthStep (1, 1) double {mustBeInRange(steeringAzimuthStep, 0, 180)} = 3
        steeringAzimuthLimit (1, 1) double {mustBeInRange(steeringAzimuthLimit, 0, 180)} = 60
    end

    properties (Dependent)
        steeringUnitVector (2, 1) double
        steeringUnitDirection (3, 1) double
        steeringPositions (:, 1) double
    end

    properties (SetAccess = {?receivingNode, ?transmittingNode}, GetAccess = public)
        node {mustBeScalarOrEmpty} = receivingNode.empty()
    end

    properties (Access = ?radarNetwork)
        spacingSpecified logical {mustBeNumericOrLogical} = true(2, 1)
    end

    properties (Access = private)
        durationSimulation (1, 1) {mustBeNonnegative} = 0
        steeringAzimuthInitial (1, 1) double {mustBeInRange(steeringAzimuthInitial, -180, 180)} = 0
    end

    methods
        function obj = uniformPlanarArray(options)
            %uniformPlanarArray Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                options.numberOfElements (2, 1) double = ones(2, 1)
                options.spacing (2, 1) double {mustBeNonnegative} = zeros(2, 1)
                options.rpm (1, 1) double {mustBeNonnegative} = 0 % rev per min
                options.scanDirection (1, 1) string {mustBeMember(options.scanDirection, ["cw", "ccw"])} = "ccw"
                options.backOfArray (1, 1) logical {mustBeNumericOrLogical} = false
				options.backOfArrayRegion (1, 1) double {mustBeInRange(options.backOfArrayRegion, 0, 360)} = 180
            end
            obj.numberOfElements = options.numberOfElements;
            obj.spacing = options.spacing;
            if ~obj.spacing(1)
                obj.spacingSpecified(1) = false;
            end
            if ~obj.spacing(2)
                obj.spacingSpecified(2) = false;
            end
            yPos = (-(obj.numberOfElements(1) - 1)/2 : (obj.numberOfElements(1) - 1)/2)*obj.spacing(1);
            zPos = (-(obj.numberOfElements(2) - 1)/2 : (obj.numberOfElements(2) - 1)/2)*obj.spacing(2);
            [yPos, zPos] = meshgrid(yPos, zPos);
            obj.positionsYZ = zeros(2, prod(options.numberOfElements));
            obj.positionsYZ(1, :) = yPos(:).';
            obj.positionsYZ(2, :) = zPos(:).';
            obj.rpm = options.rpm;
            obj.scanDirection = options.scanDirection;
            obj.backOfArray = options.backOfArray;
            obj.backOfArrayRegion = options.backOfArrayRegion;
        end

        %%% get methods

        function pos = get.positions(obj)
            pos = obj.rotationMatrix*[zeros(1, obj.numberOfTotalElements); obj.positionsYZ];
        end

        function M = get.numberOfTotalElements(obj)
            M = prod(obj.numberOfElements);
        end

        function A = get.apertureArea(obj)
            A = obj.numberOfTotalElements.*prod(obj.spacing);
        end

        function R = get.rotationMatrix(obj)
            Ryaw = [cosd(obj.yaw), -sind(obj.yaw), 0; sind(obj.yaw), cosd(obj.yaw), 0; 0, 0, 1];
            Rpitch = [cosd(obj.pitch), 0, sind(obj.pitch); 0, 1, 0; -sind(obj.pitch), 0, cosd(obj.pitch)];
            Rroll = [1, 0, 0; 0, cosd(obj.roll), -sind(obj.roll); 0, sind(obj.roll), cosd(obj.roll)];
            R = Ryaw*Rpitch*Rroll;
        end

        function n = get.normalVector(obj)
            n = obj.rotationMatrix(:, 1);
        end

        function g = get.steeringUnitVector(obj)
            g = [cosd(obj.steeringElevation).*sind(obj.steeringAzimuth); sind(obj.steeringElevation)];
        end

        function g = get.steeringUnitDirection(obj)
            g = obj.rotationMatrix*[cosd(obj.steeringElevation).*cosd(obj.steeringAzimuth); obj.steeringUnitVector];
        end

        function p = get.steeringPositions(obj)
            p = obj.positionsYZ.'*obj.steeringUnitVector;
        end

        %%% set methods

        function setorientation(obj, options)
            arguments
                obj
                options.roll (1, :) double {mustBeInRange(options.roll, -90, 90)} = [obj.roll] % degrees
                options.pitch (1, :) double {mustBeInRange(options.pitch, -90, 90)} = [obj.pitch] % degrees
                options.yaw (1, :) double {mustBeInRange(options.yaw, -180, 180)} = [obj.yaw] % degrees
                options.scanDirection (1, :) string {mustBeMember(options.scanDirection, ["cw", "ccw"])} = [obj.scanDirection]
            end
            numberOfArrays = numel(obj);
            options = checklength(options, numberOfArrays, 'number of arrays is');
            for arrayID = 1 : numberOfArrays
                obj(arrayID).roll = options.roll(arrayID);
                obj(arrayID).pitch = options.pitch(arrayID);
                obj(arrayID).yaw = options.yaw(arrayID);
                obj(arrayID).scanDirection = options.scanDirection(arrayID);
            end
        end

        function setscanparameters(obj, options)
            arguments
                obj
                options.rpm (1, :) double {mustBeNonnegative} = [obj.rpm] % rev per min
                options.steeringAzimuthStep (1, :) double {mustBeInRange(options.steeringAzimuthStep, 0, 180)} = [obj.steeringAzimuthStep]
                options.steeringAzimuthLimit (1, :) double {mustBeInRange(options.steeringAzimuthLimit, 0, 180)} = [obj.steeringAzimuthLimit]
                options.backOfArray (1, :) logical {mustBeNumericOrLogical} = [obj.backOfArray]
				options.backOfArrayRegion (1, :) double {mustBeInRange(options.backOfArrayRegion, 0, 360)} = [obj.backOfArrayRegion]
            end
            numberOfArrays = numel(obj);
            options = checklength(options, numberOfArrays, 'number of arrays is');
            for arrayID = 1 : numberOfArrays
                obj(arrayID).rpm = options.rpm(arrayID);
                obj(arrayID).steeringAzimuthStep = options.steeringAzimuthStep(arrayID);
                obj(arrayID).steeringAzimuthLimit = options.steeringAzimuthLimit(arrayID);
                obj(arrayID).backOfArray = options.backOfArray(arrayID);
                obj(arrayID).backOfArrayRegion = options.backOfArrayRegion(arrayID);
            end
        end

        function steer(obj, options)
            arguments
                obj
                options.steeringElevation (1, 1) double {mustBeInRange(options.steeringElevation, -90, 90)} = obj.steeringElevation
                options.steeringAzimuth (1, 1) double {mustBeInRange(options.steeringAzimuth, -90, 90)} = obj.steeringAzimuth
            end
            obj.steeringElevation = options.steeringElevation;
            obj.steeringAzimuth = options.steeringAzimuth;
            obj.steeringAzimuthInitial = options.steeringAzimuth;
        end

        function step(obj, timeStep)
            arguments
                obj
                timeStep (1, 1) double {mustBeNonnegative} % seconds
            end
            if obj.rpm
                switch obj.scanDirection
                    case "cw"
                        obj.yaw = mod(obj.yaw - timeStep*obj.rpm*6 + 180, 360) - 180;
                    case "ccw"
                        obj.yaw = mod(obj.yaw + timeStep*obj.rpm*6 + 180, 360) - 180;
                end
            else
                switch obj.scanDirection
                    case "cw"
                        obj.steeringAzimuth = obj.steeringAzimuth - obj.steeringAzimuthStep;
                        if obj.steeringAzimuth < -obj.steeringAzimuthLimit
                            obj.steeringAzimuth = obj.steeringAzimuthLimit;
                        end
                    case "ccw"
                        obj.steeringAzimuth = obj.steeringAzimuth + obj.steeringAzimuthStep;
                        if obj.steeringAzimuth > obj.steeringAzimuthLimit
                            obj.steeringAzimuth = -obj.steeringAzimuthLimit;
                        end
                end
            end
            obj.durationSimulation = obj.durationSimulation + timeStep;
        end

        function reset(obj)
            if obj.rpm
                switch obj.scanDirection
                    case "cw"
                        obj.yaw = mod(obj.yaw + obj.durationSimulation*obj.rpm*6 + 180, 360) - 180;
                    case "ccw"
                        obj.yaw = mod(obj.yaw - obj.durationSimulation*obj.rpm*6 + 180, 360) - 180;
                end
            else
                obj.steeringAzimuth = obj.steeringAzimuthInitial;
            end
            obj.durationSimulation = 0;
        end

        %%% visualization methods

        function visualizearray(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            plot3(obj.positions(1, :), obj.positions(2, :), obj.positions(3, :), 'o');
            xlabel('x-axis (m)'); ylabel('y-axis (m)'); zlabel('z-axis (m)');
            grid off; grid on; grid minor;
            % n = obj.normalVector.*max(obj.positions, [], 'all')/2;
            % hold on; quiver3(0, 0, 0, n(1), n(2), n(3));
            % title('sensor configuration'); hold off;
            drawnow;
        end
    end
end