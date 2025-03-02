classdef target < handle
    %target Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = private, GetAccess = public)
        position (3, :) double = zeros(3, 1) % meters
    end

    properties (SetAccess = private, GetAccess = public)
        velocity (3, :) double = zeros(3, 1) % meters per second
        meanRCS_dbsm (1, :) double = 0 % dB of meter squares
    end

    properties (Dependent)
        numberOfTargets (1, 1) double {mustBeInteger, mustBeNonnegative}
        velocityMagnitude (1, :) double {mustBeNonnegative} % meters per second
        velocityUnitDirection (3, :) double
        RCS_dbms (1, :) double {mustBePositive} % meter squares
    end

    properties (Access = private)
        x (1, :) double = nan
        y (1, :) double = nan
        z (1, :) double = nan
        durationSimulation (1, 1) {mustBeNonnegative} = 0
    end

    properties (Dependent, GetAccess = ?interface, SetAccess = private)
        firstGridLength (1, 1) double {mustBeInteger, mustBeNonnegative}
        secondGridLength (1, 1) double {mustBeInteger, mustBeNonnegative}
        firstAxisID (1, 1) double {mustBeInteger, mustBePositive, mustBeMember(firstAxisID, [1, 2, 3])}
        secondAxisID (1, 1) double {mustBeInteger, mustBePositive, mustBeMember(secondAxisID, [1, 2, 3])}
    end

    methods
        function obj = target(options)
            arguments
                options.position (3, :) double = zeros(3, 1)
                options.velocity (3, :) double = zeros(3, 1)
                options.meanRCS_dbsm (1, :) double = 0
                options.x (1, :) double = nan
                options.y (1, :) double = nan
                options.z (1, :) double = nan
            end
            %target Construct an instance of this class
            %   Detailed explanation goes here
            obj.x = options.x;
            obj.y = options.y;
            obj.z = options.z;
            if ~isscalar(options.x)
                if ~isscalar(options.y)
                    [options.x, options.y] = meshgrid(options.x, options.y);
                    options.position = [options.x(:).'; options.y(:).'; zeros(1, numel(options.x))];
                elseif ~isscalar(options.z)
                    [options.z, options.x] = meshgrid(options.z, options.x);
                    options.position = [options.x(:).'; zeros(1, numel(options.x)); options.z(:).'];
                else
                    options.position = [options.x(:).'; zeros(1, numel(options.x)); zeros(1, numel(options.x))];
                end
            elseif ~isscalar(options.y)
                if ~isscalar(options.z)
                    [options.y, options.z] = meshgrid(options.y, options.z);
                    options.position = [zeros(1, numel(options.y)); options.y(:).'; options.z(:).'];
                else
                    options.position = [zeros(1, numel(options.y)); options.y(:).'; zeros(1, numel(options.y))];
                end
            elseif ~isscalar(options.z)
                options.position = [zeros(1, numel(options.z)); zeros(1, numel(options.z)); options.z(:).';];
            end
            numberOfTargets = max([size(options.position, 2), size(options.velocity, 2), size(options.meanRCS_dbsm, 2)]);
            if size(options.position, 2) == 1
                options.position = repmat(options.position, 1, numberOfTargets);
            elseif size(options.position, 2) ~= numberOfTargets
                error('number of targets is %d', numberOfTargets);
            end
            if size(options.velocity, 2) == 1
                options.velocity = repmat(options.velocity, 1, numberOfTargets);
            elseif size(options.velocity, 2) ~= numberOfTargets
                error('number of targets is %d', numberOfTargets);
            end
            if size(options.meanRCS_dbsm, 2) == 1
                options.meanRCS_dbsm = repmat(options.meanRCS_dbsm, 1, numberOfTargets);
            elseif size(options.meanRCS_dbsm, 2) ~= numberOfTargets
                error('number of targets is %d', numberOfTargets);
            end
            obj.position = options.position;
            obj.velocity = options.velocity;
            obj.meanRCS_dbsm = options.meanRCS_dbsm;
        end

        function N = get.numberOfTargets(obj)
            N = size(obj.position, 2);
        end

        function v = get.velocityMagnitude(obj)
            v = sqrt(sum(abs(obj.velocity).^2));
        end

        function u = get.velocityUnitDirection(obj)
            u = zeros(3, obj.numberOfTargets);
            nonzeroVelocities = logical(obj.velocityMagnitude);
            u(:, nonzeroVelocities) = obj.velocity(:, nonzeroVelocities)./obj.velocityMagnitude(1, nonzeroVelocities);
        end

        function rcs = get.RCS_dbms(obj)
            rcs = obj.meanRCS_dbsm;
        end

        function step(obj, timeStep)
            arguments
                obj
                timeStep (1, 1) double {mustBeNonnegative} % seconds
            end
            obj.position = obj.position + timeStep*obj.velocity;
            obj.durationSimulation = obj.durationSimulation + timeStep;
        end

        function reset(obj)
            obj.position = obj.position - obj.durationSimulation*obj.velocity;
            obj.durationSimulation = 0;
        end

        function visualizetargets(obj)
            x0 = obj.position(1, :)/1e3;
            y0 = obj.position(2, :)/1e3;
            z0 = obj.position(3, :)/1e3;
            vx = obj.velocity(1, :)/1e3;
            vy = obj.velocity(2, :)/1e3;
            vz = obj.velocity(3, :)/1e3;
            figure; plot3(x0, y0, z0, '+k');
            hold on; quiver3(x0, y0, z0, vx, vy, vz);
            grid on; grid minor; title('target scenario');
            xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
        end

        function L = get.firstGridLength(obj)
            if ~isscalar(obj.x) && ~isscalar(obj.y)
                L = numel(obj.y);
            elseif ~isscalar(obj.x) && ~isscalar(obj.z)
                L = numel(obj.x);
            elseif ~isscalar(obj.y) && ~isscalar(obj.z)
                L = numel(obj.z);
            else
                L = 1;
            end
        end

        function L = get.secondGridLength(obj)
            if ~isscalar(obj.x)
                if ~isscalar(obj.y)
                    L = numel(obj.x);
                elseif ~isnan(obj.z)
                    L = numel(obj.z);
                else
                    L = numel(obj.x);
                end
            elseif ~isscalar(obj.y)
                L = numel(obj.y);
            else
                L = numel(obj.z);
            end
        end

        function axisID = get.firstAxisID(obj)
            if ~isscalar(obj.x)
                if ~isscalar(obj.y)
                    axisID = 1;
                elseif ~isnan(obj.z)
                    axisID = 3;
                else
                    axisID = 1;
                end
            elseif ~isscalar(obj.y)
                axisID = 2;
            else
                axisID = 3;
            end
        end

        function axisID = get.secondAxisID(obj)
            if ~isscalar(obj.x) && ~isscalar(obj.y)
                axisID = 2;
            elseif ~isscalar(obj.x) && ~isscalar(obj.z)
                axisID = 1;
            elseif ~isscalar(obj.y) && ~isscalar(obj.z)
                axisID = 3;
            else
                axisID = 1;
            end
        end
    end
end