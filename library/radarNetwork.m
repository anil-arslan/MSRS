classdef radarNetwork < handle
    %radarNetwork Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = private, GetAccess = public)
        receivingNodes (1, :) receivingNode = receivingNode.empty() % (1 x Nrx vector)
        transmittingNodes (1, :) transmittingNode = transmittingNode.empty() % (1 x Ntx vector)
        receivingNodeActivity (1, :) logical = false % (1 x Nrx vector)
        transmittingNodeActivity (1, :) logical = false % (1 x Ntx vector)
        networkCoherency (1, 1) string {mustBeMember(networkCoherency, ["coherent", "short-term coherent", "incoherent"])} = "coherent"
        carrierMode (1, 1) string {mustBeMember(carrierMode, "bandPassProcessing")} = "bandPassProcessing"
            % bandPassProcessing: transmitted signals perfectly resolved in carrier frequency
        networkMode (1, 1) string {mustBeMember(networkMode, ["multiStatic", "monoStatic"])} = "multiStatic"
        surveillanceMode (1, 1) string {mustBeMember(surveillanceMode, ["rotating", "electronicScan", "staticBeam", "custom"])} = "custom"
        beamTime (1, 1) double {mustBeNonnegative} % sec
        timeOffsetBistaticPairs double = 0 % (Ntx x Nrx matrix)
        timeDeviationBistaticPairs double {mustBeNonnegative} = 0 % (Ntx x Nrx matrix)
        frequencyOffsetBistaticPairs double = 0 % (Ntx x 1 matrix)
        frequencyDeviationBistaticPairs double {mustBeNonnegative} = 0 % (Ntx x 1 matrix)
    end

    properties (Dependent)
        timeOffsetActiveBistaticPairs double % (Ntx x Nrx matrix)
        timeDeviationActiveBistaticPairs double {mustBeNonnegative} % (Ntx x Nrx matrix)
        frequencyOffsetActiveBistaticPairs double % (Ntx x Nrx matrix)
        frequencyDeviationActiveBistaticPairs double {mustBeNonnegative} % (Ntx x Nrx matrix)
        numberOfNodes (1, 1) double {mustBeInteger, mustBeNonnegative}
        numberOfActiveNodes (1, 1) double {mustBeInteger, mustBeNonnegative}
        numberOfReceivingNodes (1, 1) double {mustBeInteger, mustBeNonnegative}
        numberOfTransmittingNodes (1, 1) double {mustBeInteger, mustBeNonnegative}
        numberOfBistaticPairs (1, 1) double {mustBeInteger, mustBeNonnegative}
        numberOfActiveReceivingNodes (1, 1) double {mustBeInteger, mustBeNonnegative}
        numberOfActiveTransmittingNodes (1, 1) double {mustBeInteger, mustBeNonnegative}
        numberOfActiveBistaticPairs (1, 1) double {mustBeInteger, mustBeNonnegative}
        multiReceiver (1, 1) logical
        multiTransmitter (1, 1) logical
        passiveNetwork (1, 1) logical
        numberOfParallelSignalFusions (1, 1) double {mustBeInteger, mustBeNonnegative}
        monoStaticTransmitterIDs (1, :) double {mustBeInteger, mustBeNonnegative} % (1 x Nrx vector)
        activeReceivingNodes (1, :) receivingNode % (1 x Nrx vector)
        activeTransmittingNodes (1, :) transmittingNode % (1 x Ntx vector)
        positionReceivingNode (3, :, :) double % (3 x Ntx x Nrx matrix)
        positionTransmittingNode (3, :, :) double % (3 x Ntx x Nrx matrix)
        distanceBaseline double {mustBePositive} % (Ntx x Nrx matrix)
        directPathDelay double {mustBePositive} % (Ntx x Nrx matrix)
        center (3, 1) double
        boundary (3, 2) double
        boundaryListened (3, 2) double
        defaultGridResolution (3, 1) double {mustBeNonnegative}
        pulseWidthSample double % (Nrx x Ntx matrix)
        dutyCycles double % (Nrx x Nrx matrix)
        matchFilter double % (L x Nrx x Ntx matrix)
        matchFilterGain_dB double % (Nrx x Ntx matrix)
        beamformingGain_dB double % (Nrx x 1 matrix)
        processingGain_dB double % (Nrx x Ntx matrix)
        noisePSDmatrix double % W/Hz (Nrx x Nrx matrix)
    end

    properties (Constant)
        speedOfLight (1, 1) double {mustBePositive} = physconst('LightSpeed');
    end

    properties (Access = private)
        arrayScanParameters struct = struct( ...
            "rpmRX", [], ...
            "rpmTX", [], ...
            "steeringAzimuthStepRX", [], ...
            "steeringAzimuthStepTX", [] ...
            )
    end

    methods
        function obj = radarNetwork(options)
            arguments
                options.receivingNodes (1, :) receivingNode = receivingNode.empty()
                options.transmittingNodes (1, :) transmittingNode = transmittingNode.empty()
                options.networkCoherency (1, 1) string {mustBeMember(options.networkCoherency, ["coherent", "short-term coherent", "incoherent"])} = "coherent"
            end
            %radarNetwork Construct an instance of this class
            %   Detailed explanation goes here
            obj.receivingNodes = options.receivingNodes;
            obj.transmittingNodes = options.transmittingNodes;
            obj.transmittingNodeActivity = true(1, obj.numberOfTransmittingNodes);
            obj.receivingNodeActivity = true(1, obj.numberOfReceivingNodes);
            obj.networkCoherency = options.networkCoherency;
            for rxID = 1 : obj.numberOfReceivingNodes
                if ~all(obj.receivingNodes(rxID).array.spacingSpecified)
                    obj.receivingNodes(rxID).array.spacing(~obj.receivingNodes(rxID).array.spacingSpecified) = obj.speedOfLight./obj.transmittingNodes(1).carrierFrequency/2;
                    N = obj.receivingNodes(rxID).array.numberOfElements;
                    yPos = (-(N(1) - 1)/2 : (N(1) - 1)/2)*obj.receivingNodes(rxID).array.spacing(1);
                    zPos = (-(N(2) - 1)/2 : (N(2) - 1)/2)*obj.receivingNodes(rxID).array.spacing(2);
                    [yPos, zPos] = meshgrid(yPos, zPos);
                    obj.receivingNodes(rxID).array.positionsYZ = zeros(2, prod(N));
                    obj.receivingNodes(rxID).array.positionsYZ(1, :) = yPos(:).';
                    obj.receivingNodes(rxID).array.positionsYZ(2, :) = zPos(:).';
                    obj.receivingNodes(rxID).array.spacingSpecified = true(2, 1);
                end
            end
            for txID = 1 : obj.numberOfTransmittingNodes
                if ~all(obj.transmittingNodes(txID).array.spacingSpecified)
                    obj.transmittingNodes(txID).array.spacing(~obj.transmittingNodes(txID).array.spacingSpecified) = obj.speedOfLight./obj.transmittingNodes(txID).carrierFrequency/2;
                    N = obj.transmittingNodes(txID).array.numberOfElements;
                    yPos = (-(N(1) - 1)/2 : (N(1) - 1)/2)*obj.transmittingNodes(txID).array.spacing(1);
                    zPos = (-(N(2) - 1)/2 : (N(2) - 1)/2)*obj.transmittingNodes(txID).array.spacing(2);
                    [yPos, zPos] = meshgrid(yPos, zPos);
                    obj.transmittingNodes(txID).array.positionsYZ = zeros(2, prod(N));
                    obj.transmittingNodes(txID).array.positionsYZ(1, :) = yPos(:).';
                    obj.transmittingNodes(txID).array.positionsYZ(2, :) = zPos(:).';
                    obj.transmittingNodes(txID).array.spacingSpecified = true(2, 1);
                end
            end
            obj.beamTime = max([obj.receivingNodes.CPI]);
            obj.timeOffsetBistaticPairs = zeros(obj.numberOfTransmittingNodes, obj.numberOfActiveReceivingNodes);
            obj.timeDeviationBistaticPairs = zeros(obj.numberOfTransmittingNodes, obj.numberOfActiveReceivingNodes);
            obj.frequencyOffsetBistaticPairs = zeros(obj.numberOfTransmittingNodes, 1);
            obj.frequencyDeviationBistaticPairs = zeros(obj.numberOfTransmittingNodes, 1);
        end

        %%% get methods

        function t = get.timeOffsetActiveBistaticPairs(obj)
            t = obj.timeOffsetBistaticPairs(obj.transmittingNodeActivity, obj.receivingNodeActivity);
        end

        function t = get.timeDeviationActiveBistaticPairs(obj)
            t = obj.timeDeviationBistaticPairs(obj.transmittingNodeActivity, obj.receivingNodeActivity);
        end

        function t = get.frequencyOffsetActiveBistaticPairs(obj)
            t = obj.frequencyOffsetBistaticPairs(obj.transmittingNodeActivity, 1);
        end

        function t = get.frequencyDeviationActiveBistaticPairs(obj)
            t = obj.frequencyDeviationBistaticPairs(obj.transmittingNodeActivity, 1);
        end

        function N = get.numberOfNodes(obj)
            N = obj.numberOfReceivingNodes + obj.numberOfTransmittingNodes;
        end

        function N = get.numberOfActiveNodes(obj)
            N = obj.numberOfActiveReceivingNodes + obj.numberOfActiveTransmittingNodes;
        end

        function N = get.numberOfReceivingNodes(obj)
            N = numel(obj.receivingNodes);
        end

        function N = get.numberOfTransmittingNodes(obj)
            N = numel(obj.transmittingNodes);
        end

        function N = get.numberOfBistaticPairs(obj)
            N = obj.numberOfReceivingNodes*obj.numberOfTransmittingNodes;
        end

        function N = get.numberOfActiveReceivingNodes(obj)
            N = nnz(obj.receivingNodeActivity);
        end

        function N = get.numberOfActiveTransmittingNodes(obj)
            N = nnz(obj.transmittingNodeActivity);
        end

        function N = get.numberOfActiveBistaticPairs(obj)
            switch obj.networkMode
                case 'multiStatic'
                    N = obj.numberOfActiveReceivingNodes*obj.numberOfActiveTransmittingNodes;
                case 'monoStatic'
                    N = 1;
            end
        end

        function tf = get.multiReceiver(obj)
            tf = logical(obj.numberOfActiveReceivingNodes - 1);
        end

        function tf = get.multiTransmitter(obj)
            tf = logical(obj.numberOfActiveTransmittingNodes - 1);
        end

        function tf = get.passiveNetwork(obj)
            tf = ~logical(obj.numberOfActiveTransmittingNodes);
        end

        function nfus = get.numberOfParallelSignalFusions(obj)
            switch obj.networkMode
                case 'multiStatic'
                    nfus = 1;
                case 'monoStatic'
                    nfus = nnz(~isnan(obj.monoStaticTransmitterIDs));
            end
        end

        function txID = get.monoStaticTransmitterIDs(obj)
            txID = nan(1, obj.numberOfActiveReceivingNodes);
            for rxID = 1 : obj.numberOfActiveReceivingNodes
                id = find(~sum(abs(obj.activeReceivingNodes(rxID).position - [obj.activeTransmittingNodes.position])));
                if ~isempty(id)
                    txID(rxID) = id;
                end
            end
        end

        function rx = get.activeReceivingNodes(obj)
            rx = obj.receivingNodes(obj.receivingNodeActivity);
        end

        function tx = get.activeTransmittingNodes(obj)
            tx = obj.transmittingNodes(obj.transmittingNodeActivity);
        end

        function pos = get.positionReceivingNode(obj)
            pos = repelem([obj.activeReceivingNodes.position], 1, obj.numberOfActiveTransmittingNodes);
            pos = reshape(pos, 3, obj.numberOfActiveTransmittingNodes, obj.numberOfActiveReceivingNodes);
        end

        function pos = get.positionTransmittingNode(obj)
            pos = repmat([obj.activeTransmittingNodes.position], 1, obj.numberOfActiveReceivingNodes);
            pos = reshape(pos, 3, obj.numberOfActiveTransmittingNodes, obj.numberOfActiveReceivingNodes);
        end

        function D = get.distanceBaseline(obj)
            D = sqrt(sum((obj.positionReceivingNode - obj.positionTransmittingNode).^2));
            D = reshape(D, obj.numberOfActiveTransmittingNodes, obj.numberOfActiveReceivingNodes);
        end

        function tau = get.directPathDelay(obj)
            tau = obj.distanceBaseline./obj.speedOfLight;
        end

        function origin = get.center(obj)
            switch obj.networkMode
                case 'monoStatic'
                    txIDs = obj.monoStaticTransmitterIDs;
                case 'multiStatic'
                    txIDs = 1 : obj.numberOfActiveTransmittingNodes;
            end
            if isempty(txIDs)
                origin = zeros(3, 1);
            else
                origin = mean([obj.activeTransmittingNodes(txIDs).position, obj.activeReceivingNodes.position], 2);
            end
        end

        function lims = get.boundary(obj)
            switch obj.networkMode
                case 'monoStatic'
                    txIDs = obj.monoStaticTransmitterIDs;
                case 'multiStatic'
                    txIDs = 1 : obj.numberOfActiveTransmittingNodes;
            end
            lims = zeros(3, 2);
            if ~isempty(txIDs)
                nodePositions = [obj.activeTransmittingNodes(txIDs).position, obj.activeReceivingNodes.position];
                lims(:, 1) = min(nodePositions, [], 2);
                lims(:, 2) = max(nodePositions, [], 2);
            end
        end

        function lims = get.boundaryListened(obj)
            lims = obj.boundary;
            offset = max([obj.activeReceivingNodes.listenedRadius]);
            if ~isempty(offset)
                lims(:, 1) = lims(:, 1) - offset;
                lims(:, 2) = lims(:, 2) + offset;
            end
        end

        function res = get.defaultGridResolution(obj)
            res = diff(obj.boundaryListened, [], 2)/100;
        end

        function N = get.pulseWidthSample(obj)
            % (Nrx x Ntx matrix)
            N = [obj.activeTransmittingNodes.pulseWidth].*[obj.activeReceivingNodes.samplingFrequency].';
        end

        function D = get.dutyCycles(obj)
            D = obj.pulseWidthSample./[obj.activeReceivingNodes.numberOfSamplesPerCPI].';
        end

        function mf = get.matchFilter(obj)
            % (L x Nrx x Ntx matrix)
            Ts = [obj.activeReceivingNodes.samplingPeriod]; % 1 x Nrx matrix
            if isempty(Ts)
                mf = 1;
                return;
            end
            N = obj.pulseWidthSample; % Nrx x Ntx matrix
            Nunique = unique(N);
            switch length(Nunique)
                case 1
                    switch obj.activeTransmittingNodes(1).transmissionType
                        case "continuous" % infinite pulse width
                            t = [obj.receivingNodes.samplingInstants];
                            mf = ones(size(t, 1), obj.numberOfActiveReceivingNodes, obj.numberOfActiveTransmittingNodes);
                        case "pulsed"
                            t = zeros(Nunique, obj.numberOfActiveReceivingNodes);
                            mf = zeros(Nunique, obj.numberOfActiveReceivingNodes, obj.numberOfActiveTransmittingNodes);
                            for txID = 1 : obj.numberOfActiveTransmittingNodes
                                for rxID = 1 : obj.numberOfActiveReceivingNodes
                                    t(:, rxID) = ((1 - N(rxID, txID) : 0)*Ts(rxID));
                                end
                                mf(:, :, txID) = obj.activeTransmittingNodes(txID).waveform(t);
                            end
                    end
                    mf = mf./sqrt(sum(abs(mf).^2, 1));
                otherwise
                    error('not implemented');
            end
        end

        function G = get.matchFilterGain_dB(obj)
            G = permute(20*log10(sum(abs(obj.matchFilter), 1)), [2 3 1]);
        end

        function G = get.beamformingGain_dB(obj)
            G = [obj.activeReceivingNodes.beamformingGain_dB].';
        end

        function G = get.processingGain_dB(obj)
            G = obj.matchFilterGain_dB + obj.beamformingGain_dB;
        end

        function B = get.noisePSDmatrix(obj)
            B = diag(10.^(.1*[obj.activeReceivingNodes.noisePowerPerSample_dB])./[obj.activeReceivingNodes.samplingFrequency]);
        end

        %%% utility methods

        function int = getinterface(obj, targets)
            int = interface("network", obj, "targets", targets);
        end

        function targets = randomtargets(obj, options)
            arguments
                obj
                options.numberOfTargets (1, 1) {mustBeInteger, mustBePositive} = 1
                options.meanRCS_dbms (1, 1) double = 0
            end
            dimensions = [1 2];
            boundaries = obj.boundaryListened/2;
            targetPositions = zeros(3, options.numberOfTargets);
            for dimID = 1 : length(dimensions)
                if boundaries(dimensions(dimID), 1)
                    targetPositions(dimensions(dimID), :) = boundaries(dimensions(dimID), 1).*((boundaries(dimensions(dimID), 2)./boundaries(dimensions(dimID), 1) - 1).*rand(1, options.numberOfTargets) + 1);
                else
                    targetPositions(dimensions(dimID), :) = boundaries(dimensions(dimID), 2).*rand(1, options.numberOfTargets);
                end
            end
            targets = target("position", targetPositions, "meanRCS_dbms", options.meanRCS_dbms);
        end

        %%% set methods

        function setsynchronization(obj, options)
            arguments
                obj
                options.timeOffsetBistaticPairs double = obj.timeOffsetBistaticPairs
                options.timeDeviationBistaticPairs double {mustBeNonnegative} = obj.timeDeviationBistaticPairs
                options.frequencyOffsetBistaticPairs (:, 1) double = obj.frequencyOffsetBistaticPairs
                options.frequencyDeviationBistaticPairs (:, 1) double {mustBeNonnegative} = obj.frequencyDeviationBistaticPairs
            end
            obj.timeOffsetBistaticPairs = sizeupdate(options.timeOffsetBistaticPairs);
            obj.timeDeviationBistaticPairs = sizeupdate(options.timeDeviationBistaticPairs);
            options.frequencyOffsetBistaticPairs = sizeupdate(options.frequencyOffsetBistaticPairs);
            options.frequencyDeviationBistaticPairs = sizeupdate(options.frequencyDeviationBistaticPairs);
            obj.frequencyOffsetBistaticPairs = options.frequencyOffsetBistaticPairs(:, 1);
            obj.frequencyDeviationBistaticPairs = options.frequencyDeviationBistaticPairs(:, 1);
            function output = sizeupdate(input)
                if isscalar(input)
                    output = repmat(input, obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
                elseif size(input, 1) == 1
                    if size(input, 2) == obj.numberOfReceivingNodes
                        output = repmat(input, obj.numberOfTransmittingNodes, 1);
                    elseif size(input, 2) == obj.numberOfTransmittingNodes
                        output = repmat(input.', 1, obj.numberOfReceivingNodes);
                    else
                        error('size of the input must be [%d x %d]', obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
                    end
                elseif size(input, 2) == 1
                    if size(input, 1) == obj.numberOfReceivingNodes
                        output = repmat(input.', obj.numberOfTransmittingNodes, 1);
                    elseif size(input, 1) == obj.numberOfTransmittingNodes
                        output = repmat(input, 1, obj.numberOfReceivingNodes);
                    else
                        error('size of the input must be [%d x %d]', obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
                    end
                else
                    if ~isequal(size(input), [obj.numberOfTransmittingNodes obj.numberOfReceivingNodes])
                        error('size of the input must be [%d x %d]', obj.numberOfTransmittingNodes, obj.numberOfReceivingNodes);
                    else
                        output = input;
                    end
                end
            end
        end

        function setsurveillance(obj, surveillanceMode)
            arguments
                obj
                surveillanceMode (1, 1) string {mustBeMember(surveillanceMode, ["rotating", "electronicScan", "staticBeam", "custom"])} = obj.surveillanceMode
            end
            arraysRX = [obj.receivingNodes.array];
            arraysTX = [obj.transmittingNodes.array];
            if ~strcmpi(obj.surveillanceMode, "staticBeam")
                if any(strcmpi(obj.surveillanceMode, ["custom", "rotating"]))
                    obj.arrayScanParameters.rpmRX = [arraysRX.rpm];
                    obj.arrayScanParameters.rpmTX = [arraysTX.rpm];
                end
                if any(strcmpi(obj.surveillanceMode, ["custom", "electronicScan"]))
                    obj.arrayScanParameters.steeringAzimuthStepRX = [arraysRX.steeringAzimuthStep];
                    obj.arrayScanParameters.steeringAzimuthStepTX = [arraysTX.steeringAzimuthStep];
                end
            end
            switch surveillanceMode
                case "rotating"
                    arraysRX.setscanparameters("steeringAzimuthStep", 0, "rpm", obj.arrayScanParameters.rpmRX);
                    arraysTX.setscanparameters("steeringAzimuthStep", 0, "rpm", obj.arrayScanParameters.rpmTX);
                case "electronicScan"
                    arraysRX.setscanparameters("steeringAzimuthStep", obj.arrayScanParameters.steeringAzimuthStepRX, "rpm", 0);
                    arraysTX.setscanparameters("steeringAzimuthStep", obj.arrayScanParameters.steeringAzimuthStepTX, "rpm", 0);
                case "staticBeam"
                    arraysRX.setscanparameters("steeringAzimuthStep", 0, "rpm", 0);
                    arraysTX.setscanparameters("steeringAzimuthStep", 0, "rpm", 0);
                case "custom"
                    arraysRX.setscanparameters("steeringAzimuthStep", obj.arrayScanParameters.steeringAzimuthStepRX, "rpm", obj.arrayScanParameters.rpmRX);
                    arraysTX.setscanparameters("steeringAzimuthStep", obj.arrayScanParameters.steeringAzimuthStepTX, "rpm", obj.arrayScanParameters.rpmTX);
            end
            obj.surveillanceMode = surveillanceMode;
        end

        function activitycontrol(obj, options)
            arguments
                obj
                options.transmittingNodeActivity (1, :) logical = obj.transmittingNodeActivity
                options.receivingNodeActivity (1, :) logical = obj.receivingNodeActivity
            end
            obj.transmittingNodeActivity = options.transmittingNodeActivity;
            obj.receivingNodeActivity = options.receivingNodeActivity;
            if strcmpi(obj.networkMode, "monoStatic")
                obj.receivingNodeActivity(obj.receivingNodeActivity) = ~isnan(obj.monoStaticTransmitterIDs);
            end
        end

        function settingsnetwork(obj, options)
            arguments
                obj
                options.networkCoherency (1, 1) string {mustBeMember(options.networkCoherency, ["coherent", "short-term coherent", "incoherent"])} = obj.networkCoherency
                options.carrierMode (1, 1) string {mustBeMember(options.carrierMode, "bandPassProcessing")} = obj.carrierMode
                options.networkMode (1, 1) string {mustBeMember(options.networkMode, ["multiStatic", "monoStatic"])} = obj.networkMode
                options.beamTime (1, 1) double {mustBeNonnegative} = obj.beamTime % sec
            end
            obj.networkCoherency = options.networkCoherency;
            obj.carrierMode = options.carrierMode;
            obj.networkMode = options.networkMode;
            if strcmpi(obj.networkMode, "monoStatic")
                obj.receivingNodeActivity = obj.receivingNodeActivity & ~isnan(obj.monoStaticTransmitterIDs);
            end
            obj.beamTime = options.beamTime;
        end

        function setdirectionfinders(obj)
            for rxID = 1 : obj.numberOfActiveReceivingNodes
                lambda = obj.activeTransmittingNodes(obj.monoStaticTransmitterIDs(rxID)).carrierWavelength;
                azimuthBeamWidth = 70*lambda/(obj.activeReceivingNodes(rxID).array.numberOfElements(1, :).*obj.activeReceivingNodes(rxID).array.spacing(1, :));
                elevationBeamWidth = 70*lambda/(obj.activeReceivingNodes(rxID).array.numberOfElements(2, :).*obj.activeReceivingNodes(rxID).array.spacing(2, :));
                azimuthRes = azimuthBeamWidth/10;
                elevationRes = elevationBeamWidth/10;
                switch obj.activeReceivingNodes(rxID).beamformingMode
                    case 'bypass'
                        azimuthBoundary = [-90 90];
                        elevationBoundary = [-90 90];
                        azimuthDF = azimuthBoundary(1) : azimuthRes : azimuthBoundary(2);
                        elevationDF = elevationBoundary(1) : elevationRes : elevationBoundary(2);
                        [azimuthDF, elevationDF] = meshgrid(azimuthDF, elevationDF);
                        unitDirectionTargets = [cosd(elevationDF(:)).*sind(azimuthDF(:)), sind(elevationDF(:))];
                        tableDF = exp(-1j*2*pi*unitDirectionTargets*obj.activeReceivingNodes(rxID).array.positionsYZ./lambda);
                        numberOfTargets = size(unitDirectionTargets, 1);
                        for targetID = 1 : numberOfTargets
                            tableDF(targetID, :) = tableDF(targetID, :)./norm(tableDF(targetID, :));
                        end
                    otherwise
                        azimuthBoundary = [min(obj.activeReceivingNodes(rxID).beamCentersAzimuth), max(obj.activeReceivingNodes(rxID).beamCentersAzimuth)];
                        elevationBoundary = [min(obj.activeReceivingNodes(rxID).beamCentersElevation), max(obj.activeReceivingNodes(rxID).beamCentersElevation)];
                        if diff(azimuthBoundary) > azimuthRes
                            azimuthDF = azimuthBoundary(1) - azimuthBeamWidth : azimuthRes : azimuthBoundary(2) + azimuthBeamWidth;
                        else
                            azimuthDF = 0;
                            azimuthRes = 0;
                        end
                        if diff(elevationBoundary) > elevationRes
                            elevationDF = elevationBoundary(1) - elevationBeamWidth : elevationRes : elevationBoundary(2) + elevationBeamWidth;
                        else
                            elevationDF = 0;
                            elevationRes = 0;
                        end
                        [azimuthDF, elevationDF] = meshgrid(azimuthDF, elevationDF);
                        unitDirectionTargets = [cosd(elevationDF(:)).*sind(azimuthDF(:)), sind(elevationDF(:))];
                        a = exp(1j*2*pi*unitDirectionTargets*obj.activeReceivingNodes(rxID).array.positionsYZ./lambda);
                        H = obj.activeReceivingNodes(rxID).beamformer.*exp(1j*2*pi*obj.activeReceivingNodes(rxID).beamCenterPositions./lambda);
                        tableDF = conj(a)*H;
                        numberOfTargets = size(unitDirectionTargets, 1);
                        for targetID = 1 : numberOfTargets
                            tableDF(targetID, :) = tableDF(targetID, :)./norm(tableDF(targetID, :));
                        end
                end
                obj.receivingNodes(rxID).directionFinder.table = tableDF;
                obj.receivingNodes(rxID).directionFinder.elevation = elevationDF;
                obj.receivingNodes(rxID).directionFinder.azimuth = azimuthDF;
                obj.receivingNodes(rxID).directionFinder.azimuthResolution = azimuthRes;
                obj.receivingNodes(rxID).directionFinder.elevationResolution = elevationRes;
            end
        end

        %%% visualizeton methods

        function visualizewaveformsampled(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
                options.domain (1, 1) string {mustBeMember(options.domain, ["ambiguity", "frequency", "time"])} = "ambiguity"
                options.axisAmbiguity (1, 1) string {mustBeMember(options.axisAmbiguity, ["3D", "zeroDoppler", "zeroRange"])} = "zeroDoppler"
                options.plot (1, 1) string {mustBeMember(options.plot, ["magnitude", "phase", "real", "imaginary"])} = "magnitude"
            end
            mfAll = obj.matchFilter;
            if size(mfAll, 1) == 1
                warning('match filter is single sample');
                return;
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            minMag = -30; % dB
            N = unique(obj.pulseWidthSample);
            nfft = 2.^(nextpow2(N) + 1);
            bins = -nfft/2 : nfft/2 - 1;
            delays = -N + 1 : N - 1;
            if ~strcmpi(options.domain, "ambiguity")
                options.axisAmbiguity = "";
            end
            if strcmpi(options.axisAmbiguity, "3D")
                txIDs = obj.numberOfActiveTransmittingNodes;
            else
                txIDs = 1 : obj.numberOfActiveTransmittingNodes;
            end
            for txID = txIDs
                for rxID = 1 : obj.numberOfActiveReceivingNodes
                    if strcmpi(options.axisAmbiguity, "3D") && rxID ~= 1
                        break;
                    end
                    mf = mfAll(:, rxID, txID);
                    switch options.domain
                        case "ambiguity"
                            options.plot = "magnitude";
                            switch options.axisAmbiguity
                                case "zeroDoppler"
                                    af = xcorr(mf);
                                    plot(delays, 20*log10(abs(af)));
                                case "zeroRange"
                                    af = fftshift(fft(abs(mf).^2, nfft));
                                    plot(bins, 20*log10(abs(af)));
                                case "3D"
                                    af = zeros(2*N - 1, nfft);
                                    for k = 1 : nfft(1)
                                        af(:, k) = xcorr(mf.*exp(1j*2*pi*bins(k)*(0 : N - 1).'/nfft), mf);
                                    end
                                    af = 20*log10(abs(af));
                                    af(af < minMag) = minMag;
                                    m = mesh(bins, delays, af);
                                    colorbar;
                            end
                        case "time"
                            switch options.plot
                                case "magnitude"
                                    plot(0 : N - 1, abs(mf));
                                case "phase"
                                    plot(0 : N - 1, 180*unwrap(angle(mf))/pi);
                                case "real"
                                    plot(0 : N - 1, real(mf));
                                case "imaginary"
                                    plot(0 : N - 1, imag(mf));
                            end
                        case "frequency"
                            switch options.plot
                                case "magnitude"
                                    plot(bins, 20*log10(fftshift(abs(fft(mf, nfft)))));
                                case "phase"
                                    plot(bins, 180*unwrap(angle(fftshift(fft(mf, nfft))))/pi);
                                case "real"
                                    plot(bins, 20*log10(fftshift(abs(fft(mf, nfft)))));
                                case "imaginary"
                                    plot(bins, 20*log10(fftshift(abs(fft(mf, nfft)))));
                            end
                    end
                    hold on;
                end
            end
            switch options.domain
                case "ambiguity"
                    title('ambiguity functions');
                case "frequency"
                    title('waveforms in time');
                    xlabel('frequency bins');
                case "time"
                    title('frequency responses of waveform');
                    xlabel('time samples');
            end
            grid off; grid on; grid minor;
            if ~strcmpi(options.axisAmbiguity, "3D")
                if strcmpi(options.axisAmbiguity, "zeroDoppler")
                    xlabel('delay samples');
                elseif strcmpi(options.domain, "ambiguity")
                    xlabel('doppler bins');
                end
                switch options.plot
                    case "magnitude"
                        ylabel('power (dB)');
                    case "phase"
                        ylabel('phase (°)');
                    case "real"
                        ylabel('real part');
                    case "imaginary"
                        ylabel('imaginary part');
                end
                [idTx, idRx] = meshgrid(1 : obj.numberOfActiveTransmittingNodes, 1 : obj.numberOfActiveReceivingNodes);
                leg = legend(num2str([idTx(:) idRx(:)]), 'Location', 'best');
                title(leg, 'TX | RX ID');
            else
                m.FaceColor = 'flat';
                xlabel('doppler bins'); ylabel('delay samples');
                switch options.plot
                    case "magnitude"
                        zlabel('power (dB)');
                    case "phase"
                        zlabel('phase (°)');
                    case "real"
                        zlabel('real part');
                    case "imaginary"
                        zlabel('imaginary part');
                end
                zlim([minMag 0]);
            end
            hold off; drawnow;
        end

        function visualizenetwork(obj, options)
            arguments
                obj
                options.figureID double {mustBeInteger, mustBeNonnegative} = []
            end
            if isempty(options.figureID)
                figure;
            else
                figure(options.figureID);
            end
            posRX = [obj.receivingNodes.position]/1e3;
            posTX = [obj.transmittingNodes.position]/1e3;
            plot3(posRX(1, :), posRX(2, :), posRX(3, :), 'vb', 'LineWidth', 3);
            hold on; plot3(posTX(1, :), posTX(2, :), posTX(3, :), 'vr', 'LineWidth', 3);
            text(posRX(1, :), posRX(2, :), posRX(3, :), num2str((1 : obj.numberOfActiveReceivingNodes).'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            text(posTX(1, :), posTX(2, :), posTX(3, :), num2str((1 : obj.numberOfActiveTransmittingNodes).'), "FontSize", 20, "FontWeight", "bold", "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
            for rxID = 1 : obj.numberOfReceivingNodes
                n = obj.receivingNodes(rxID).array.normalVector;
                quiver3(posRX(1, rxID), posRX(2, rxID), posRX(3, rxID), n(1), n(2), n(3), 'b');
            end
            for txID = 1 : obj.numberOfTransmittingNodes
                n = obj.transmittingNodes(txID).array.normalVector;
                quiver3(posTX(1, txID), posTX(2, txID), posTX(3, txID), n(1), n(2), n(3), 'r');
            end
            posRX = repelem(posRX, 1, obj.numberOfTransmittingNodes);
            posTX = repmat(posTX, 1, obj.numberOfReceivingNodes);
            line([posRX(1, :); posTX(1, :)], [posRX(2, :); posTX(2, :)], [posRX(3, :); posTX(3, :)], 'lineStyle', '--', 'Color', 'k');
            grid off; grid on; grid minor; title('radar network configuration');
            xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
            legend('RX', 'TX', 'Location', 'best');
            ax = gca; lims = [ax.XLim; ax.YLim; ax.ZLim];
            lims = [min(lims(:)), max(lims(:))];
            xlim(lims); ylim(lims); zlim(lims);
            hold off; drawnow;
        end
    end
end