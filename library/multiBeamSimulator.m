classdef multiBeamSimulator < handle

    properties (SetAccess = private, GetAccess = public)
        fusionCenter (1, 1) fusionCenter = fusionCenter.empty()
    end

    properties (Access = private)
        isWait (1, 1) logical {mustBeNumericOrLogical, mustBeMember(isWait, [0, 1])} = 0;
        isStep (1, 1) logical {mustBeNumericOrLogical, mustBeMember(isStep, [0, 1])} = 0;
        isRestart (1, 1) logical {mustBeNumericOrLogical, mustBeMember(isRestart, [0, 1])} = 0;
        isSimulation (1, 1) logical {mustBeNumericOrLogical, mustBeMember(isSimulation, [0, 1])} = 1;
    end

    methods
        function obj = multiBeamSimulator(options)
            arguments
                options.fusionCenter (1, 1) fusionCenter
            end
            %fusionCenter Construct an instance of this class
            %   Detailed explanation goes here
            obj.fusionCenter = options.fusionCenter;
        end

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
        
        % function detection = startsimulation(obj, options)
        %     arguments
        %         obj
        %         options.duration (1, 1) double {mustBePositive} = 1
        %     end
        %     stepTime = obj.network.beamTime;
        %     numberOfSteps = ceil(options.duration ./stepTime);
        %     currentNumberOfTrials = obj.monteCarlo.numberOfTrials;
        %     currentNumberOfTrialsParallel = obj.monteCarlo.numberOfTrialsParallel;
        %     cleanup = onCleanup(@() recoverstate(obj, currentNumberOfTrials, currentNumberOfTrialsParallel));
        %     obj.configure("numberOfTrials", 1, "numberOfTrialsParallel", 1);
        %     obj.resetsimulation;
        %     figID = 19234;
        %     h = figure(figID);
        %     if exist('h', 'var') && ishandle(h)
        %         close(19234);
        %     end
        %     obj.isWait = 0;
        %     obj.isStep = 0;
        %     obj.isRestart = 0;
        %     obj.isSimulation = 1;
        %     S.f = figure(figID); hold on; view(0, 90);
        %     set(S.f, 'WindowState', 'maximized');
        %     S.pb = uicontrol( ...
        %         'Style', 'pushbutton', ...
        %         'Units', 'normalized', ...
        %         'Position', [0.92 0.85 0.05 0.05], ...
        %         'FontSize', 14);
        %     S.sb = uicontrol( ...
        %         'Style', 'pushbutton', ...
        %         'Units', 'normalized', ...
        %         'Position', [0.92 0.78 0.05 0.05], ...
        %         'FontSize', 14);
        %     S.rb = uicontrol( ...
        %         'Style', 'pushbutton', ...
        %         'Units', 'normalized', ...
        %         'Position', [0.92 0.71 0.05 0.05], ...
        %         'FontSize', 14);
        %     S.eb = uicontrol( ...
        %         'Style', 'pushbutton', ...
        %         'Units', 'normalized', ...
        %         'Position', [0.92 0.64 0.05 0.05], ...
        %         'FontSize', 14);
        %     set(S.pb, 'Callback', {@obj.pb_call, S});
        %     set(S.sb, 'Callback', {@obj.sb_call, S});
        %     set(S.rb, 'Callback', {@obj.rb_call, S});
        %     set(S.eb, 'Callback', {@obj.eb_call, S});
        %     set(S.pb, 'String', char(9208));
        %     set(S.sb, 'String', char(8618));
        %     set(S.rb, 'String', char(9198));
        %     set(S.eb, 'String', char(10006));
        %     posRX = [obj.network.activeReceivingNodes.position]/1e3;
        %     posTX = [obj.network.activeTransmittingNodes.position]/1e3;
        %     posRXrep = repelem(posRX, 1, obj.interfaces.numberOfTransmittingNodes);
        %     posTXrep = repmat(posTX, 1, obj.interfaces.numberOfReceivingNodes);
        % 
        %     numberOfActiveTransmittingNodes = obj.network.numberOfActiveTransmittingNodes;
        %     numberOfActiveReceivingNodes = obj.network.numberOfActiveReceivingNodes;
        % 
        %     colors = {'m', 'c'};
        %     phi = linspace(-180, 180, 1000);
        %     d = max(sqrt(sum(obj.network.boundary.^2)));
        %     interfacesTX = interface.empty;
        %     interfacesRX = interface.empty;
        %     visualizationNormalizationTX = zeros(1, numberOfActiveTransmittingNodes);
        %     for txID = 1 : numberOfActiveTransmittingNodes
        %         targetsTemp = target( ...
        %             'position', [ ...
        %             d*cosd(phi) + posTX(1, txID)*1e3; ...
        %             d*sind(phi) + posTX(2, txID)*1e3; ...
        %             zeros(1, 1000) + posTX(3, txID)*1e3]);
        %         interfacesTX(txID) = obj.network.getinterface(targetsTemp);
        %         visualizationNormalizationTX(txID) = max(interfacesTX(txID).distanceTX(txID, :, :), [], 'all')/1e3;
        %     end
        %     visualizationNormalizationRX = zeros(1, numberOfActiveReceivingNodes);
        %     for rxID = 1 : numberOfActiveReceivingNodes
        %         targetsTemp = target( ...
        %             'position', [ ...
        %             d*cosd(phi) + posRX(1, rxID)*1e3; ...
        %             d*sind(phi) + posRX(2, rxID)*1e3; ...
        %             zeros(1, 1000) + posRX(3, rxID)*1e3]);
        %         interfacesRX(rxID) = obj.network.getinterface(targetsTemp);
        %         visualizationNormalizationRX(rxID) = max(interfacesRX(rxID).distanceRX(:, rxID, :), [], 'all')/1e3;
        %     end
        % 
        %     plot3(posRX(1, :), posRX(2, :), posRX(3, :), 'xb', 'LineWidth', 2, 'MarkerSize', 5);
        %     plot3(posTX(1, :), posTX(2, :), posTX(3, :), '+r', 'LineWidth', 2, 'MarkerSize', 5);
        %     line([posRXrep(1, :); posTXrep(1, :)], [posRXrep(2, :); posTXrep(2, :)], [posRXrep(3, :); posTXrep(3, :)], 'lineStyle', '--', 'Color', 'k');
        % 
        %     ax = gca;
        %     lines = ax.Children;
        % 
        %     grid off; grid on; grid minor; title('scenario');
        %     xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
        %     xlim(obj.network.boundaryListened(1, :)/1e3);
        %     ylim(obj.network.boundaryListened(2, :)/1e3);
        %     zlim(obj.network.boundaryListened(3, :)/1e3);
        % 
        %     switch obj.network.networkMode
        %         case "multiStatic"
        %             detection = struct( ...
        %                 'power', [], ...
        %                 'x', [], ...
        %                 'y', [], ...
        %                 'z', [], ...
        %                 'numberOfDetections', [], ...
        %                 'position', []);
        %             detection = repmat(detection, 1, numberOfSteps);
        %         case "monoStatic"
        %             detection = struct( ...
        %                 'power', [], ...
        %                 'range', [], ...
        %                 'elevation', [], ...
        %                 'azimuth', [], ...
        %                 'numberOfDetections', [], ...
        %                 'position', []);
        %             detection = repmat({repmat(detection, 1, numberOfSteps)}, obj.network.numberOfActiveReceivingNodes, 1);
        %     end
        %     numberOfLinesInitial = length(ax.Children);
        %     stepID = 0;
        %     start = true;
        %     while obj.isSimulation
        %         if stepID == numberOfSteps
        %             set(S.pb, 'String', char(9205));
        %             uiwait(S.f);
        %         end
        %         if obj.isRestart && stepID
        %             obj.isRestart = 0;
        %             obj.resetsimulation;
        %             stepID = 0;
        %             if ~obj.isStep
        %                 obj.isWait = 0;
        %                 set(S.pb, 'String', char(9208));
        %             end
        %             continue
        %         end
        %         if stepID == numberOfSteps
        %             continue
        %         end
        %         figure(figID);
        %         lineID = numberOfLinesInitial + 1;
        %         if ~start
        %             if ~isempty(obj.interfaces.targets)
        %                 obj.interfaces.targets.step(stepTime);
        %             end
        %             for rxID = 1 : numberOfActiveReceivingNodes
        %                 obj.network.activeReceivingNodes(rxID).array.step(stepTime);
        %             end
        %             for txID = 1 : numberOfActiveTransmittingNodes
        %                 obj.network.transmittingNodes(txID).array.step(stepTime);
        %             end
        %         end
        %         obj.applyspatialprocessing;
        %         switch obj.network.networkMode
        %             case "multiStatic"
        %                 detections = obj.detectionFromIntegration;
        %             case "monoStatic"
        %                 detections = obj.detectionFromMatchFiltration;
        %         end
        %         if ~start && sum([detections.numberOfDetections]) == 0
        %             continue;
        %         end
        %         if (obj.isWait || obj.isStep) && stepID
        %             uiwait(S.f);
        %             if obj.isRestart
        %                 obj.isRestart = 0;
        %                 obj.resetsimulation;
        %                 stepID = 0;
        %                 continue
        %             end
        %         end
        % 
        %         if ~isempty(obj.interfaces.targets)
        %             x = obj.interfaces.targets.position(1, :)/1e3;
        %             y = obj.interfaces.targets.position(2, :)/1e3;
        %             z = obj.interfaces.targets.position(3, :)/1e3;
        % 
        %             if start
        %                 if obj.interfaces.numberOfTargets < 11
        %                     plot3(x, y, z, '+k', 'LineWidth', 2, 'MarkerSize', 10);
        %                 elseif obj.interfaces.numberOfTargets < 21
        %                     plot3(x, y, z, '*k', 'LineWidth', 1, 'MarkerSize', 10);
        %                 else
        %                     plot3(x, y, z, '.k');
        %                 end
        %             else
        %                 lines(lineID).XData = x;
        %                 lines(lineID).YData = y;
        %                 lines(lineID).ZData = z;
        %                 lineID = lineID + 1;
        %             end
        %         end
        % 
        %         switch obj.network.networkMode
        %             case "monoStatic"
        %                 for rxID = 1 : obj.network.numberOfActiveReceivingNodes
        %                     detection{rxID}(stepID + 1) = detections(rxID);
        %                     pos = detection{rxID}(stepID + 1).position;
        %                     if isempty(pos)
        %                         if start
        %                             % dummy to create line
        %                             if obj.network.numberOfActiveReceivingNodes == 2
        %                                 plot3(nan, nan, nan, ['.', colors{rxID}], 'LineWidth', 1, 'MarkerSize', 10);
        %                             else
        %                                 plot3(nan, nan, nan, '.c', 'LineWidth', 1, 'MarkerSize', 10);
        %                             end
        %                             continue;
        %                         else
        %                             lines(lineID).XData = [];
        %                             lines(lineID).YData = [];
        %                             lines(lineID).ZData = [];
        %                             lineID = lineID + 1;
        %                             continue
        %                         end
        %                     end
        %                     xEstimation = pos(1, :)/1e3;
        %                     yEstimation = pos(2, :)/1e3;
        %                     zEstimation = pos(3, :)/1e3;
        %                     if start
        %                         if obj.network.numberOfActiveReceivingNodes == 2
        %                             plot3(xEstimation, yEstimation, zEstimation, ['.', colors{rxID}], 'LineWidth', 1, 'MarkerSize', 10);
        %                         else
        %                             plot3(xEstimation, yEstimation, zEstimation, '.c', 'LineWidth', 1, 'MarkerSize', 10);
        %                         end
        %                     else
        %                         lines(lineID).XData = xEstimation;
        %                         lines(lineID).YData = yEstimation;
        %                         lines(lineID).ZData = zEstimation;
        %                         lineID = lineID + 1;
        %                     end
        %                 end
        %             case "multiStatic"
        %                 detection(stepID + 1) = detections;
        %                 pos = detection(stepID + 1).position;
        %                 if isempty(pos)
        %                     if start
        %                         % dummy to create line
        %                         plot3(nan, nan, nan, '.r', 'LineWidth', 1, 'MarkerSize', 10);
        %                         continue;
        %                     else
        %                         lines(lineID).XData = [];
        %                         lines(lineID).YData = [];
        %                         lines(lineID).ZData = [];
        %                         continue
        %                     end
        %                 end
        %                 xEstimation = pos(1, :)/1e3;
        %                 yEstimation = pos(2, :)/1e3;
        %                 zEstimation = pos(3, :)/1e3;
        %                 if start
        %                     plot3(xEstimation, yEstimation, zEstimation, '.r', 'LineWidth', 1, 'MarkerSize', 10);
        %                 else
        %                     lines(lineID).XData = xEstimation;
        %                     lines(lineID).YData = yEstimation;
        %                     lines(lineID).ZData = zEstimation;
        %                 end
        %         end
        % 
        %         if (~strcmpi(obj.network.surveillanceMode, "staticBeam") && ~start) || start
        %             for txID = 1 : numberOfActiveTransmittingNodes
        %                 G = visualizationNormalizationTX(txID).*permute(abs(interfacesTX(txID).transmittedBeam(txID, :, :).^2./interfacesTX(txID).network.activeTransmittingNodes(txID).array.numberOfTotalElements), [1 3 2]);
        %                 u = permute(interfacesTX(txID).unitDirectionTX, [1 3 2]);
        %                 if start
        %                     plot3(G.*u(1, :, txID) + posTX(1, txID), G.*u(2, :, txID) + posTX(2, txID), u(3, :, txID) + posTX(3, txID), 'r');
        %                 else
        %                     lines(lineID).XData = G.*u(1, :, txID) + posTX(1, txID);
        %                     lines(lineID).YData = G.*u(2, :, txID) + posTX(2, txID);
        %                     lines(lineID).ZData = u(3, :, txID) + posTX(3, txID);
        %                     lineID = lineID + 1;
        %                 end
        %             end
        %             for rxID = 1 : numberOfActiveReceivingNodes
        %                 switch obj.network.networkMode
        %                     case "monoStatic"
        %                         txID = obj.network.monoStaticTransmitterIDs(rxID);
        %                         if isnan(txID)
        %                             break;
        %                         end
        %                     case "multiStatic"
        %                         txID = 1;
        %                 end
        %                 G = visualizationNormalizationRX(rxID).*permute(abs(interfacesRX(rxID).receivedBeamSpaceObservations{rxID}(txID, :, :, :).^2./interfacesRX(rxID).network.activeReceivingNodes(rxID).array.numberOfTotalElements), [3 4 1 2]);
        %                 u = permute(interfacesRX(rxID).unitDirectionRX, [3 1 2]);
        %                 if start
        %                     if obj.network.numberOfActiveReceivingNodes == 2
        %                         plot3(G.*u(:, 1, rxID) + posRX(1, rxID), G.*u(:, 2, rxID) + posRX(2, rxID), u(:, 3, rxID) + posRX(3, rxID), colors{rxID}, 'LineWidth', 2);
        %                     else
        %                         plot3(G.*u(:, 1, rxID) + posRX(1, rxID), G.*u(:, 2, rxID) + posRX(2, rxID), u(:, 3, rxID) + posRX(3, rxID), 'b', 'LineWidth', 2);
        %                     end
        %                 else
        %                     numberOfTotalChannels = obj.network.activeReceivingNodes(rxID).numberOfTotalChannels;
        %                     for chID = 1 : numberOfTotalChannels
        %                         lines(lineID).XData = G(:, chID).*u(:, 1, rxID) + posRX(1, rxID);
        %                         lines(lineID).YData = G(:, chID).*u(:, 2, rxID) + posRX(2, rxID);
        %                         lines(lineID).ZData = u(:, 3, rxID) + posRX(3, rxID);
        %                         lineID = lineID + 1;
        %                     end
        %                 end
        %             end
        %             for txID = 1 : numberOfActiveTransmittingNodes
        %                 n = visualizationNormalizationTX(txID).*obj.network.activeTransmittingNodes(txID).array.normalVector/2;
        %                 if start
        %                     quiver3(posTX(1, txID), posTX(2, txID), posTX(3, txID), n(1), n(2), n(3), 'Color', 'k');
        %                 else
        %                     lines(lineID).XData = posTX(1, txID);
        %                     lines(lineID).YData = posTX(2, txID);
        %                     lines(lineID).ZData = posTX(3, txID);
        %                     lines(lineID).UData = n(1);
        %                     lines(lineID).VData = n(2);
        %                     lines(lineID).WData = n(3);
        %                     lineID = lineID + 1;
        %                 end
        %             end
        %             for rxID = 1 : numberOfActiveReceivingNodes
        %                 n = visualizationNormalizationRX(rxID).*obj.network.activeReceivingNodes(rxID).array.normalVector/2;
        %                 if start
        %                     quiver3(posRX(1, rxID), posRX(2, rxID), posRX(3, rxID), n(1), n(2), n(3), 'Color', 'k');
        %                 else
        %                     lines(lineID).XData = posRX(1, rxID);
        %                     lines(lineID).YData = posRX(2, rxID);
        %                     lines(lineID).ZData = posRX(3, rxID);
        %                     lines(lineID).UData = n(1);
        %                     lines(lineID).VData = n(2);
        %                     lines(lineID).WData = n(3);
        %                     lineID = lineID + 1;
        %                 end
        %             end
        %         end
        %         % for txID = 1 : numberOfActiveTransmittingNodes
        %         %     s = 1.5*visualizationNormalizationTX(txID).*obj.network.activeTransmittingNodes(txID).array.steeringUnitDirection;
        %         %     if start
        %         %         quiver3(posTX(1, txID), posTX(2, txID), posTX(3, txID), s(1), s(2), s(3), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
        %         %     else
        %         %         lines(lineID).XData = posTX(1, txID);
        %         %         lines(lineID).YData = posTX(2, txID);
        %         %         lines(lineID).ZData = posTX(3, txID);
        %         %         lines(lineID).UData = s(1);
        %         %         lines(lineID).VData = s(2);
        %         %         lines(lineID).WData = s(3);
        %         %         lineID = lineID + 1;
        %         %     end
        %         % end
        %         % for rxID = 1 : numberOfActiveReceivingNodes
        %         %     s = 1.5*visualizationNormalizationRX(rxID).*obj.network.activeReceivingNodes(rxID).array.steeringUnitDirection;
        %         %     if start
        %         %         if obj.network.numberOfActiveReceivingNodes == 2
        %         %             quiver3(posRX(1, rxID), posRX(2, rxID), posRX(3, rxID), s(1), s(2), s(3), 'Color', colors{rxID}, 'LineStyle', '--', 'LineWidth', 2);
        %         %         else
        %         %             quiver3(posRX(1, rxID), posRX(2, rxID), posRX(3, rxID), s(1), s(2), s(3), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
        %         %         end
        %         %     else
        %         %         lines(lineID).XData = posRX(1, rxID);
        %         %         lines(lineID).YData = posRX(2, rxID);
        %         %         lines(lineID).ZData = posRX(3, rxID);
        %         %         lines(lineID).UData = s(1);
        %         %         lines(lineID).VData = s(2);
        %         %         lines(lineID).WData = s(3);
        %         %         lineID = lineID + 1;
        %         %     end
        %         % end
        % 
        %         if start
        %             lines = flipud(ax.Children);
        %             start = false;
        %         end
        %         if stepID < numberOfSteps
        %             stepID = stepID + 1;
        %         end
        %         drawnow;
        %     end
        %     function recoverstate(obj, currentNumberOfTrials, currentNumberOfTrialsParallel)
        %         obj.configure( ...
        %             "numberOfTrials", currentNumberOfTrials, ...
        %             "numberOfTrialsParallel", currentNumberOfTrialsParallel);
        %         obj.resetsimulation;
        %     end
        % end
        % 
        % function resetsimulation(obj)
        %     if ~isempty(obj.interfaces.targets)
        %         obj.interfaces.targets.reset;
        %     end
        %     for rxID = 1 : obj.network.numberOfActiveReceivingNodes
        %         obj.network.activeReceivingNodes(rxID).array.reset;
        %     end
        %     for txID = 1 : obj.network.numberOfActiveTransmittingNodes
        %         obj.network.transmittingNodes(txID).array.reset;
        %     end
        % end
    end
end