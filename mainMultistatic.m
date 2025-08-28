%%% Anıl ARSLAN 2303980
% clc; clear; % close all;
addpath(genpath([pwd '/library']));
set(groot, "defaultFigurePosition", [680 458 560 420]);

d = 10e3;

arrayRX = uniformPlanarArray( ...
    "numberOfElements", [1 16], ...
    "rpm", 0);
arrayTX = uniformPlanarArray( ...
    "numberOfElements", [1 16], ...
    "rpm", 0);

% M = 9;
% theta = [-90 -67.5 -45 -22.5 0 22.5 45 67.5 90] - 90;
% x = d.*cosd(theta);
% y = d.*sind(theta) + d;
% Rm = [x; y; zeros(1, M)];

Rm = (-3 : 3)*d;
M = size(Rm, 2);
Rm = [Rm; -d*ones(1, M); zeros(1, M)];
Rm(2, 1) = 0;
Rm(2, 7) = 0;
Rm(2, 4) = -1.5*d;
Rm(2, 3) = d;
Rm(2, 5) = d;
Rm(1, 3) = -2*d;
Rm(1, 5) = 2*d;

receivers = receivingNode( ...
    'position', Rm, ...
    'samplingFrequency', 1e6, ...
    'CPIsecond', 0.25*1e-3, ...
    'array', arrayRX ...
    );

transmitters = transmittingNode( ...
    'position', [0; 1.5*d; 0], ...
    'array', arrayTX, ...
    'peakPower_W', 40e3, ...
    'pulseWidth', 10e-6, ...
    'carrierFrequency', 1e9);
transmitters.setLFM( ...
    "bandWidth", 1e6, ...
    "frequencyDirection", "symmetric");

receiverArrays = [receivers.array];
% receiverArrays.setorientation("yaw", [30 45 90 135 150]);
% receiverArrays.setscanparameters("backOfArrayRegion", 90);
% transmitters.array.setorientation("yaw", 90);
% transmitters.array.setscanparameters("backOfArrayRegion", 240);


% transmitters.array.steer("steeringAzimuth", atand(-14050/10100));
% receiverArrays(4).steer("steeringAzimuth", atand(-14050/10100));
% receiverArrays.setscanparameters("backOfArray", true);
% transmitters.array.setscanparameters("backOfArray", true);


% previous
% receiverArrays.setorientation("yaw", [15 -150 165 -90 -75 -15 -105]);
% receiverArrays.setscanparameters("backOfArrayRegion", 150);
% transmitters.array.setorientation("yaw", -90);
% transmitters.array.setscanparameters("backOfArrayRegion", 150);
% transmitters.array.steer("steeringAzimuth", atand(5/25));
% receiverArrays(4).steer("steeringAzimuth", atand(5/25));

% just previous
% receiverArrays.setscanparameters("backOfArrayRegion", 180);
% receiverArrays.setorientation("yaw", [0 -180 180 -90 -90 0 -90]);
% transmitters.array.setscanparameters("backOfArrayRegion", 180);
% transmitters.array.setorientation("yaw", -90);


receiverArrays.setscanparameters("backOfArray", true);
transmitters.array.setscanparameters("backOfArray", true);


network = radarNetwork( ...
    'receivingNodes', receivers, ...
    'transmittingNodes', transmitters);
network.settingsnetwork( "networkMode", "multiStatic");
% network.settingsnetwork( "networkMode", "monoStatic");
network.activitycontrol( ...
    "transmittingNodeActivity", 1, ...
    "receivingNodeActivity", ones(1, M));
network.setsurveillance("staticBeam");
network.setsynchronization( ...
    "timeOffsetBistaticPairs", 0, ...
    ..."timeDeviationBistaticPairs", 5e-8, ... Ts/10
    ..."timeDeviationBistaticPairs", 5e-7, ... Ts
    "frequencyOffsetBistaticPairs", 0, ...
    "frequencyDeviationBistaticPairs", 0 ... 5e8
    );

targets = target( ...
    'position', [10200; -15100; 0], ...
    'meanRCS_dbsm', 10*log10(5));
int = interface( ...
    'network', network, ...
    'targets', targets);
fc = fusionCenter( ...
    "interfaces", int, ...
    "gridResolution", [150 150 inf]);

%%% Configuration %%%

int.configure( ...
    'noise', 1, ...
    'directPath', 0, ...
    "pathLoss", 1, ...
    'fractionalDelay', 'off', ...
    ...'fractionalDelay', 'sinc-based', ...
    "timeSynchronization", 1, ...
    "frequencySynchronization", 1, ...
    'spatialCoherency', 'noncoherent', ...
    ...'spatialCoherency', 'coherent', ...
    'swerling', 1);
fc.configure( ...
    "globalPFA", 1e-6, ...
    "globalIntegrationRule", "SLC", ...
    "localPFA", 1, ...
    "seed", 0, ...
    "seedShuffle", 0, ...
    "numberOfTrials", 1, ...
    'numberOfTrialsParallel', 1, ...
    'numberOfGuardCells', 2);

%% Dictionary
fc.setdictionary;

%% Visualization %%%
arrayTX.visualizearray;
transmitters.visualizewaveform;
fig = network.visualizenetwork;
xlim([-50, 50]); ylim([-40, 50]);
figureName = 'scenario_1';
savefig(fig, ['C:\GitRepo\MSRS\figuresSim\' figureName '.fig']);
saveas(fig, ['C:\GitRepo\MSRS\figuresSim\' figureName '.eps'], 'epsc');
network.visualizewaveformsampled;
fig = int.visualizescenario("showPattern", 1, "figureID", 578);
int.visualizeellipses("ellipseType", "target");
int.visualizeellipses("ellipseType", "resolution"); xlim([9 11.5]); ylim([-16.5 -14]);
% figureName = 'scenario_1_resolution'; fig2 = gcf;
% savefig(fig2, ['C:\GitRepo\MSRS\figuresSim\' figureName '.fig']);
% saveas(fig2, ['C:\GitRepo\MSRS\figuresSim\' figureName '.eps'], 'epsc');
int.visualizereceivedsignals("receivingNodeIDs", 1 : 7, "trialID", 1); % draws different trial all the time
int.visualizebeamformedsignals("receivingNodeIDs", 1, "trialID", 1); % draws different trial all the time
fc.visualizehypothesizedvariables("receivingNodeID", 1, "variable", "averageSNR");

%% Single Target Simulation

targets = target( ...
    'position', [10200; -15100; 0], ...
    'meanRCS_dbsm', 10*log10(1) ...
    );
% % sampling straddle loss is not known
int.settargets(targets)
fc.configure( ...
    "seed", 0, ...
    "globalIntegrationRule", "BC", ...
    "localPFA", 1e-2);
fc.applyspatialprocessing("saveSpectrum", 1, "doesPrint", true);

% Visualization %%%
trialID = 1; % coh kontrol edilmeli

% fc.visualizefilteredsignals("receivingNodeIDs", 3 : 4, "trialID", trialID); % draws saved data
% xlim([-9 inf]); ylim([-20 20*log10(10)]);
% fc.visualizeestimation("trialID", trialID);
fig2 = fc.visualizeintegratedsignals( ...
    "trialID", trialID, ...
    "plotMode", "image");
% clim([0, 20]);
% figureName = 'scenario_1_response_localPFA_BC';
% savefig(fig2, ['C:\GitRepo\MSRS\figuresSim\' figureName '.fig']);
% saveas(fig2, ['C:\GitRepo\MSRS\figuresSim\' figureName '.eps'], 'epsc');

%% Resolution simulation

% targets = target( ...
%     'position', [-120; 50; 0], ... straddle loss
%     'position', [0; 1; 0], ... straddle loss
%     'position', [-150; 550; 0], ... straddle loss
%     'position', [500; 500; 0], ... on cell w/out straddle
%     ...'position', [0; 600; 0], ... on cell w/out straddle
%     ...'position', [-250; 750; 0], ... straddle loss
%     'meanRCS_dbsm', -4, ...
%     'meanRCS_dbsm', 0 ...
%     );
    %%% sampling straddle loss is not known
% int.settargets(targets)
% clc;
fc.configure( ...
    "seed", 0, ...
    'localPFA', 1, ... 0.1
    "globalIntegrationRule", "SLC", ...
    "numberOfTrials", 1, ...
    "numberOfTrialsParallel", 1000 ...
    );
fc.simulatedetection( ...
    "randomOnCell", 1, ...
    "numberOfTargets", 1 ...
    );
% close all;
fc.visualizedetectionsimulation;

%% analytical coverage

clc;
fc.configure( ...
    "globalIntegrationRule", "BC", ...
    'localPFA', 1e-3);

fc.analyticalcoverage("meanRCS_dbsm", 10*log10(1));
close all;
fc.visualizecoveragesimulation( ...
    "contourLevelDetection", 0.85, ...
    "saveFigure", 1, "header", 'scenario_1_coverage_BC');
% fc.visualizecoveragesimulation( ...
%     "contourLevelDetection", 0.7, ...
%     "saveFigure", 0, "header", 'scenario_1_coverage_BC');

%% analytical coverage BC

fc.configure("globalIntegrationRule", "SLC");
localPFA = logspace(0, -8, 100);
localPFA = 1;

gridScan = fc.gridPointsMesh;
cellPositions = reshape(permute(cat(4, gridScan.x, gridScan.y, gridScan.z), [4 1 2 3]), [3 prod(fc.gridSize)]);
cellPositions = cellPositions(:, ~fc.blindZone);
globalPDmodel = zeros(prod(fc.gridSize), length(localPFA));
fc.interfaces.settargets(target( ...
    "position", cellPositions, ...
    "meanRCS_dbsm", 0));
dtc = fc.configuration.detector;
dtc.localPFA{1} = repmat(permute(localPFA, [1 3 4 2]), 7, 1);
globalPDmodel(~fc.blindZone, :) = permute(dtc.globalPD, [2 4 1 3]);


dimensions = {"x", "y", "z"};
dimensions = dimensions(size(gridScan.x) ~= 1);
xLabel = dimensions{1} + " (km)";
yLabel = dimensions{2} + " (km)";
x1 = fc.gridPoints{1}/1e3;
x2 = fc.gridPoints{2}/1e3;

globalPDmodel = reshape(globalPDmodel, [fc.gridSize([2 1]), length(localPFA)]);

%%
colors = rgb2hex([ ...
    0.0000    0.4470    0.7410 % blue
    0.8500    0.3250    0.0980 % orange
    0.9290    0.6940    0.1250 % yellow
    0.4940    0.1840    0.5560 % magenta
    0.4660    0.6740    0.1880 % green
    0.3010    0.7450    0.9330 % cyan
    0.6350    0.0780    0.1840]);
load('globalPDmodel_BC.mat', 'globalPDmodel');
globalPDmodel_BC = globalPDmodel;
load('globalPDmodel_SLC_localPFA.mat', 'globalPDmodel');
globalPDmodel_SLC_localPFA = globalPDmodel;
load('globalPDmodel_SLC.mat', 'globalPDmodel');
globalPDmodel_SLC = globalPDmodel;

fig85 = figure(1985); hold on;
% for i = 1 : length(localPFA)
%     contour(x1, x2, globalPDmodel(:, :, i), [-1 0.85], 'LineWidth', 2, 'EdgeColor', [0.93 0.69 0.13]);
% end
contour(x1, x2, globalPDmodel_SLC, [-1 0.85], 'LineWidth', 2, 'EdgeColor', colors(1));
contour(x1, x2, globalPDmodel_SLC_localPFA, [-1 0.85], 'LineWidth', 2, 'EdgeColor', colors(2));
contour(x1, x2, globalPDmodel_BC(:, :, 86), [-1 0.85], 'LineWidth', 2, 'EdgeColor', colors(3));
contour(x1, x2, globalPDmodel_BC(:, :, 47), [-1 0.85], 'LineWidth', 2, 'EdgeColor', colors(4));
contour(x1, x2, globalPDmodel_BC(:, :, 13), [-1 0.85], 'LineWidth', 2, 'EdgeColor', colors(5));
grid off; grid on; grid minor;
xlabel(xLabel); ylabel(yLabel); zlabel('p_D');
legend("centralized SLC", "distributed SLC", "BC-\lambda^{local} = " + scinot(localPFA(86)), "BC-\lambda^{local} = " + scinot(localPFA(47)), "BC-\lambda^{local} = " + scinot(localPFA(13)), "Location", "best");
figureName = 'scenario_1_coverage_contour_85';
savefig(fig85, ['C:\GitRepo\MSRS\figuresSim\' figureName '.fig']);
saveas(fig85, ['C:\GitRepo\MSRS\figuresSim\' figureName '.eps'], 'epsc');


fig70 = figure(1970); hold on;
% for i = 1 : length(localPFA)
%     contour(x1, x2, globalPDmodel(:, :, i), [-1 0.70], 'LineWidth', 2, 'EdgeColor', [0.93 0.69 0.13]);
% end
contour(x1, x2, globalPDmodel_SLC, [-1 0.70], 'LineWidth', 2, 'EdgeColor', colors(1));
contour(x1, x2, globalPDmodel_SLC_localPFA, [-1 0.75], 'LineWidth', 2, 'EdgeColor', colors(2));
contour(x1, x2, globalPDmodel_BC(:, :, 86), [-1 0.70], 'LineWidth', 2, 'EdgeColor', colors(3));
contour(x1, x2, globalPDmodel_BC(:, :, 47), [-1 0.70], 'LineWidth', 2, 'EdgeColor', colors(4));
contour(x1, x2, globalPDmodel_BC(:, :, 13), [-1 0.70], 'LineWidth', 2, 'EdgeColor', colors(5));
grid off; grid on; grid minor;
xlabel(xLabel); ylabel(yLabel); zlabel('p_D');
legend("centralized SLC", "distributed SLC", "BC-\lambda^{local} = " + scinot(localPFA(86)), "BC-\lambda^{local} = " + scinot(localPFA(47)), "BC-\lambda^{local} = " + scinot(localPFA(13)), "Location", "best");
figureName = 'scenario_1_coverage_contour_70';
savefig(fig70, ['C:\GitRepo\MSRS\figuresSim\' figureName '.fig']);
saveas(fig70, ['C:\GitRepo\MSRS\figuresSim\' figureName '.eps'], 'epsc');

fig1 = figure;
img = imagesc(x1, x2, globalPDmodel_BC(:, :, 86));
colorbar; colormap('gray'); clim([0 1]);
ax = gca; set(ax, 'Ydir', 'Normal');
% set(img, 'AlphaData', visibleZone);
xlabel(xLabel); ylabel(yLabel);
figureName = 'scenario_1_coverage_BC_max';
savefig(fig1, ['C:\GitRepo\MSRS\figuresSim\' figureName '.fig']);
saveas(fig1, ['C:\GitRepo\MSRS\figuresSim\' figureName '.eps'], 'epsc');

fig1 = figure;
img = imagesc(x1, x2, globalPDmodel_BC(:, :, 47));
colorbar; colormap('gray'); clim([0 1]);
ax = gca; set(ax, 'Ydir', 'Normal');
% set(img, 'AlphaData', visibleZone);
xlabel(xLabel); ylabel(yLabel);
figureName = 'scenario_1_coverage_BC_submax';
savefig(fig1, ['C:\GitRepo\MSRS\figuresSim\' figureName '.fig']);
saveas(fig1, ['C:\GitRepo\MSRS\figuresSim\' figureName '.eps'], 'epsc');


fig1 = figure;
img = imagesc(x1, x2, globalPDmodel_SLC_localPFA);
colorbar; colormap('gray'); clim([0 1]);
ax = gca; set(ax, 'Ydir', 'Normal');
% set(img, 'AlphaData', visibleZone);
xlabel(xLabel); ylabel(yLabel);
figureName = 'scenario_1_coverage_SLC_local_PFA';
savefig(fig1, ['C:\GitRepo\MSRS\figuresSim\' figureName '.fig']);
saveas(fig1, ['C:\GitRepo\MSRS\figuresSim\' figureName '.eps'], 'epsc');

fig1 = figure;
img = imagesc(x1, x2, globalPDmodel_SLC);
colorbar; colormap('gray'); clim([0 1]);
ax = gca; set(ax, 'Ydir', 'Normal');
% set(img, 'AlphaData', visibleZone);
xlabel(xLabel); ylabel(yLabel);
figureName = 'scenario_1_coverage_SLC';
savefig(fig1, ['C:\GitRepo\MSRS\figuresSim\' figureName '.fig']);
saveas(fig1, ['C:\GitRepo\MSRS\figuresSim\' figureName '.eps'], 'epsc');

%% simulated coverage

% clc;
% fc.configure( ...
%     "seed", 0, ...
%     'localPFA', 1e-3, ... 0.1
%     "numberOfTrials", 1, ...
%     "numberOfTrialsParallel", 1);
% 
% fc.simulatecoverage("meanRCS_dbsm", 10*log10(5), "onCellCenters", 1, "neighbourOffset", 100);
% close all;
% fc.visualizecoveragesimulation( ...
%     "contourLevelDetection", 0.85, ...
%     "saveFigure", 0, "header", 'scenario_1');