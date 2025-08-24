%%% Anıl ARSLAN 2303980
clc; %clear; close all;
addpath(genpath([pwd '/library']));
set(groot, "defaultFigurePosition", [680 458 560 420]);

receivers = receivingNode( ...
    'position', [-physconst('LightSpeed')*3e-5; -5e3; 0], ...
    'CPIsecond', 0.25*1e-3, ...
    'samplingFrequency', 1e6);

transmitters = transmittingNode( ...
    'position', [physconst('LightSpeed')*3e-5; -5e3; 0]);
transmitters.setLFM( ...
    "bandWidth", 1e6, ...
    "frequencyDirection", "symmetric", ...
    "frequencyOffset", 0);

network = radarNetwork( ...
    'receivingNodes', receivers, ...
    'transmittingNodes', transmitters);

targets = target('position', [0; -5e3 + 1500; 0]);
int = interface( ...
    'network', network, ...
    'targets', targets);

int.visualizeellipses("ellipseType", "resolution");
posRX = receivers.position/1e3;
posTX = transmitters.position/1e3;
hold on; plot(posRX(1, :), posRX(2, :), '.b', 'LineWidth', 2, 'MarkerSize', 5);
plot(posTX(1, :), posTX(2, :), '.r', 'LineWidth', 2, 'MarkerSize', 5);
text(posRX(1, :), posRX(2, :), 'RX', 'Color', 'b', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(posTX(1, :), posTX(2, :), 0, 'TX', 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
xlim([-10 10]); ylim([-15 5]);
figureName = 'resolution'; fig2 = gcf;
savefig(fig2, ['C:\GitRepo\MSRS\figuresSim\' figureName '.fig']);
saveas(fig2, ['C:\GitRepo\MSRS\figuresSim\' figureName '.eps'], 'epsc');