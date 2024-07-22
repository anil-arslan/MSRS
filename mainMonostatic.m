%%% Anıl ARSLAN 2303980
clc; %clear; close all;
addpath(genpath([pwd '/library']));

receivers = receivingNode( ...
    'samplingFrequency', 2e7, ...
    'CPIsecond', 10e-6);
transmitters = transmittingNode( ...
    'inputPower_W', 1e3, ...
    'carrierFrequency', 1e9, ...
    'pulseWidth', 1e-6);
transmitters.setLFM( ...
    "bandWidth", 5e6, ...
    "frequencyDirection", "decreasing");
% transmitters.visualizewaveform("domain", "frequency");
network = radarNetwork( ...
    'receivingNodes', receivers, ...
    'transmittingNodes', transmitters);
% network.visualizenetwork;
% network.visualizewaveformsampled("axisAmbiguity", "3D");
targets = target( ...
    'position', [300; -300; 0], ...
    'meanRCS_dbms', -inf);
% targets = target( ...
%     'position', [d/4; -d/8; d/4], ...
%     'meanRCS_dbms', 1e3);
% targets.visualizetargets;
int = interface( ...
    'network', network, ...
    'targets', targets);
%%
sp = spu("interface", int, "gridResolution", [20 22 6000]);
sp.setintegrationindices;
% sp.visualizeintegrationindices;
%%
% transmitters.settaper("taperType", "rectwin");
int.configure('noise', 1, 'directPath', 0);
sp.configure("PFA", 1e-2);
%%
% int.visualizescenario;
int.visualizereceivedsignals;
sp.visualizefilteredsignals;
% sp.visualizehypothesizedtimedelays("dimension", "x-y-z");
% sp.visualizeintegrationindices;
%%
% sp.visualizeintegratedsignals("plot", "real");
sp.configure("processingAlgorithm", 1, "seedShuffle", 1, "numberOfTrials", 1, 'numberOfTrialsParallel', 1);
sp.visualizeintegratedsignals("plot", "magnitude");
% sp.visualizeestimation;
%%
% [PD, PFA] = sp.simulatedetection;