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
    'position', [-318.96229; -310.96229; 0], ...
    'meanRCS_dbms', -10);
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
int.configure( ...
    'noise', 1, ...
    'spatialCoherency', 'incoherent');
sp.configure( ...
    "PFA", 1e-3, ...
    "processingAlgorithm", 6, ...
    "seedShuffle", 1, ...
    "numberOfTrials", 200, ...
    'numberOfTrialsParallel', 200);
%%
% int.visualizescenario;
int.visualizereceivedsignals;
sp.visualizefilteredsignals;
% sp.visualizehypothesizedtimedelays("dimension", "x-y-z");
% sp.visualizeintegrationindices;
%%
% sp.visualizeintegratedsignals("plot", "real");
% sp.setmatchfilteredsignals;
sp.visualizeintegratedsignals;
% sp.visualizeestimation;
%%
[PD, PFA] = sp.simulatedetection;
%%
% sp.visualizereceiveroperatingcharacteristics;
% sp.visualizedetectioncharacteristics;
% sp.visualizedetectioncharacteristics("PFA", 1e-3);
% sp.visualizereceiveroperatingcharacteristics("scanAlgorithms", 1, "snr_dB", 10);