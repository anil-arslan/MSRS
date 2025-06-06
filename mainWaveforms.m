%%% Anıl ARSLAN 2303980
clc; %clear; close all;
addpath(genpath([pwd '/library']));

d = 1e3;
receivers = receivingNode( ...
    'samplingFrequency', 4e7, ...
    'CPIsecond', 20e-6);
transmitters = transmittingNode( ...
    'position', [-d 0; 0 -d; 0 0], ...
    'inputPower', 1e3, ...
    'carrierFrequency', 1e9, ...
    'pulseWidth', 1e-6);
network = radarNetwork( ...
    'receivingNodes', receivers, ...
    'transmittingNodes', transmitters);
transmitters(1).setunmodulation;
transmitters(2).setLFM( ...
    "bandWidth", 5e6, ...
    "frequencyDirection", "decreasing", ...
    "frequencyOffset", 0);
transmitters.visualizewaveform( ...
    "domain", "time", ...
    "plot", "real");
transmitters.visualizewaveform( ...
    "domain", "frequency", ...
    "plot", "magnitude");
network.visualizewaveformsampled( ...
    "domain", "ambiguity", ...
    "axisAmbiguity", "3D");
network.visualizewaveformsampled( ...
    "domain", "frequency", ...
    "plot", "magnitude");