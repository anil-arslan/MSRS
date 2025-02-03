%%% Anıl ARSLAN 2303980
clc; %clear; close all;
addpath(genpath([pwd '/library']));

arrayRX = planarArray("numberOfElements", [11 11]);
arrayTX = planarArray("numberOfElements", [5 5]);
arrayRX.steer("steeringAzimuth", 30);
d = 1e3;
posRX = [
    d -d;
    0 0;
    0 0];
receivers = receivingNode( ...
    'position', posRX, ...
    'array', arrayRX, ...
    'samplingFrequency', 2e7, ...
    'CPIsecond', 10e-6);
transmitters = transmittingNode( ...
    'position', [0; -d; 0], ...
    'array', arrayTX, ...
    'inputPower_W', 1e3, ...
    'carrierFrequency', 1e9, ...
    'pulseWidth', 1e-6);
transmitters.setLFM( ...
    "bandWidth", 5e6, ...
    "frequencyDirection", "decreasing");
network = radarNetwork( ...
    'receivingNodes', receivers, ...
    'transmittingNodes', transmitters);
targets = target( ...
    ... 'position', [1.03770999999961; 1.03770999999961; 0], ...
    'position', [ ...
    1.03770999999961 0; ...
    -1389.6 0; ...
    0 0], ...
    'meanRCS_dbms', -25);
int = interface( ...
    'network', network, ...
    'targets', targets);
% arrayRX.visualizearray;
% arrayTX.visualizearray;
% network.visualizenetwork;
int.visualizescenario("showPattern", 1);
%%
sp = spu("interface", int, "gridResolution", [20 20 6000]);
sp.setintegrationindices;
%%
int.configure( ...
    'noise', 1, ...
    'directPath', 0, ...
    'spatialCoherency', 'deterministic');
sp.configure( ...
    "PFA", 1e-4, ...
    "processingAlgorithm", 6, ...
    "seedShuffle", 1, ...
    "numberOfTrials", 200, ...
    'numberOfTrialsParallel', 200);
%%
int.visualizereceivedsignals;
sp.visualizefilteredsignals;
%%
% sp.setmatchfilteredsignals;
sp.visualizeintegratedsignals("trialID", 31); % 17 20
% sp.visualizeestimation;
%%
% [PD, PFA] = sp.simulatedetection;
%%
sp.visualizereceiveroperatingcharacteristics;
sp.visualizedetectioncharacteristics("PFA", 1e-4);