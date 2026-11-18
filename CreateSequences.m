%% Pulseq example sequences 
% including triggers and synchronization pre-scans for field-monitoring

% (c) 2026 Skope Magnetic Resonance Technologies AG

%% Clean up
clear all
close all
clc

%% Check if Pulseq and SAFE PNS prediction modules have been added
if not(isfolder('pulseq/matlab'))
    error("Please run 'git submodule init' and 'git submodule update' to get the latest Pulseq scripts.")
end

if not(isfolder('safe_pns_prediction'))
    error("Please run 'git submodule init' and 'git submodule update' to get SAFE PNS Prediction scripts.")
end

%% Add Pulseq, sequences, methods and SAFE PNS
addpath('pulseq/matlab')
addpath(genpath('methods'))
addpath('sequences')
addpath('safe_pns_prediction')

%% Define scanner type
% ../methods/GetMRSystemSpecs.m
% 'Siemens 3T Cima.X', 'Siemens 3T Connectom', 'Siemens 7T Terra SC72CD', 'Siemens 9.4T SC72CD'
scannerType = 'Siemens 3T Cima.X';

%% Acquire more slices for in vivo scans and activate fat saturation
invivo = false;

%% Create the gradient heat sequence
paramsGradHeat = SequenceParams('grad_heat', scannerType);
paramsGradHeat.triggerOutput = 'ext1';

gradHeat = skope_grad_heat(paramsGradHeat);

% Plot the sequence and inspect the first 12 s to confirm the 10 s pause and trigger placement.
gradHeat.plot([0 30]);
gradHeat.test();
