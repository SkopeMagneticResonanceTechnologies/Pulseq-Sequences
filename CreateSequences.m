%% Pulseq example sequences 
% including triggers and synchronization pre-scans for field-monitoring

% (c) 2025 Skope Magnetic Resonance Technologies AG

%% Clean up
clear all
close all
clc

%% Check if Pulseq module has been added
if not(isfolder('pulseq/matlab'))
    error("Please run 'git submodule init' and 'git submodule update' to get the latest Pulseq scripts.")
end

%% Add Pulseq, sequences and methods
addpath('pulseq/matlab')
addpath('methods')
addpath('sequences')

%% Create a 3D monpolar dual-echo gradient-echo sequence
% The first 10 scans are sync scans
% After a pause of 4s, 50 dummy pulses are played out
% Thereafter, the actual scans are being played out
paramsGre3d = SequenceParams('gre3d');

paramsGre3d.Nz = paramsGre3d.Nx;
paramsGre3d.Nz = 4; % remove - only for quick testing

gre3d = skope_gre_3d(paramsGre3d);

% Plot sequence information after 5 sync scans and 4s pause
timeRange = [0 100e-3] + paramsGre3d.nDummy*paramsGre3d.TR + 5*(paramsGre3d.TR+1) + 4;
gre3d.plot(timeRange);

% Test sequence
gre3d.test();

%% Load Pulseq file (included in skope-i)
[grad, tx, rx, trig, labels, flags, defs] = mexSequenceSimulator('exports/skope_gre_3d_TRA_AP.seq', 'gradient');
%Labels: [SLC, SEG, REP, AVG, ECO, PHS, SET, LIN, PAR]
% Flags: [NAV, REV, SMS, PMC]

figure
stairs(labels), legend({'SLC', 'SEG', 'REP', 'AVG', 'ECO', 'PHS', 'SET', 'LIN', 'PAR'})

figure
stairs(flags), legend({'NAV', 'REV', 'SMS', 'PMC'})