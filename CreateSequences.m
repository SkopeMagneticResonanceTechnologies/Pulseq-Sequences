%% Pulseq example sequences 
% including triggers and synchronization pre-scans for field-monitoring

% (c) 2024 Skope Magnetic Resonance Technologies AG

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

%% Create a 2D mono-polar dual-echo gradient-echo sequence
% Get default sequence parameters
paramsGre2d = SequenceParams('gre2d');

gre2d = skope_gre_2d(paramsGre2d);

% Plot sequence information after sync 
timeRange = [4.23 4.4];
gre2d.plot(timeRange);

% Test sequence
gre2d.test();

%% Create a 2D echo-planar imaging sequence
paramsEpi2d = SequenceParams('epi2d');
epi2d = skope_epi_2d(paramsEpi2d);

% Plot sequence information
timeRange = [5.5 5.536];
epi2d.plot(timeRange);

% Test sequence
epi2d.test();

%% Create off-resonance and position calibration sequence
paramsOpc = SequenceParams('opc');

opc = skope_offresAndPosCalib(paramsOpc);

% Plot sequence
timeRange = [0 5];
opc.plot(timeRange);

% Test sequence
opc.test();

%% Create local eddy current calibration sequence
paramsLec = SequenceParams('lec');

lec = skope_localEddyCalib(paramsLec);

% Plot sequence
timeRange = [0 22e-3];
lec.plot(timeRange);

% Test sequence
lec.test();

%% Create a series of blips
paramsGtf = SequenceParams('gtf');

gtf = skope_gtf(paramsGtf);

% Plot sequence information
timeRange = [0 300];
gtf.plot(timeRange);

% Test sequence
gtf.test();

%% Create two interleaved series of blips with half and double amplitude
paramsGtf = SequenceParams('gtf','linearityCheck');

gtf = skope_gtf(paramsGtf);

% Plot sequence information
timeRange = [0 300];
gtf.plot(timeRange);

% Test sequence
gtf.test();

%% Create a 3D monpolar dual-echo gradient-echo sequence
paramsGre3d = SequenceParams('gre3d');

gre3d = skope_gre_3d(paramsGre3d);

% Plot sequence information after sync 
timeRange = [0 100e-3];
gre3d.plot(timeRange);

% Test sequence
gre3d.test();
