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

%% Create a 2D mono-polar dual-echo gradient-echo sequence for all directions
% Get default sequence parameters
paramsGre2d = SequenceParams('gre2d');

paramsGre2d.sliceOrientation = SliceOrientation.SAG;
paramsGre2d.phaseEncDir = PhaseEncodingDirection.AP;
gre2d = skope_gre_2d(paramsGre2d);

paramsGre2d.sliceOrientation = SliceOrientation.SAG;
paramsGre2d.phaseEncDir = PhaseEncodingDirection.HF;
gre2d = skope_gre_2d(paramsGre2d);

paramsGre2d.sliceOrientation = SliceOrientation.COR;
paramsGre2d.phaseEncDir = PhaseEncodingDirection.RL;
gre2d = skope_gre_2d(paramsGre2d);

paramsGre2d.sliceOrientation = SliceOrientation.COR;
paramsGre2d.phaseEncDir = PhaseEncodingDirection.HF;
gre2d = skope_gre_2d(paramsGre2d);

paramsGre2d.sliceOrientation = SliceOrientation.TRA;
paramsGre2d.phaseEncDir = PhaseEncodingDirection.RL;
gre2d = skope_gre_2d(paramsGre2d);

paramsGre2d.sliceOrientation = SliceOrientation.TRA;
paramsGre2d.phaseEncDir = PhaseEncodingDirection.AP;
gre2d = skope_gre_2d(paramsGre2d);

% Plot sequence information after sync 
timeRange = [4.25 4.27];
gre2d.plot(timeRange);

%% Test sequence
gre2d.test();

%% Create a 2D echo-planar imaging sequence for all directions
paramsEpi2d = SequenceParams('epi2d');

paramsEpi2d.sliceOrientation = SliceOrientation.SAG;
paramsEpi2d.phaseEncDir = PhaseEncodingDirection.AP;
epi2d = skope_epi_2d(paramsEpi2d);

paramsEpi2d.sliceOrientation = SliceOrientation.SAG;
paramsEpi2d.phaseEncDir = PhaseEncodingDirection.HF;
epi2d = skope_epi_2d(paramsEpi2d);

paramsEpi2d.sliceOrientation = SliceOrientation.COR;
paramsEpi2d.phaseEncDir = PhaseEncodingDirection.RL;
epi2d = skope_epi_2d(paramsEpi2d);

paramsEpi2d.sliceOrientation = SliceOrientation.COR;
paramsEpi2d.phaseEncDir = PhaseEncodingDirection.HF;
epi2d = skope_epi_2d(paramsEpi2d);

paramsEpi2d.sliceOrientation = SliceOrientation.TRA;
paramsEpi2d.phaseEncDir = PhaseEncodingDirection.RL;
epi2d = skope_epi_2d(paramsEpi2d);

paramsEpi2d.sliceOrientation = SliceOrientation.TRA;
paramsEpi2d.phaseEncDir = PhaseEncodingDirection.AP;
epi2d = skope_epi_2d(paramsEpi2d);

% Plot sequence information
timeRange = [5.999 6.070];
epi2d.plot(timeRange);

% Test sequence
epi2d.test();

%% EPI with acceleration factor 2
paramsEpi2d.sliceOrientation = SliceOrientation.TRA;
paramsEpi2d.phaseEncDir = PhaseEncodingDirection.AP;
paramsEpi2d.accFacPE = 2;
paramsEpi2d.TE  = 18.1e-3;
epi2d = skope_epi_2d(paramsEpi2d);

% Plot sequence information
timeRange = [5.999 6.070];
epi2d.plot(timeRange);

% Test sequence
epi2d.test();

%% EPI with acceleration factor 2
paramsEpi2d.sliceOrientation = SliceOrientation.COR;
paramsEpi2d.phaseEncDir = PhaseEncodingDirection.RL;
paramsEpi2d.accFacPE = 2;
paramsEpi2d.TE  = 25.1e-3;

paramsEpi2d.Nx = 120;
paramsEpi2d.Ny = 120;
paramsEpi2d.readoutTime = 650e-6;
epi2d = skope_epi_2d(paramsEpi2d);

% Plot sequence information
timeRange = [5.999 6.070];
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

% Run as well with one average for nominal gradient simulation
paramsGtf.nAve = 1; 
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

%% Create off-resonance and position calibration sequence
load('./waveforms/sweepWaveform.mat')
paramsSweep = SequenceParams('sweep');
sweep = skope_sweep(paramsSweep,sweepWaveform);

% Run as well with one average for nominal gradient simulation
paramsSweep.nAve = 1;   
sweep = skope_sweep(paramsSweep,sweepWaveform);

% Plot sequence
timeRange = [0 8];
sweep.plot(timeRange);

% Test sequence
sweep.test();

%% Create a 2D spiral gradient-echo sequence
% Get default sequence parameters
paramsSpiral2d = SequenceParams('spiral2d');
load('./waveforms/spiralGrad_FOV192_RES1mm_minRise6_maxAmp40_nitlv16.mat'); % [Hz/m]
spiral2d = skope_spiral_2d(paramsSpiral2d,spiralWaveform);

% Plot sequence information after sync 
timeRange = [5.4 5.42];
spiral2d.plot(timeRange);

% Test sequence
spiral2d.test();