%% Pulseq example sequences 
% including triggers and synchronization pre-scans for field-monitoring

% (c) 2026 Skope Magnetic Resonance Technologies AG

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
addpath(genpath('methods'))
addpath('sequences')

%% Define scanner type
% ../methods/GetMRSystemSpecs.m
% 'Siemens 3T Cima.X', 'Siemens 3T Connectom', 'Siemens 7T Terra SC72CD', 'Siemens 9.4T SC72CD'
scannerType = 'Siemens 3T Cima.X';

%% Acquire more slices for in vivo scans and activate fat saturation
invivo = true;

%% Create a 2D mono-polar dual-echo gradient-echo (GRE) sequence 
% Get default sequence parameters
paramsGre2d = SequenceParams('gre2d',scannerType);
paramsGre2d.fov = 0.22;
paramsGre2d.Nx = 80;
paramsGre2d.Ny = 80;
paramsGre2d.thickness = 3e-3;

paramsGre2d.nSlices = 12;
paramsGre2d.distanceFactorPercentage = 150;
if invivo
    paramsGre2d.nSlices = 44;
    paramsGre2d.distanceFactorPercentage = 10;
    paramsGre2d.seqSpecName = 'invivo';
    % GRE has currently no fat saturation pulse 
end

% Generate the sequence
gre2d = skope_gre_2d(paramsGre2d);

% Plot sequence diagram for the first 10 s
timeRange = [0 10];
gre2d.plot(timeRange);

% Test sequence
gre2d.test();

%% Create a 2D echo-planar imaging (EPI) sequence
paramsEpi2d = SequenceParams('epi2d',scannerType);
paramsEpi2d.fov = 0.22;
paramsEpi2d.Nx = 74;
paramsEpi2d.Ny = 74;
paramsEpi2d.TE = 25e-3;
paramsEpi2d.TR = 60e-3;
paramsEpi2d.readoutTime = 500e-6; 
paramsEpi2d.thickness = 3e-3;

paramsEpi2d.nSlices = 12;
paramsEpi2d.distanceFactorPercentage = 150; 
if invivo
    paramsEpi2d.nSlices = 44;
    paramsEpi2d.distanceFactorPercentage = 10;
    paramsEpi2d.seqSpecName = 'invivo';
    paramsEpi2d.doPlayFatSat = true;
end

% Generate the sequence
epi2d = skope_epi_2d(paramsEpi2d);

% Plot sequence diagram for the first 10 s
timeRange = [0 10];
epi2d.plot(timeRange);

% Test sequence
epi2d.test();

%% Create a 2D echo-planar imaging (EPI) sequence with 2-fold in-plane acceleration
paramsEpi2d = SequenceParams('epi2d',scannerType);
paramsEpi2d.fov = 0.22;
paramsEpi2d.accFacPE = 2;
paramsEpi2d.Nx = 120;
paramsEpi2d.Ny = 120;
paramsEpi2d.TE = 25e-3;
paramsEpi2d.TR = 60e-3;
paramsEpi2d.readoutTime = 500e-6; 
paramsEpi2d.thickness = 3e-3;

paramsEpi2d.nSlices = 12;
paramsEpi2d.distanceFactorPercentage = 150; 
if invivo
    paramsEpi2d.nSlices = 44;
    paramsEpi2d.distanceFactorPercentage = 10;
    paramsEpi2d.seqSpecName = 'invivo';
    paramsEpi2d.doPlayFatSat = true;
end

% Generate the sequence
epi2d = skope_epi_2d(paramsEpi2d);

% Plot sequence diagram for the first 10 s
timeRange = [0 10];
epi2d.plot(timeRange);

% Test sequence
epi2d.test();

%% Create a 2D echo-planar imaging (EPI) sequence with no acceleration and 2-fold multi-band excitation
paramsEpi2d = SequenceParams('epi2d',scannerType);
paramsEpi2d.fov = 0.22;
paramsEpi2d.accFacPE = 1;
paramsEpi2d.multiBandFactor = 2;
paramsEpi2d.Nx = 74;
paramsEpi2d.Ny = 74;
paramsEpi2d.TE = 30e-3;
paramsEpi2d.TR = 80e-3;
paramsEpi2d.readoutTime = 500e-6; 
paramsEpi2d.thickness = 3e-3;

paramsEpi2d.nSlices = 12;
paramsEpi2d.distanceFactorPercentage = 150; 
if invivo
    paramsEpi2d.nSlices = 44;
    paramsEpi2d.distanceFactorPercentage = 10;
    paramsEpi2d.seqSpecName = 'invivo';
    paramsEpi2d.doPlayFatSat = true;
end

% Generate the sequence
epi2d = skope_epi_2d(paramsEpi2d);

% Plot sequence diagram for the first 10 s
timeRange = [0 10];
epi2d.plot(timeRange);

% Test sequence
epi2d.test();

%% Create a 2D echo-planar imaging (EPI) sequence with no acceleration and 4-fold multi-band excitation
paramsEpi2d = SequenceParams('epi2d',scannerType);
paramsEpi2d.fov = 0.22;
paramsEpi2d.accFacPE = 1;
paramsEpi2d.multiBandFactor = 4;
paramsEpi2d.Nx = 74;
paramsEpi2d.Ny = 74;
paramsEpi2d.TE = 35e-3;
paramsEpi2d.TR = 80e-3;
paramsEpi2d.readoutTime = 500e-6; 
paramsEpi2d.thickness = 3e-3;

paramsEpi2d.nSlices = 12;
paramsEpi2d.distanceFactorPercentage = 150; 
if invivo
    paramsEpi2d.nSlices = 44;
    paramsEpi2d.distanceFactorPercentage = 10;
    paramsEpi2d.seqSpecName = 'invivo';
    paramsEpi2d.doPlayFatSat = true;
end

% Generate the sequence
epi2d = skope_epi_2d(paramsEpi2d);

% Plot sequence diagram for the first 10 s
timeRange = [0 10];
epi2d.plot(timeRange);

% Test sequence
epi2d.test();

%% Create a 2D echo-planar imaging (EPI) sequence with 2-fold in-plane acceleration and 2-fold multi-band excitation
paramsEpi2d = SequenceParams('epi2d',scannerType);
paramsEpi2d.fov = 0.22;
paramsEpi2d.accFacPE = 2;
paramsEpi2d.multiBandFactor = 2;
paramsEpi2d.Nx = 120;
paramsEpi2d.Ny = 120;
paramsEpi2d.TE = 25e-3;
paramsEpi2d.TR = 60e-3;
paramsEpi2d.readoutTime = 500e-6; 
paramsEpi2d.thickness = 3e-3;

paramsEpi2d.nSlices = 12;
paramsEpi2d.distanceFactorPercentage = 150; 
if invivo
    paramsEpi2d.nSlices = 44;
    paramsEpi2d.distanceFactorPercentage = 10;
    paramsEpi2d.seqSpecName = 'invivo';
    paramsEpi2d.doPlayFatSat = true;
end

% Generate the sequence
epi2d = skope_epi_2d(paramsEpi2d);

% Plot sequence diagram for the first 10 s
timeRange = [0 10];
epi2d.plot(timeRange);

% Test sequence
epi2d.test();

%% EPI with acceleration factor 3 and higher resolution
paramsEpi2d = SequenceParams('epi2d',scannerType);
paramsEpi2d.fov = 0.22;
paramsEpi2d.accFacPE = 3;
paramsEpi2d.Nx = 130;
paramsEpi2d.Ny = 130;
paramsEpi2d.TE = 25e-3;
paramsEpi2d.TR = 60e-3;
paramsEpi2d.readoutTime = 500e-6; 
paramsEpi2d.thickness = 3e-3;

paramsEpi2d.nSlices = 12;
paramsEpi2d.distanceFactorPercentage = 150; 
if invivo
    paramsEpi2d.nSlices = 44;
    paramsEpi2d.distanceFactorPercentage = 10;
    paramsEpi2d.seqSpecName = 'invivo';
    paramsEpi2d.doPlayFatSat = true;
end

% Generate the sequence
epi2d = skope_epi_2d(paramsEpi2d);

% Plot sequence diagram for the first 10 s
timeRange = [0 10]; 
epi2d.plot(timeRange);

% Test sequence
epi2d.test();

%% Create a 2D spin-echo EPI sequence with diffusion encoding
% navigator is by default disabled here
paramsSeEpi2dDiff = SequenceParams('se_epi2d_diff',scannerType);
paramsSeEpi2dDiff.fov = 0.22;
paramsSeEpi2dDiff.accFacPE = 3;
paramsSeEpi2dDiff.Nx = 128;
paramsSeEpi2dDiff.Ny = 128;
% Reduce slew rate (To be corrected by independent slew rate for diffusion gradients)
paramsSeEpi2dDiff.maxSlew = 120; 

paramsSeEpi2dDiff.nSlices = 12;
paramsSeEpi2dDiff.distanceFactorPercentage = 150; 
if invivo
    paramsSeEpi2dDiff.nSlices = 44;
    paramsSeEpi2dDiff.distanceFactorPercentage = 10;
    paramsSeEpi2dDiff.seqSpecName = 'invivo';
    paramsSeEpi2dDiff.doPlayFatSat = true;
end

% Generate the sequence
seepi2d = skope_se_epi_2d_diff(paramsSeEpi2dDiff);

% Plot sequence diagram for the first 5 to 20 s
timeRange = [5 20];
seepi2d.plot(timeRange);

% Test sequence
seepi2d.test();

%% Create a 3D mono-polar dual-echo gradient-echo sequence
paramsGre3d = SequenceParams('gre3d',scannerType);

% Generate the sequence
gre3d = skope_gre_3d(paramsGre3d);

% Plot sequence information after sync 
timeRange = [0 100e-3];
gre3d.plot(timeRange);

% Test sequence
gre3d.test();

%% Create off-resonance and position calibration sequence for all possible trigger output channels
paramsOpc = SequenceParams('opc',scannerType);

paramsOpc.triggerOutput = 'ext1'; % Default optical output
opc = skope_offresAndPosCalib(paramsOpc);

paramsOpc.triggerOutput = 'osc0';
opc = skope_offresAndPosCalib(paramsOpc);

% Plot sequence
timeRange = [0 5];
opc.plot(timeRange);

% Test sequence
opc.test();

%% Create local eddy current calibration sequence
paramsLec = SequenceParams('lec',scannerType);

% Generate the sequence
lec = skope_localEddyCalib(paramsLec);

% Plot sequence
timeRange = [0 22e-3];
lec.plot(timeRange);

% Test sequence
lec.test();

%% Create a series of blips
paramsGtf = SequenceParams('gtf',scannerType);

% paramsGtf.nAve = 1; 
gtf = skope_gtf(paramsGtf);

% Plot sequence information
timeRange = [0 300];
gtf.plot(timeRange);

% Test sequence
gtf.test();

%% Create two interleaved series of blips with half and double amplitude
paramsGtf = SequenceParams('gtf',scannerType,'linearityCheck');
gtf = skope_gtf(paramsGtf);

% Plot sequence information
timeRange = [0 300];
gtf.plot(timeRange);

% Test sequence
gtf.test();

%% Create off-resonance and position calibration sequence
load('./waveforms/sweepWaveform.mat')
paramsSweep = SequenceParams('sweep',scannerType);

%paramsSweep.nAve = 1;   
sweep = skope_sweep(paramsSweep,sweepWaveform);

% Plot sequence
timeRange = [0 8];
sweep.plot(timeRange);

% Test sequence
sweep.test();