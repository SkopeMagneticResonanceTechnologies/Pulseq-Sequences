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

%% Create a 2D mono-polar dual-echo gradient-echo (GRE) sequence 
% Get default sequence parameters
paramsGre2d = SequenceParams('gre2d',scannerType);

% Set slice orientation and encoding direction
paramsGre2d.sliceOrientation = SliceOrientation.TRA; %TRA
paramsGre2d.phaseEncDir = PhaseEncodingDirection.AP; %AP
% Generate the sequence
gre2d = skope_gre_2d(paramsGre2d);

% % Plot first 10 s
% timeRange = [0 10];
% gre2d.plot(timeRange);
% 
% % Plot sequence information after sync 
% timeRange = [4.25 4.27]; %s
% gre2d.plot(timeRange);

% Test sequence
% gre2d.test();

%% Create a 2D echo-planar imaging (EPI) sequence
paramsEpi2d = SequenceParams('epi2d',scannerType);
paramsEpi2d.sliceOrientation = SliceOrientation.TRA;
paramsEpi2d.phaseEncDir = PhaseEncodingDirection.AP;
paramsEpi2d.accFacPE = 2;
paramsEpi2d.nRep = 1;
paramsEpi2d.TE = 39e-3;
paramsEpi2d.TR = 130e-3;
paramsEpi2d.Nx = 96; 
paramsEpi2d.Ny = 96;
paramsEpi2d.multiBandFactor = 2;
switch scannerType
    case 'Siemens 3T Cima.X'
        paramsEpi2d.readoutTime = 500e-6;
    otherwise
        paramsEpi2d.readoutTime = 800e-6;
end

clc
% paramsEpi2d.ro_os = 1;
paramsEpi2d.nSlices = 10; %testing=1, otherwise = 15
paramsEpi2d.nDummy = 1; %testing=1, otherwise = 10
paramsEpi2d.seqSpecName = ['os' num2str(paramsEpi2d.ro_os) ''];
epi2d = skope_epi_2d(paramsEpi2d);

% Plot sequence information
% timeRange = [25.05 25.13];
% epi2d.plot([5.3 5.5]);
% epi2d.plot(timeRange);
epi2d.plot([0 10]);
% Test sequence
epi2d.test();

%% EPI with acceleration factor 3 and higher resolution
paramsEpi2d = SequenceParams('epi2d',scannerType);
paramsEpi2d.sliceOrientation = SliceOrientation.TRA;
paramsEpi2d.phaseEncDir = PhaseEncodingDirection.AP;
paramsEpi2d.accFacPE = 3;
paramsEpi2d.nRep = 1;
paramsEpi2d.TE = 26e-3;
paramsEpi2d.TR = 130e-3;
paramsEpi2d.Nx = 130;
paramsEpi2d.Ny = 130;
paramsEpi2d.maxSlew = 170;
switch scannerType
    case 'Siemens 3T Cima.X'
        paramsEpi2d.readoutTime = 800e-6;
    otherwise
        paramsEpi2d.readoutTime = 680e-6;
end

paramsEpi2d.nSlices = 15; %testing=1, otherwise = 15
paramsEpi2d.nDummy = 10; %testing=1, otherwise = 10
paramsEpi2d.seqSpecName = ['os' num2str(paramsEpi2d.ro_os)];
epi2d = skope_epi_2d(paramsEpi2d);

% Plot sequence information
% timeRange = [5 6];
% timeRange = [26.08 26.18];
% epi2d.plot(timeRange);

% Test sequence
epi2d.test();

%% Create a 2D spin-echo EPI sequence with diffusion encoding
% navigator is by default disabled here
paramsSeEpi2dDiff = SequenceParams('se_epi2d_diff',scannerType);
paramsSeEpi2dDiff.sliceOrientation = SliceOrientation.TRA;
paramsSeEpi2dDiff.phaseEncDir = PhaseEncodingDirection.AP;
paramsSeEpi2dDiff.accFacPE = 3;
paramsSeEpi2dDiff.nRep = 1;
paramsSeEpi2dDiff.TE = 60e-3;
paramsSeEpi2dDiff.TR = 130e-3;
paramsSeEpi2dDiff.Nx = 128;
paramsSeEpi2dDiff.Ny = 128;
paramsSeEpi2dDiff.maxSlew = 120; %to be corrected by independent slew rate of diff gradients
switch scannerType
    case 'Siemens 3T Cima.X'
        paramsSeEpi2dDiff.readoutTime = 800e-6;
    otherwise
        paramsSeEpi2dDiff.readoutTime = 680e-6;
end

paramsSeEpi2dDiff.nSlices = 15; %testing=1, otherwise = 15
paramsSeEpi2dDiff.nDummy = 3; %testing=1, otherwise = 3
paramsSeEpi2dDiff.seqSpecName = ['os' num2str(paramsSeEpi2dDiff.ro_os)];
seepi2d = skope_se_epi_2d_diff(paramsSeEpi2dDiff);

% Plot sequence information
timeRange = [5 20];
seepi2d.plot(timeRange);

% Test sequence
% seepi2d.test();

%% Create a 3D monpolar dual-echo gradient-echo sequence
paramsGre3d = SequenceParams('gre3d',scannerType);
gre3d = skope_gre_3d(paramsGre3d);

% Plot sequence information after sync 
timeRange = [0 100e-3];
gre3d.plot(timeRange);

% Test sequence
gre3d.test();

%% Create off-resonance and position calibration sequence for all possible trigger output channels

paramsOpc = SequenceParams('opc',scannerType);
paramsOpc.triggerOutput = 'ext1'; % Default optical ouput
opc = skope_offresAndPosCalib(paramsOpc);

paramsOpc.triggerOutput = 'osc0';
opc = skope_offresAndPosCalib(paramsOpc);

paramsOpc.triggerOutput = 'osc1';
opc = skope_offresAndPosCalib(paramsOpc);

% Plot sequence
timeRange = [0 5];
opc.plot(timeRange);

% Test sequence
opc.test();

%% Create local eddy current calibration sequence
paramsLec = SequenceParams('lec',scannerType);
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

% paramsSweep.nAve = 1;   
sweep = skope_sweep(paramsSweep,sweepWaveform);

% Plot sequence
timeRange = [0 8];
sweep.plot(timeRange);

% Test sequence
sweep.test();