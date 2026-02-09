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
addpath('methods')
addpath('sequences')

%% Define scanner type
% 'Siemens 3T Cima.X', 'Siemens 7T Terra SC72CD', 'Siemens 9.4T SC72CD'
scannerType = 'Siemens 3T Cima.X';

coeffsX.Amp = [-0.0014    0.0752    0.0071];
coeffsX.Tau = [1.2498    0.1236    0.0106];
coeffsY.Amp = [-0.0262    0.1384    0.0180];
coeffsY.Tau = [0.1943    0.1190    0.0026];
coeffsZ.Amp = [0.0249    0.2048    0.0455];
coeffsZ.Tau = [0.1872    0.1231    0.0445];

%% Create a 2D mono-polar dual-echo gradient-echo (GRE) sequence 
% Get default sequence parameters
paramsGre2d = SequenceParams('gre2d',scannerType);

% Set slice orientation and encoding direction
paramsGre2d.sliceOrientation = SliceOrientation.TRA;
paramsGre2d.phaseEncDir = PhaseEncodingDirection.AP;

paramsGre2d.seqSpecName = 'os2';
paramsGre2d.nDummy = 0;
% Generate the sequence
gre2d = skope_gre_2d(paramsGre2d);

% Plot first 10 s
timeRange = [0 10];
gre2d.plot(timeRange);

% Plot sequence information after sync 
timeRange = [4.25 4.27]; %s
gre2d.plot(timeRange);

% Test sequence
gre2d.test();

%% Simulate ECC phase
phase = mexSequenceSimulator('exports/Siemens 3T Cima.X/skope_gre_2d_TRA_AP_os2.seq', 'eddyPhase',coeffsX, coeffsY, coeffsZ);
plot(phase)

%% Create a 2D echo-planar imaging (EPI) sequence
paramsEpi2d = SequenceParams('epi2d',scannerType);
paramsEpi2d.sliceOrientation = SliceOrientation.TRA;
paramsEpi2d.phaseEncDir = PhaseEncodingDirection.AP;
paramsEpi2d.accFacPE = 1;
paramsEpi2d.nRep = 1;
paramsEpi2d.TE = 33.3e-3;
paramsEpi2d.TR = 130e-3;
paramsEpi2d.Nx = 100; 
paramsEpi2d.Ny = 100;
switch scannerType
    case 'Siemens 3T Cima.X'
        paramsEpi2d.readoutTime = 520e-6;
    otherwise
        paramsEpi2d.readoutTime = 800e-6;
end

paramsEpi2d.seqSpecName = 'os2';
paramsEpi2d.alpha = 7;
paramsEpi2d.nDummy = 0;

epi2d = skope_epi_2d(paramsEpi2d);

% Plot sequence information
timeRange = [0 10];
epi2d.plot(timeRange);

% Test sequence
epi2d.test();

%% Simulate ECC phase
phase = mexSequenceSimulator('exports/Siemens 3T Cima.X/skope_epi_2d_TRA_AP_R1_os2.seq', 'eddyPhase',coeffsX, coeffsY, coeffsZ);
plot(phase)

%% Create a 2D spin-echo EPI sequence with diffusion encoding
% navigator is by default disabled here
paramsSeEpi2dDiff = SequenceParams('se_epi2d_diff',scannerType);
paramsSeEpi2dDiff.sliceOrientation = SliceOrientation.TRA;
paramsSeEpi2dDiff.phaseEncDir = PhaseEncodingDirection.AP;
paramsSeEpi2dDiff.accFacPE = 1;
paramsSeEpi2dDiff.nRep = 1;
paramsSeEpi2dDiff.TE = 100e-3;
paramsSeEpi2dDiff.TR = 140e-3;
paramsSeEpi2dDiff.Nx = 100;
paramsSeEpi2dDiff.Ny = 100;
paramsSeEpi2dDiff.maxSlew = 120;

paramsSeEpi2dDiff.seqSpecName = 'os2';

switch scannerType
    case 'Siemens 3T Cima.X'
        paramsSeEpi2dDiff.readoutTime = 600e-6;
    otherwise
        paramsSeEpi2dDiff.readoutTime = 680e-6;
end
paramsSeEpi2dDiff.nSlices = 15; %testing
paramsSeEpi2dDiff.nDummy = 0; %testing
seepi2d = skope_se_epi_2d_diff(paramsSeEpi2dDiff);

% Plot sequence information
timeRange = [0 10];
seepi2d.plot(timeRange);

% Test sequence
seepi2d.test();

%% Simulate ECC phase
phase = mexSequenceSimulator('exports/Siemens 3T Cima.X/skope_se_epi_2d_diff_TRA_AP_os2.seq', 'eddyPhase',coeffsX, coeffsY, coeffsZ);
plot(phase,'*')
