%% Pulseq example sequences tested on Siemens 3T Cima.X
% for phantom
% GRE; EPI, accelerated EPI high res
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

%% 2D GRE
paramsGre2d = SequenceParams('gre2d');

paramsGre2d.sliceOrientation = SliceOrientation.TRA;
paramsGre2d.phaseEncDir = PhaseEncodingDirection.AP;
paramsGre2d.scannerType = 'Siemens 3T Cima.X';
paramsGre2d.maxGrad = 190;
paramsGre2d.maxSlew = 190;

% paramsGre2d.doMonitoringDuringRF = 1;
% paramsGre2d.seqSpecName = 'FMduringRF';   

% Generate the sequence
gre2d = skope_gre_2d(paramsGre2d);

% Plot sequence information after sync 
timeRange = [4.25 4.27]; %s
gre2d.plot(timeRange);

% Test sequence
gre2d.test();

%% 2D GRE with monitoring during RFs
paramsGre2d = SequenceParams('gre2d');

paramsGre2d.sliceOrientation = SliceOrientation.TRA;
paramsGre2d.phaseEncDir = PhaseEncodingDirection.AP;
paramsGre2d.scannerType = 'Siemens 3T Cima.X';
paramsGre2d.maxGrad = 190;
paramsGre2d.maxSlew = 190;
                    
% Generate the sequence
gre2d = skope_gre_2d(paramsGre2d);

% Plot sequence information after sync 
timeRange = [4.25 4.27]; %s
gre2d.plot(timeRange);

% Test sequence
gre2d.test();

%% 2D EPI
paramsEpi2d = SequenceParams('epi2d');
paramsEpi2d.sliceOrientation = SliceOrientation.TRA;
paramsEpi2d.phaseEncDir = PhaseEncodingDirection.AP;
paramsEpi2d.Nx = 100; 
paramsEpi2d.Ny = 100;
paramsEpi2d.TE = 31e-3;
paramsEpi2d.TR = 130e-3;
paramsEpi2d.readoutTime = 0.500e-3;
paramsEpi2d.nRep = 1;
paramsEpi2d.accFacPE = 1;
paramsEpi2d.scannerType = 'Siemens 3T Cima.X';
paramsEpi2d.maxGrad = 190;
paramsEpi2d.maxSlew = 190;

% paramsEpi2d.nDummy = 1; %testing

paramsEpi2d.seqSpecName = ['ROtime_' num2str(paramsEpi2d.readoutTime*1e6)];

epi2d = skope_epi_2d(paramsEpi2d);

% Test sequence
epi2d.test();

%% high resolution 2D EPI with acceleration factor R=4
paramsEpi2d = SequenceParams('epi2d');
paramsEpi2d.sliceOrientation = SliceOrientation.TRA;
paramsEpi2d.phaseEncDir = PhaseEncodingDirection.AP;
paramsEpi2d.Nx = 200; 
paramsEpi2d.Ny = 200;
paramsEpi2d.TE = 25e-3;
paramsEpi2d.TR = 130e-3;
paramsEpi2d.readoutTime = 0.800e-3;
paramsEpi2d.nRep = 1;
paramsEpi2d.accFacPE = 4;
paramsEpi2d.scannerType = 'Siemens 3T Cima.X';
paramsEpi2d.maxGrad = 190;
paramsEpi2d.maxSlew = 190;

% paramsEpi2d.nDummy = 1; %testing

paramsEpi2d.seqSpecName = ['ROtime_' num2str(paramsEpi2d.readoutTime*1e6)];

epi2d = skope_epi_2d(paramsEpi2d);
%%
% Test sequence
epi2d.test();
