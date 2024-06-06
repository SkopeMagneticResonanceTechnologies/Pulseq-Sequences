% (c) 2022 Skope Magnetic Resonance Technologies AG

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

%% Select scanner
scannerType = "Siemens Terra 7T SC72CD";

%% Create a series of blips
gtf = skope_gtf(scannerType);

% Plot sequence information after sync 
timeRange = [0 18*3*20e-3];
gtf.plot(timeRange);

% Test sequence
gtf.test();

%% Create a 2D monpolar dual-echo gradient-echo sequence
gre2d = skope_gre_2d(scannerType);

% Plot sequence information after sync 
timeRange = [4.23 4.4];
gre2d.plot(timeRange);

% Test sequence
gre2d.test();

%% Create a 3D monpolar dual-echo gradient-echo sequence
close all
gre3d = skope_gre_3d(scannerType);

% Plot sequence information after sync 
timeRange = [0 100e-3];
gre3d.plot(timeRange);

% Test sequence
gre3d.test();

%% Create a 2D echo-planar imaging sequence
epi = skope_epi_2d(scannerType);

% Plot sequence information
timeRange = [5.410 5.480];
epi.plot(timeRange);

% Test sequence
epi.test();

%% Create off-resonance and position calibration sequence
opc = skope_offresAndPosCalib(scannerType);

% Plot sequence
timeRange = [0 5];
opc.plot(timeRange);

% Test sequence
opc.test();

%% Create local eddy current calibration sequence
lec = skope_localEddyCalib(scannerType);

% Plot sequence
timeRange = [0 22e-3];
lec.plot(timeRange);

% Test sequence
lec.test();
