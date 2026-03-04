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

%% Create a multi-shot 2D spiral gradient-echo sequence with pre-generated spiral waveform
% Get default sequence parameters
paramsSpiral2d = SequenceParams('spiral2d',scannerType);
load('./waveforms/spiralGrad_FOV192_RES1mm_minRise6_maxAmp40_nitlv16.mat'); % [Hz/m]
obj.sys.gamma = 42576000;
spiralWaveform = spiralWaveform / obj.sys.gamma * 1000;
paramsSpiral2d.seqSpecName = 'hardcoded';

spiral2d = skope_spiral_2d(paramsSpiral2d,spiralWaveform);

% Plot sequence information after sync 
timeRange = [5.4 5.42];
spiral2d.plot(timeRange);

% Test sequence
spiral2d.test();

%% Create a multi-shot 2D spiral gradient-echo sequence with specific spiral waveform
minTimeGradientDir = fullfile(pwd, 'minTimeGradient', 'Matlab');
if not(isfolder(minTimeGradientDir))
    error('Download https://people.eecs.berkeley.edu/~mlustig/software/tOptGrad_V0.2.tar.gz')
end
addpath(genpath(minTimeGradientDir))

%-------------------------------------------------------------------------
% Compute trajectory
%-------------------------------------------------------------------------
Nitlv = 16;             % Number of interleves
r = 0;                  % rv/riv Indicates type of solution
res	= 1;                % Resolution (in mm)
fov	= [25,25];          % Vector of fov (in cm)
radius = [0,1];         % Vector of radius corresponding to the fov
Gmax = 4;               % Max gradient (default 3 G/CM = 30 mT/m)
Smax = 15;              % Max slew (default 10 G/cm/ms = 100 mT/m/ms)
T = 10e-3;              % Sampling rate (in ms) - 10 us on Siemens systems
ds = [];                % Step size for integration
interpType = 'linear';   % Type of interpolation used to interpolate the fov accept: linear, cubic, spline

[k_rv,g_rv,s_rv,time_rv,Ck_rv] = vdSpiralDesign(Nitlv, r, res,fov,radius,Gmax,Smax,T,ds,interpType);

% Convert gradient to mT/m
g_rv = [0,0; g_rv(:,1:2); 0,0] * 10;
plot(g_rv)
plot(g_rv(:,1), g_rv(:,2))
% size(g_rv,1)
%-------------------------------------------------------------------------


%%
paramsSpiral2d = SequenceParams('spiral2d',scannerType);
paramsSpiral2d.Ny = Nitlv;
paramsSpiral2d.Nx = 248; %from fov/res
paramsSpiral2d.fov = 250e-3; %from fov(1)
paramsSpiral2d.readoutTime = 4.1e-3; %from size(k_rv,1) * 1e-3
paramsSpiral2d.mode = 'multiShot';

% Create sequence
spiral2d = skope_spiral_2d(paramsSpiral2d,g_rv);

% Plot sequence information after sync 
timeRange = [5.4 5.42];
spiral2d.plot(timeRange);

% Test sequence
spiral2d.test();

%% Create a single-shot 2D spiral gradient-echo with specific spiral waveform
%-------------------------------------------------------------------------
% Compute trajectory
%-------------------------------------------------------------------------
Nitlv = 1;              % Number of interleves
r = 0;                  % rv/riv Indicates type of solution
res	= 1.8;              % Resolution (in mm)
fov	= [25 24 22];       % Vector of fov (in cm)
radius = [0,0.5,1];     % Vector of radius corresponding to the fov
Gmax = 4;               % Max gradient (default 3 G/CM = 30 mT/m)
Smax = 10;              % Max slew (default 10 G/cm/ms = 100 mT/m/ms)
T = 10e-3;              % Sampling rate (in ms) - 10 us on Siemens systems
ds = [];                % Step size for integration
interpType = 'cubic';   % Type of interpolation used to interpolate the fov accept: linear, cubic, spline

[k_rv,g_rv,s_rv,time_rv,Ck_rv] = vdSpiralDesign(Nitlv, r, res,fov,radius,Gmax,Smax,T,ds,interpType);

% Convert gradient to mT/m
g_rv = [0,0; g_rv(:,1:2); 0,0] * 10;
%-------------------------------------------------------------------------

paramsSpiral2d = SequenceParams('spiral2d',scannerType);
paramsSpiral2d.Ny = Nitlv;
paramsSpiral2d.Nx = 140; %from fov/res
paramsSpiral2d.fov = 250e-3; %from fov(1)
paramsSpiral2d.readoutTime = 70.1e-3; %from size(k_rv,1) * 1e-3
paramsSpiral2d.mode = 'singleShot';
paramsSpiral2d.seqSpecName = '';

% Create sequence
spiral2d = skope_spiral_2d(paramsSpiral2d,g_rv);

% Plot sequence information after sync 
timeRange = [5.4 5.49];
spiral2d.plot(timeRange);

% Test sequence
spiral2d.test();