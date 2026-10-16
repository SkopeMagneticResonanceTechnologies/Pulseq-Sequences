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
addpath('dependencies/pulseq_latest') % necessary Gradient Spectrum computation (pulseq 1.5.1+)

%% Define scanner type
% ../methods/GetMRSystemSpecs.m
% 'Siemens 3T Cima.X', 'Siemens 3T Connectom', 'Siemens 7T Terra SC72CD', 'Siemens 9.4T SC72CD'
scannerType = 'Siemens 3T Cima.X';

%% Acquire more slices for in vivo scans and activate fat saturation
invivo = false;

%% Create a 2D mono-polar dual-echo gradient-echo (GRE) sequence 
% Get default sequence parameters
paramsGre2d = SequenceParams('gre2d',scannerType);

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

%% Create a 2D echo-planar imaging (EPI) sequence with 2-fold in-plane acceleration and 2-fold multi-band excitation
% Get default sequence parameters
paramsEpi2d = SequenceParams('epi2d',scannerType);

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

%% Create a 2D spin-echo EPI sequence with diffusion encoding and 2-fold in-plane acceleration
% Get default sequence parameters
paramsSeEpi2dDiff = SequenceParams('se_epi2d_diff',scannerType); 

if invivo
    paramsSeEpi2dDiff.nSlices = 44;
    paramsSeEpi2dDiff.distanceFactorPercentage = 10;
    paramsSeEpi2dDiff.seqSpecName = 'invivo';
    paramsSeEpi2dDiff.doPlayFatSat = true;
end

% ---- b-encoding from external text file ----
% paramsSeEpi2dDiff.bDir = readmatrix('C:\My folders\Pulseq-Sequences\dependencies\bencoding\decompressed_dwepi_R2_MB2_PF_1.5mm_64b2000_20260724095017_6_bvec.txt');
% paramsSeEpi2dDiff.bFactor = readmatrix('C:\My folders\Pulseq-Sequences\dependencies\bencoding\decompressed_dwepi_R2_MB2_PF_1.5mm_64b2000_20260724095017_6_bval.txt');
% % obj.bDir must be N x 3
% if size(paramsSeEpi2dDiff.bDir,2) ~= 3
%     if size(paramsSeEpi2dDiff.bDir,1) == 3
%         paramsSeEpi2dDiff.bDir = paramsSeEpi2dDiff.bDir.';
%     else
%         error('bDir must have size N x 3 or 3 x N.');
%     end
% end
% % obj.bFactor must be N x 1
% paramsSeEpi2dDiff.bFactor = paramsSeEpi2dDiff.bFactor(:);
% paramsSeEpi2dDiff.nbValues = size(paramsSeEpi2dDiff.bDir,1);

% Generate the sequence
seepi2d = skope_se_epi_2d_diff(paramsSeEpi2dDiff);

% Plot sequence diagram for the first 4 to 20 s
timeRange = [4 20];
seepi2d.plot(timeRange);

% Test sequence
seepi2d.test();

%% Create a multi-shot 2D spin-echo spiral sequence with diffusion encoding and spiral arms undersampling with specific spiral waveform
minTimeGradientDir = fullfile(fileparts(pwd), 'minTimeGradient', 'Matlab');
if not(isfolder(minTimeGradientDir))
    error('Download https://people.eecs.berkeley.edu/~mlustig/software/tOptGrad_V0.2.tar.gz')
end
addpath(genpath(minTimeGradientDir))
%-------------------------------------------------------------------------
% Compute trajectory
%-------------------------------------------------------------------------
Nitlv = 2;                  % Number of interleves
r = 0;                      % rv/riv Indicates type of solution
res	= 1.2;                  % Resolution (in mm)
fov	= [22 22];              % Vector of fov (in cm)
radius = [0 1];             % Vector of radius corresponding to the fov
Gmax = 20;                  % Max gradient (default 3 G/CM = 30 mT/m)
Smax = 10;                  % Max slew (default 10 G/cm/ms = 100 mT/m/ms)
T = 10e-3;                  % Sampling rate (in ms) - 10 us on Siemens systems
ds = [];                    % Step size for integration
interpType = 'linear';       % Type of interpolation used to interpolate the fov accept: linear, cubic, spline % cubic for variable density

[k_rv,g_rv,s_rv,time_rv,Ck_rv] = vdSpiralDesign(Nitlv, r, res,fov,radius,Gmax,Smax,T,ds,interpType);

% Convert gradient to mT/m
g_rv = [0,0; g_rv(:,1:2); 0,0] * 10;
figure, plot(g_rv), xlabel('datapoints'), ylabel('gradients [mt/m]'), title('Gradient waveforms')
figure, plot(k_rv(:,1), k_rv(:,2)), xlabel('kx'), ylabel('ky'), title('Spiral trajectory')
% size(g_rv,1)
%-------------------------------------------------------------------------

paramsSpiral2d = SequenceParams('se_spiral_2d_diff',scannerType);
paramsSpiral2d.Ny = Nitlv;
paramsSpiral2d.Nx = ceil(fov(1)./res);
paramsSpiral2d.fov = fov(1)*1e-2; 
paramsSpiral2d.readoutTime = time_rv; 
paramsSpiral2d.mode = 'MS';
% paramsSpiral2d.nSlices = 2;
paramsSpiral2d.accFac = 2;
paramsSpiral2d.TE = 31e-3; 
paramsSpiral2d.TR = 105e-3; 

% ---- b-encoding from external file ----
% paramsSpiral2d.bDir = readmatrix('C:\My folders\Pulseq-Sequences\dependencies\bencoding\decompressed_dwepi_R2_MB2_PF_1.5mm_64b2000_20260724095017_6_bvec.txt');
% paramsSpiral2d.bFactor = readmatrix('C:\My folders\Pulseq-Sequences\dependencies\bencoding\decompressed_dwepi_R2_MB2_PF_1.5mm_64b2000_20260724095017_6_bval.txt');
% % obj.bDir must be N x 3
% if size(paramsSpiral2d.bDir,2) ~= 3
%     if size(paramsSpiral2d.bDir,1) == 3
%         paramsSpiral2d.bDir = paramsSpiral2d.bDir.';
%     else
%         error('bDir must have size N x 3 or 3 x N.');
%     end
% end
% % obj.bFactor must be N x 1
% paramsSpiral2d.bFactor = paramsSpiral2d.bFactor(:);
% paramsSpiral2d.nbValues = size(paramsSpiral2d.bDir,1);

% Create sequence
spiral2d = skope_se_spiral_2d_diff(paramsSpiral2d,g_rv);

% Plot sequence information after sync 
% timeRange = [5.4 5.42];
% spiral2d.plot(timeRange);

% Test sequence
spiral2d.test();

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
% paramsLec.triggerOutput = 'osc0';
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