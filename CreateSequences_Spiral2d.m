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
    error("SAFE PNS Prediction submodule not found. Please run 'git submodule init' and 'git submodule update'.")
end

%% Add Pulseq, sequences and methods
addpath('pulseq/matlab')
addpath('methods')
addpath('sequences')
addpath('safe_pns_prediction')
addpath('dependencies/pulseq_latest') % necessary for CNS and Gradient Spectrum computation

%% Define scanner type
% 'Siemens 3T Cima.X', 'Siemens 7T Terra SC72CD', 'Siemens 9.4T SC72CD'
scannerType = 'Siemens 3T Cima.X';

% %% Create a multi-shot 2D spiral gradient-echo sequence with pre-generated spiral waveform
% % Get default sequence parameters
% paramsSpiral2d = SequenceParams('spiral2d',scannerType);
% load('./waveforms/spiralGrad_FOV192_RES1mm_minRise6_maxAmp40_nitlv16.mat'); % [Hz/m]
% obj.sys.gamma = 42576000;
% spiralWaveform = spiralWaveform / obj.sys.gamma * 1000;
% paramsSpiral2d.seqSpecName = 'hardcoded';
% 
% spiral2d = skope_spiral_2d(paramsSpiral2d,spiralWaveform);
% 
% % Plot sequence information after sync 
% timeRange = [5.4 5.42];
% spiral2d.plot(timeRange);
% 
% % Test sequence
% spiral2d.test();

%% Create a multi-shot 2D spiral spin-echo sequence for diffusion with specific spiral waveform
minTimeGradientDir = fullfile(fileparts(pwd), 'minTimeGradient', 'Matlab');
if not(isfolder(minTimeGradientDir))
    error('Download https://people.eecs.berkeley.edu/~mlustig/software/tOptGrad_V0.2.tar.gz')
end
addpath(genpath(minTimeGradientDir))


%%
%-------------------------------------------------------------------------
% Compute trajectory
%-------------------------------------------------------------------------
Nitlv = 16;             % Number of interleves
r = 0;                  % rv/riv Indicates type of solution
res	= 1.2;              % Resolution (in mm)
fov	= [22,22];          % Vector of fov (in cm)
radius = [0,1];         % Vector of radius corresponding to the fov
Gmax = 10;               % Max gradient (default 3 G/CM = 30 mT/m)
Smax = 10;               % Max slew (default 10 G/cm/ms = 100 mT/m/ms)
T = 10e-3;              % Sampling rate (in ms) - 10 us on Siemens systems
ds = [];                % Step size for integration
interpType = 'linear';  % Type of interpolation used to interpolate the fov accept: linear, cubic, spline

[k_rv,g_rv,s_rv,time_rv,Ck_rv] = vdSpiralDesign(Nitlv, r, res,fov,radius,Gmax,Smax,T,ds,interpType);

% Convert gradient to mT/m
g_rv = [0,0; g_rv(:,1:2); 0,0] * 10;
figure, plot(g_rv), xlabel('datapoints'), ylabel('gradients [mt/m]'), title('Gradient waveforms')
figure, plot(k_rv(:,1), k_rv(:,2)), xlabel('kx'), ylabel('ky'), title('Spiral trajectory')
% size(g_rv,1)
%-------------------------------------------------------------------------

%% 
% see the difference between spiral arms = 8 and 16 (for undersampling purposes)
% Nitlv = 8;
% [k8,g8,s_rv,time_rv,Ck_rv] = vdSpiralDesign(Nitlv, r, res,fov,radius,Gmax,Smax,T,ds,interpType);
% 
% Nitlv = 16;
% [k16,g16,s_rv,time_rv,Ck_rv] = vdSpiralDesign(Nitlv, r, res,fov,radius,Gmax,Smax,T,ds,interpType);
% plot(k8(:,1),k8(:,2))
% hold on
% plot(k16(:,1),k16(:,2))

%%
paramsSpiral2d = SequenceParams('se_spiral_2d_diff',scannerType);
paramsSpiral2d.Ny = Nitlv;
paramsSpiral2d.Nx = ceil(fov(1)./res);
paramsSpiral2d.fov = fov(1)*1e-2; 
paramsSpiral2d.readoutTime = time_rv; 
paramsSpiral2d.mode = 'MS';
% paramsSpiral2d.nSlices = 2;
paramsSpiral2d.nDummy = 1;
paramsSpiral2d.accFac = 1;
paramsSpiral2d.doPlayFatSat = true;
paramsSpiral2d.TE = 34e-3;
paramsSpiral2d.TR = 120e-3;

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

% Sequence name
moreName = '';
paramsSpiral2d.seqSpecName = sprintf( ...
    '%ddir_b%d_te%d_tr%d_%.1fmm%s', ...
    paramsSpiral2d.nbValues - 1, ...
    max(paramsSpiral2d.bFactor), ...
    round(paramsSpiral2d.TE*1e3), ...
    round(paramsSpiral2d.TR*1e3), ...
    res, ...
    moreName);

fprintf('seqSpecName = %s\n', paramsSpiral2d.seqSpecName);

% Create sequence
spiral2d = skope_se_spiral_2d_diff(paramsSpiral2d,g_rv);

% Plot sequence information after sync 
% timeRange = [5.4 5.42];
% spiral2d.plot(timeRange);

% Test sequence
% spiral2d.test();

%% Create a single-shot 2D spiral spin-echo sequence for diffusion with specific spiral waveform
%-------------------------------------------------------------------------
% Compute trajectory
%-------------------------------------------------------------------------
Nitlv = 1;              % Number of interleves
r = 0;                  % rv/riv Indicates type of solution
res	= 1;                % Resolution (in mm)
fov	= [22 22 22];       % Vector of fov (in cm)
radius = [0,0.5,1];     % Vector of radius corresponding to the fov
Gmax = 15;              % Max gradient (default 3 G/CM = 30 mT/m)
Smax = 9.5;              % Max slew (default 10 G/cm/ms = 100 mT/m/ms)
T = 10e-3;              % Sampling rate (in ms) - 10 us on Siemens systems
ds = [];                % Step size for integration
interpType = 'cubic';   % Type of interpolation used to interpolate the fov accept: linear, cubic, spline

[k_rv,g_rv,s_rv,time_rv,Ck_rv] = vdSpiralDesign(Nitlv, r, res,fov,radius,Gmax,Smax,T,ds,interpType);

% Convert gradient to mT/m
g_rv = [0,0; g_rv(:,1:2); 0,0] * 10;
figure, plot(g_rv), xlabel('datapoints'), ylabel('gradients [mt/m]'), title('Gradient waveforms')
figure, plot(k_rv(:,1), k_rv(:,2)), xlabel('kx'), ylabel('ky'), title('Spiral trajectory')


%%
paramsSpiral2d = SequenceParams('se_spiral_2d_diff',scannerType);
paramsSpiral2d.Ny = Nitlv;
paramsSpiral2d.Nx = ceil(fov(1)./res);
paramsSpiral2d.fov = fov(1)*1e-2; 
paramsSpiral2d.readoutTime = time_rv; %from size(k_rv,1) * 1e-3
paramsSpiral2d.mode = 'SS';
paramsSpiral2d.nDummy = 1;
% paramsSpiral2d.nSlices = 2;
paramsSpiral2d.doPlayFatSat = true;
paramsSpiral2d.TE = 35e-3;
paramsSpiral2d.TR = 210e-3;

% Sequence name
moreName = '';
paramsSpiral2d.seqSpecName = sprintf( ...
    '%ddir_b%d_te%d_tr%d_%.1fmm%s', ...
    paramsSpiral2d.nbValues - 1, ...
    max(paramsSpiral2d.bFactor), ...
    round(paramsSpiral2d.TE*1e3), ...
    round(paramsSpiral2d.TR*1e3), ...
    res, ...
    moreName);

fprintf('seqSpecName = %s\n', paramsSpiral2d.seqSpecName);

% Create sequence
spiral2d = skope_se_spiral_2d_diff(paramsSpiral2d,g_rv);

% Plot sequence information after sync 
% timeRange = [5.4 5.49];
% spiral2d.plot(timeRange);

% Test sequence
% spiral2d.test();