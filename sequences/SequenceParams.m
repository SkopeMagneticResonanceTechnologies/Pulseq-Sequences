classdef SequenceParams
    %SEQUENCEPARAMS Default sequence parameters

    properties        
        fov         % Field of view [Unit: m]
        Nx          % Number of readout samples
        Ny          % Number of phase encoding steps     
        Nz          % Number of 3D phase encoding steps 
        alpha       % Flip angle [Unit: deg]
        thickness   % Slice thickness [Unit: m]
        nSlices     % Number of slices
        TE          % Echo times [Unit: s]
        TR          % Excitation repetition time [Unit: s]                 
        readoutTime % ADC duration [Unit: s]
        maxGrad     % Used gradient amplitude by sequence
        maxSlew     % Used slew rate by sequence
        
        % Defaults
        scannerType = 'Siemens Terra 7T SC72CD';
        nRep = 1;               % Number of repetitions
        nAve = 1;               % Number of averages  
        mode = 'default';       % Allow to switch between different versions
        doPlayFatSat = true;    % Play out fat-saturation pulse (for EPI)
        nDummy = 0;             % Number of dummy pulses to reach steady state
        accFacPE = 1;           % Acceleration factor [Phase] (only used for EPI at the moment)
        
        sliceOrientation = SliceOrientation.TRA;  % Slice orientation        
        phaseEncDir = PhaseEncodingDirection.AP;  % Phase encoding direction
    end

    methods
        function obj = SequenceParams(seqName,mode)

            if not(exist('mode','var'))
                mode = 'default';
            end

            switch lower(seqName)
                case 'gre2d'
                    obj.fov = 200e-3; 
                    obj.Nx = 128; 
                    obj.Ny = obj.Nx; 
                    obj.alpha = 7;   
                    obj.thickness = 3e-3; 
                    obj.nSlices = 50;
                    obj.TE = [6 12] * 1e-3;
                    obj.TR = 25e-3;       
                    obj.readoutTime = 3.2e-3;
                    obj.maxGrad = 28;
                    obj.maxSlew = 150;
                    obj.nDummy = 10;
                case 'epi2d'
                    obj.TE = 33e-3;
                    obj.TR = 200e-3;
                    obj.readoutTime = 0.680e-3;
                    obj.alpha = 90;
                    obj.fov = 200e-3;
                    obj.Nx = 80;
                    obj.Ny = 80;
                    obj.thickness = 3e-3;
                    obj.nSlices = 15;
                    obj.maxGrad = 32;
                    obj.maxSlew = 130;
                    obj.nDummy = 5;   
                    obj.nRep = 10;
                case 'epi2d_diff'
                    obj.TE = 68e-3;
                    obj.TR = 10; % Volume TR
                    obj.readoutTime = 5.2e-4;
                    obj.alpha = 90;
                    obj.fov = 220e-3;
                    obj.Nx = 140;
                    obj.Ny = 140;
                    obj.thickness = 1.5e-3;
                    obj.nSlices = 100;
                    obj.maxGrad = 38;
                    obj.maxSlew = 180;
                    obj.nDummy = 1;   
                    obj.nRep = 1;
                    obj.accFacPE = 3;
                case 'gre3d'
                    obj.fov = [0.22 0.22 0.22]; 
                    obj.Nx = 100; % Divisible by 4
                    obj.Ny = obj.Nx; 
                    obj.Nz = obj.Nx; 
                    obj.alpha = 1;     
                    obj.TE = [5.2000    9.7000   14.2000   18.7000   23.2000   27.7000   32.2000   36.7000] * 1e-3;
                    obj.TR = 50e-3;   
                    obj.readoutTime = 14e-6*obj.Nx;  % 14us dwell time
                    obj.maxGrad = 35;
                    obj.maxSlew = 150;
                    obj.nDummy = 50;
                case 'spiral2d'
                    % spiral-trajectory not adaptive to input params (hard coded)
                    % don't change!
                    obj.fov = 192e-3; 
                    obj.Nx = 192; 
                    obj.Ny = 16; 
                    obj.alpha = 15;   
                    obj.thickness = 3e-3; 
                    obj.nSlices = 15;
                    obj.TE = 2.5 * 1e-3;
                    obj.TR = 140e-3;       
                    obj.readoutTime = 8e-3;
                    obj.maxGrad = 40;
                    obj.maxSlew = 150;
                case 'gtf'
                    if strcmpi(mode,'default')
                        obj.TR = 1;  
                        obj.maxGrad = 40;
                        obj.maxSlew = 200;
                        obj.nAve = 5;
                        obj.mode = mode;
                    elseif strcmpi(mode,'linearityCheck')
                        obj.TR = 1;  
                        obj.maxGrad = 40;
                        obj.maxSlew = 200;
                        obj.nAve = 1;
                        obj.mode = mode;
                    else
                        error('Unknown sequence mode.')
                    end
                case 'opc'
                    obj.TR = 200e-3;
                    obj.maxGrad = 40;
                    obj.maxSlew = 200;
                case 'sweep'
                    obj.TR = 1;
                    obj.maxGrad = 40;
                    obj.maxSlew = 200;
                    obj.nAve = 100;
                case 'lec'
                    obj.TR = 200e-3;
                    obj.maxGrad = 40;
                    obj.maxSlew = 200;
                otherwise
                    error(['Unknown sequence: ', seqName])
            end
        end
    end
end