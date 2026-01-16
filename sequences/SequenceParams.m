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
        seqSpecName % String name appended to .seq file
        doMonitoringDuringRF %boolean to enable monitoring during RFs (functionality active for GRE only)

        %diffusion properties
        bFactor 
        bDir
        nbValues
        
        % Defaults
        scannerType = 'Siemens 9.4T SC72CD';
        nRep = 1;               % Number of repetitions
        nAve = 1;               % Number of averages  
        mode = 'default';       % Allow to switch between different versions
        doPlayFatSat = false;   % Play out fat-saturation pulse (for EPI)
        nDummy = 0;             % Number of dummy pulses to reach steady state
        accFacPE = 1;           % Acceleration factor [Phase] (only used for EPI at the moment)
        
        sliceOrientation = SliceOrientation.TRA;  % Slice orientation        
        phaseEncDir = PhaseEncodingDirection.AP;  % Phase encoding direction

        triggerOutput = 'ext1'; % Default optical output (other ouputs: 'osc0', 'osc1')
    end

    methods
        function obj = SequenceParams(seqName,scannerType,mode)

            if not(exist('scannerType','var'))
                scannerType = 'Siemens 7T Terra SC72CD';
                warning('Using default scanner type "Siemens 7T Terra SC72CD"')
            end
            obj.scannerType = scannerType;
            
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
                    obj.nSlices = 15;
                    obj.TE = [6 12] * 1e-3;
                    obj.TR = 25e-3;       
                    obj.readoutTime = 3.2e-3;
                    obj.maxGrad = 28;
                    obj.maxSlew = 150;
                    obj.nDummy = 10;
                    obj.doMonitoringDuringRF = 0;
                case 'epi2d'
                    obj.TE = 36e-3;
                    obj.TR = 200e-3;
                    switch scannerType
                        case 'Siemens 3T Cima.X'
                            obj.readoutTime = 0.500e-3;
                        otherwise
                            obj.readoutTime = 0.680e-3;
                    end
                    obj.alpha = 90;
                    obj.fov = 200e-3;
                    obj.Nx = 80;
                    obj.Ny = 80;
                    obj.thickness = 3e-3;
                    obj.nSlices = 15;
                    obj.maxGrad = 32;
                    obj.maxSlew = 180;
                    obj.nDummy = 4;   % totalNofDummy=nDummy*nSlices (without FM trigger)
                    obj.nRep = 10;
                case 'se_epi2d_diff'
                    obj.TE = 110e-3;
                    obj.TR = 200e-3;
                    obj.readoutTime = 0.680e-3;
                    obj.alpha = 90;
                    obj.fov = 200e-3;
                    obj.Nx = 80;
                    obj.Ny = 80;
                    obj.thickness = 3e-3;
                    obj.nSlices = 1;
                    obj.maxGrad = 32;
                    obj.maxSlew = 130;
                    obj.nDummy = 5;         % totalNofDummy=nDummy*nSlices*bEncoding (without FM trigger)
                    obj.accFacPE = 1;       % Acceleration factor [Phase] (only used for EPI at the moment)
                    obj.doPlayFatSat = 0;
                    obj.bFactor=[0, 1000, 1000, 1000]; %bencoding
                    obj.bDir = [0, 1, 2, 3]; %axis
                    obj.nbValues = length(obj.bDir);     
                    obj.seqSpecName = '';
                case 'gre3d'
                    obj.fov = [0.56 0.56 0.56]*1e-2*40053000/42577481; 
                    obj.Nx = 56; 
                    obj.Ny = obj.Nx; 
                    obj.Nz = obj.Nx; 
                    obj.alpha = 1;     
                    obj.TE = [12.3 28.16] * 1e-3 + 1e-3; % one millisecond for phase estimation
                    obj.TR = 100e-3;   
                    obj.readoutTime = 7.84e-3;  
                    obj.maxGrad = 35;
                    obj.maxSlew = 150;
                    obj.nDummy = 50;
                case 'spiral2d' %to be changed.
                    % spiral-trajectory not adaptive to input params (hard coded)
                    % don't change!
                    obj.fov = 192e-3; %to be changed.
                    obj.Nx = 192; %to be changed.
                    obj.Ny = 16; 
                    obj.alpha = 15;   
                    obj.thickness = 3e-3; 
                    obj.nSlices = 15;
                    obj.TE = 2.5 * 1e-3;
                    obj.TR = 140e-3;       
                    obj.readoutTime = 8e-3; %to be changed.
                    obj.mode = 'multiShot';
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