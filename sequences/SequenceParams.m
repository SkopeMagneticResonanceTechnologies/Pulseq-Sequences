classdef SequenceParams
    %SEQUENCEPARAMS Default sequence parameters.
    %   SequenceParams(seqName, scannerType, mode) handles scanner-agnostic
    %   sequences (gtf, opc, sweep, lec) directly, and delegates all
    %   imaging sequences to the appropriate scanner-specific parameter
    %   class (SequenceParams_Cima_X, SequenceParams_Terra_X, ...).

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
        addPhaseCorrLines
        distanceFactorPercentage = 0;

        %diffusion properties
        bFactor 
        bDir
        nbValues
        maxDiffSlew
        
        % Defaults
        scannerType = 'Siemens 9.4T SC72CD';
        nRep = 1;               % Number of repetitions
        nAve = 1;               % Number of averages  
        mode = 'default';       % Allow to switch between different versions
        doPlayFatSat = false;   % Play out fat-saturation pulse (for EPI)
        nDummy = 0;             % Number of dummy pulses to reach steady state
        accFacPE = 1;           % Acceleration factor [Phase] (only used for EPI at the moment)
        multiBandFactor = 1
        partFourierFactor = 1
        
        sliceOrientation = SliceOrientation.TRA;  % Slice orientation        
        phaseEncDir = PhaseEncodingDirection.AP;  % Phase encoding direction

        triggerOutput = 'ext1'; % Default optical output (other ouputs: 'osc0', 'osc1')
    end

    methods
        function obj = SequenceParams(seqName, scannerType, mode)

            if not(exist('scannerType','var'))
                scannerType = 'Siemens 7T Terra SC72CD';
                warning('Using default scanner type "Siemens 7T Terra SC72CD"')
            end
            obj.scannerType = scannerType;
            
            if not(exist('mode','var'))
                mode = 'default';
            end

            % ---- Scanner-agnostic sequences (same on all hardware) ------
            switch lower(seqName)
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
                    return
                case 'opc'
                    obj.TR = 200e-3;
                    obj.maxGrad = 40;
                    obj.maxSlew = 200;
                    return
                case 'sweep'
                    obj.TR = 1;
                    obj.maxGrad = 40;
                    obj.maxSlew = 200;
                    obj.nAve = 100;
                    return
                case 'lec'
                    obj.TR = 200e-3;
                    obj.maxGrad = 40;
                    obj.maxSlew = 200;
                    return
            end

            % ---- Imaging sequences: delegate fully to scanner class -----
            switch scannerType
                case 'Siemens 3T Cima.X'
                    obj = SequenceParams_Cima_X.apply(obj, seqName, mode);
                case {'Siemens 7T Terra.X', 'Siemens 7T Terra SC72CD'}
                    obj = SequenceParams_Terra_X.apply(obj, seqName, mode);
                case {'United 5T uMR Jupiter'}
                    obj = SequenceParams_uM_Jupiter_5T.apply(obj, seqName, mode);
                otherwise
                    error('No parameters defined for scanner "%s" and sequence "%s".', scannerType, seqName);
            end
        end
    end
end