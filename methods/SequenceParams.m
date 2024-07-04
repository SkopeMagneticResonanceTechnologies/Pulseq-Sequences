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
        scannerType = "Siemens Terra 7T SC72CD";
        nRep = 1;           % Number of repetitions
        nAve = 1;           % Number of averages  
        signFlip = -1;      % Bug fix for Pulseq error in version 1.4.0.
    end

    methods
        function obj = SequenceParams(seqName)
            switch lower(seqName)
                case 'gre2d'
                    obj.fov = 200e-3; 
                    obj.Nx = 128; 
                    obj.Ny = obj.Nx; 
                    obj.alpha = 7;   
                    obj.thickness = 3e-3; 
                    obj.nSlices = 5;
                    obj.TE = [6 12] * 1e-3;
                    obj.TR = 25e-3;       
                    obj.readoutTime = 3.2e-3;
                    obj.maxGrad = 28;
                    obj.maxSlew = 150;
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
                case 'epi2d'
                    obj.TE = 30e-3;
                    obj.TR = 150e-3;
                    obj.readoutTime = 0.8e-3;
                    obj.alpha = 90;
                    obj.fov = 256e-3;
                    obj.Nx = 64;
                    obj.Ny = 64;
                    obj.thickness = 4e-3;
                    obj.nSlices = 5;
                    obj.maxGrad = 32;
                    obj.maxSlew = 130;
                case 'gtf'
                    obj.TR = 1;  
                    obj.maxGrad = 40;
                    obj.maxSlew = 200;
                    obj.nAve = 5;
                case 'opc'
                    obj.TR = 200e-3;
                    obj.maxGrad = 40;
                    obj.maxSlew = 200;
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