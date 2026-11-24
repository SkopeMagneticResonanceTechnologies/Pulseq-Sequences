classdef SequenceParams_Cima_X
    %SEQUENCEPARAMS_CIMA_X  All sequence parameters for the Siemens 3T Cima.X.
    %   Called by SequenceParams for imaging sequences (gre2d, epi2d, etc.).
    %   Edit this file to adjust defaults specifically for this scanner.

    methods (Static)
        function obj = apply(obj, seqName, mode) %#ok<INUSD>
            switch lower(seqName)
                case 'gre2d'
                    obj.fov          = 240e-3;
                    obj.Nx           = 120;
                    obj.Ny           = 120;
                    obj.alpha        = 7;
                    obj.thickness    = 2e-3;
                    obj.nSlices      = 36;
                    obj.TE           = [7 14] * 1e-3;
                    obj.TR           = 25e-3;
                    obj.readoutTime  = 4.8e-3;
                    obj.maxGrad      = 28;
                    obj.maxSlew      = 150;
                    obj.nDummy       = 5;
                    obj.doMonitoringDuringRF = false;
                    obj.distanceFactorPercentage = 150;
                    obj.gradFreeTime = 200e-6;
                case 'epi2d'
                    obj.TE           = 25e-3;
                    obj.TR           = 60e-3;
                    obj.alpha        = 90;
                    obj.fov          = 220e-3;
                    obj.Nx           = 120; %130
                    obj.Ny           = 120; %130
                    obj.thickness    = 3e-3;
                    obj.nSlices      = 12;
                    obj.nDummy       = 5;   % totalNofDummy=nDummy*nSlices (without FM trigger)
                    obj.nRep         = 1;
                    obj.addPhaseCorrLines = 1;
                    obj.accFacPE     = 2; %3
                    obj.multiBandFactor = 2;
                    obj.readoutTime  = 0.500e-3; %660
                    obj.maxGrad      = 190;
                    obj.maxSlew      = 190;
                    obj.distanceFactorPercentage = 150;
                    obj.doPlayFatSat = true;
                    obj.gradFreeTime = 200e-6;
                case 'se_epi2d_diff'
                    obj.TE           = 84e-3;
                    obj.TR           = 122e-3;
                    obj.alpha        = 90;
                    obj.fov          = 220e-3;
                    obj.Nx           = 146;
                    obj.Ny           = 146;
                    obj.thickness    = 1.5e-3;
                    obj.nSlices      = 30;
                    obj.nDummy       = 1;   % no trigger played out. Total number of readout dummies = nDummy*nSlices.
                    obj.nRep         = 1;
                    obj.accFacPE     = 2;
                    obj.partFourierFactor = 7/8;
                    obj.addPhaseCorrLines = true;
                    obj.doPlayFatSat = true;
                    %bvals = [2000 4000 6000 8000 10000];
                    bvals = [10000];
                    dirs = [
                        1,0,0;                  % x 
                        0,1,0;                  % y
                        0,0,1;                  % z
                        [1,1,0]/sqrt(2);        % xy
                        [0,1,1]/sqrt(2);        % yz
                        [1,0,1]/sqrt(2)];       % xz
                    %obj.bFactor = [0; repelem(bvals(:), size(dirs,1), 1)]; % b = 0 and the other values
                    obj.bFactor = [repelem(bvals(:), size(dirs,1), 1)];
                    %obj.bDir = [0,0,0; repmat(dirs,numel(bvals),1)]; % b = 0 direction and the others
                    obj.bDir = [repmat(dirs,numel(bvals),1)];
                    obj.nbValues     = size(obj.bDir,1);
                    obj.seqSpecName  = '';
                    obj.readoutTime  = 0.600e-3;
                    obj.maxGrad      = 190;
                    obj.maxSlew      = 190;
                    obj.maxDiffSlew  = 63; 
                    obj.maxDiffGrad  = 172;
                    obj.distanceFactorPercentage = 150;
                    obj.gradFreeTime = 300e-6;
                case 'se_spiral_2d_diff'
                    obj.TE           = 53e-3;
                    obj.TR           = 95e-3;
                    obj.alpha        = 90;
                    obj.thickness    = 1.5e-3;
                    obj.nSlices      = 30;
                    obj.mode         = 'MS';
                    obj.maxGrad      = 190;
                    obj.maxSlew      = 190;
                    obj.distanceFactorPercentage = 150;
                    obj.nDummy       = 1;
                    obj.accFac       = 3;
                    obj.doPlayFatSat = true;
                    %bvals = [2000 4000 6000 8000 10000];
                    bvals = [10000 10000 10000 10000 10000 10000];
                    dirs = [
                        1,0,0;                  % x 
                        0,1,0;                  % y
                        0,0,1;                  % z
                        [1,1,0]/sqrt(2);        % xy
                        [0,1,1]/sqrt(2);        % yz
                        [1,0,1]/sqrt(2)];       % xz
                    %obj.bFactor = [0; repelem(bvals(:), size(dirs,1), 1)]; % b = 0 and the other values
                    obj.bFactor = [repelem(bvals(:), size(dirs,1), 1)];
                    %obj.bDir = [0,0,0; repmat(dirs,numel(bvals),1)]; % b = 0 direction and the others
                    obj.bDir = [repmat(dirs,numel(bvals),1)];
                    obj.nbValues     = size(obj.bDir,1);
                    obj.maxDiffSlew  = 63;
                    obj.maxDiffGrad  = 172;
                    obj.gradFreeTime = 300e-6;
                    obj.seqSpecName = '';
                case 'gre3d'
                    obj.fov          = [0.256 0.224 0.176];
                    obj.Nx           = 256;
                    obj.Ny           = 224;
                    obj.Nz           = 160;
                    obj.alpha        = 15;
                    obj.TE           = [6.7 11.9 17.2 22.4 27.6 32.9]*1e-3;
                    obj.TR           = 43e-3;
                    obj.readoutTime  = 12.5e-6 * obj.Nx;
                    obj.maxGrad      = 190;
                    obj.maxSlew      = 190;
                    obj.nDummy       = 1;               %Total number of readout dummies = nDummy*Ny
                    obj.distanceFactorPercentage = 200;
                    obj.accFacPE     = 2;
                    obj.accFacZ      = 2;
                    obj.fullySampledZ = 0.2;
                    obj.doPlayFatSat = false;
                    obj.gradFreeTime = 200e-6;
                otherwise
                    error('SequenceParams_Cima_X: unknown sequence "%s".', seqName)
            end
        end
    end
end