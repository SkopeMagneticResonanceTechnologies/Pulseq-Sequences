classdef SequenceParams_Terra_X
    %SEQUENCEPARAMS_TERRA_X  All sequence parameters for the Siemens 7T Terra.X / Terra SC72CD.
    %   Called by SequenceParams for imaging sequences (gre2d, epi2d, etc.).
    %   Edit this file to adjust defaults specifically for this scanner.

    methods (Static)
        function obj = apply(obj, seqName, mode) %#ok<INUSD>
            switch lower(seqName)
                case 'gre2d'
                    obj.fov          = 220e-3;
                    obj.Nx           = 80;
                    obj.Ny           = 80;
                    obj.alpha        = 7;
                    obj.thickness    = 3e-3;
                    obj.nSlices      = 12;
                    obj.TE           = [6 12] * 1e-3;
                    obj.TR           = 25e-3;
                    obj.readoutTime  = 3.2e-3;
                    obj.maxGrad      = 130;
                    obj.maxSlew      = 190;
                    obj.nDummy       = 10;
                    obj.doMonitoringDuringRF = 0;
                    obj.distanceFactorPercentage = 150;
                case 'epi2d'
                    obj.fov          = 220e-3;
                    obj.TE           = 25e-3;
                    obj.TR           = 80e-3;
                    obj.alpha        = 90;
                    obj.Nx           = 120;
                    obj.Ny           = 120;
                    obj.thickness    = 3e-3;
                    obj.nSlices      = 12;
                    obj.nDummy       = 5;   % totalNofDummy=nDummy*nSlices (without FM trigger)
                    obj.nRep         = 1;
                    obj.addPhaseCorrLines = 1;
                    obj.readoutTime  = 500e-6;
                    obj.maxGrad      = 130;
                    obj.maxSlew      = 190;
                    obj.distanceFactorPercentage = 150;
                    obj.accFacPE = 2;
                    obj.multiBandFactor = 2;
                case 'se_epi2d_diff'
                    obj.TE           = 70e-3;
                    obj.TR           = 105e-3;
                    obj.alpha        = 90;
                    obj.fov          = 220e-3;
                    obj.Nx           = 146;
                    obj.Ny           = 146;
                    obj.thickness    = 1.5e-3;
                    obj.nSlices      = 30;
                    obj.nDummy       = 1;   % totalNofDummy=nDummy*nSlices*bEncoding (without FM trigger)
                    obj.nRep         = 1;
                    obj.accFacPE     = 2;
                    obj.addPhaseCorrLines = true;
                    obj.partFourierFactor = 7/8;
                    obj.doPlayFatSat = 0;
                    bvals = [2000];
                    dirs = [
                        1,0,0;                  % x 
                        0,1,0;                  % y
                        0,0,1;                  % z
                        [1,1,0]/sqrt(2);        % xy
                        [0,1,1]/sqrt(2);        % yz
                        [1,0,1]/sqrt(2)];       % xz
                    obj.bFactor = [0; repelem(bvals(:), size(dirs,1), 1)]; % b = 0 and the other values
                    obj.bDir = [0,0,0; repmat(dirs,numel(bvals),1)]; % b = 0 direction and the others
                    obj.nbValues     = size(obj.bDir,1);
                    obj.seqSpecName  = '';
                    obj.readoutTime  = 700e-6;
                    obj.maxGrad      = 130;
                    obj.maxSlew      = 170;
                    obj.maxDiffSlew  = 100; 
                    obj.maxDiffGrad  = 130;
                    obj.distanceFactorPercentage = 150;
                    obj.gradFreeTime = 300e-6;
                case 'se_spiral_2d_diff'
                    obj.TE           = 35e-3;
                    obj.TR           = 75e-3;
                    obj.alpha        = 90;
                    obj.thickness    = 1.5e-3;
                    obj.nSlices      = 30;
                    obj.mode         = 'MS';
                    obj.maxGrad      = 130;
                    obj.maxSlew      = 190;
                    obj.distanceFactorPercentage = 150;
                    obj.nDummy       = 1;
                    obj.accFac       = 3;
                    obj.doPlayFatSat = true;
                    bvals = [2000];
                    dirs = [
                        1,0,0;                  % x 
                        0,1,0;                  % y
                        0,0,1;                  % z
                        [1,1,0]/sqrt(2);        % xy
                        [0,1,1]/sqrt(2);        % yz
                        [1,0,1]/sqrt(2)];       % xz
                    obj.bFactor = [0; repelem(bvals(:), size(dirs,1), 1)]; % b = 0 and the other values
                    obj.bDir = [0,0,0; repmat(dirs,numel(bvals),1)]; % b = 0 direction and the others
                    obj.nbValues     = size(obj.bDir,1);
                    obj.maxGrad      = 130;
                    obj.maxSlew      = 160;
                    obj.maxDiffSlew  = 100;
                    obj.maxDiffGrad  = 130;
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
                    obj.maxGrad      = 130;
                    obj.maxSlew      = 190;
                    obj.nDummy       = 1;               %Total number of readout dummies = nDummy*Ny
                    obj.distanceFactorPercentage = 200;
                    obj.accFacPE     = 2;
                    obj.accFacZ      = 2;
                    obj.fullySampledZ = 0.2;
                    obj.doPlayFatSat = false;
                    obj.gradFreeTime = 200e-6;
                otherwise
                    error('SequenceParams_Terra_X: unknown sequence "%s".', seqName)
            end
        end
    end
end