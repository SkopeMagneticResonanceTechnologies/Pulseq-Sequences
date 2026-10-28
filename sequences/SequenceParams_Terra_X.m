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
                    obj.TE           = 65e-3;
                    obj.TR           = 140e-3;
                    obj.alpha        = 90;
                    obj.fov          = 220e-3;
                    obj.Nx           = 120;
                    obj.Ny           = 120;
                    obj.thickness    = 3e-3;
                    obj.nSlices      = 12;
                    obj.nDummy       = 5;   % totalNofDummy=nDummy*nSlices*bEncoding (without FM trigger)
                    obj.nRep         = 1;
                    obj.accFacPE     = 3;
                    obj.addPhaseCorrLines = true;
                    obj.partFourierFactor = 7/8;
                    obj.doPlayFatSat = 0;
                    obj.bFactor      = [0, 2000, 2000, 2000, 2000./sqrt(2), 2000./sqrt(2), 2000./sqrt(2)];
                    obj.bDir         = [0,0,0;... % b0
                                        1,0,0;... % x
                                        0,1,0;... % y
                                        0,0,1;... % z
                                        1,1,0;... % xy
                                        0,1,1;... % yz
                                        1,0,1;];  % xz
                    obj.nbValues     = size(obj.bDir,1);
                    obj.seqSpecName  = '';
                    obj.readoutTime  = 500e-6;
                    obj.maxGrad      = 130;
                    obj.maxSlew      = 190;
                    obj.maxDiffSlew  = 15;
                    obj.distanceFactorPercentage = 150;
                case 'gre3d'
                    obj.fov          = [0.56 0.56 0.56]*1e-2*40053000/42577481;
                    obj.Nx           = 56;
                    obj.Ny           = 56;
                    obj.Nz           = 56;
                    obj.alpha        = 1;
                    obj.TE           = [12.3 28.16]*1e-3 + 1e-3;
                    obj.TR           = 100e-3;
                    obj.readoutTime  = 7.84e-3;
                    obj.maxGrad      = 130;
                    obj.maxSlew      = 190;
                    obj.nDummy       = 50;
                case 'spiral2d'
                    obj.fov          = 192e-3;
                    obj.Nx           = 192;
                    obj.Ny           = 16;
                    obj.alpha        = 15;
                    obj.thickness    = 3e-3;
                    obj.nSlices      = 15;
                    obj.TE           = 2.5e-3;
                    obj.TR           = 140e-3;
                    obj.readoutTime  = 8e-3;
                    obj.mode         = 'multiShot';
                    obj.maxGrad      = 130;
                    obj.maxSlew      = 190;
                    obj.distanceFactorPercentage = 200;
                otherwise
                    error('SequenceParams_Terra_X: unknown sequence "%s".', seqName)
            end
        end
    end
end