classdef SequenceParams_Cima_X
    %SEQUENCEPARAMS_CIMA_X  All sequence parameters for the Siemens 3T Cima.X.
    %   Called by SequenceParams for imaging sequences (gre2d, epi2d, etc.).
    %   Edit this file to adjust defaults specifically for this scanner.

    methods (Static)
        function obj = apply(obj, seqName, mode) %#ok<INUSD>
            switch lower(seqName)
                case 'gre2d'
                    obj.fov          = 220e-3;
                    obj.Nx           = 128;
                    obj.Ny           = 128;
                    obj.alpha        = 7;
                    obj.thickness    = 3e-3;
                    obj.nSlices      = 12;
                    obj.TE           = [6 12] * 1e-3;
                    obj.TR           = 25e-3;
                    obj.readoutTime  = 3.2e-3;
                    obj.maxGrad      = 28;
                    obj.maxSlew      = 150;
                    obj.nDummy       = 5;
                    obj.doMonitoringDuringRF = 0;
                    obj.distanceFactorPercentage = 150;
                case 'epi2d'
                    obj.TE           = 33e-3;
                    obj.TR           = 130e-3;
                    obj.alpha        = 90;
                    obj.fov          = 220e-3;
                    obj.Nx           = 96;
                    obj.Ny           = 96;
                    obj.thickness    = 3e-3;
                    obj.nSlices      = 12;
                    obj.nDummy       = 5;   % totalNofDummy=nDummy*nSlices (without FM trigger)
                    obj.nRep         = 1;
                    obj.addPhaseCorrLines = 1;
                    obj.readoutTime  = 0.500e-3;
                    obj.maxGrad      = 190;
                    obj.maxSlew      = 190;
                    obj.distanceFactorPercentage = 150;
                case 'se_epi2d_diff'
                    obj.TE           = 64e-3;
                    obj.TR           = 130e-3;
                    obj.alpha        = 90;
                    obj.fov          = 220e-3;
                    obj.Nx           = 80;
                    obj.Ny           = 80;
                    obj.thickness    = 3e-3;
                    obj.nSlices      = 12;
                    obj.nDummy       = 5;   % totalNofDummy=nDummy*nSlices*bEncoding (without FM trigger)
                    obj.nRep         = 1;
                    obj.accFacPE     = 1;
                    obj.addPhaseCorrLines = 1;
                    obj.doPlayFatSat = 0;
                    obj.bFactor      = [0, 2000, 2000, 2000, 2000, 2000, 2000];
                    obj.bDir         = [0,0,0;... % b0
                                        1,0,0;... % x
                                        0,1,0;... % y
                                        0,0,1;... % z
                                        [1,1,0]./sqrt(2);... % xy
                                        [0,1,1]./sqrt(2);... % yz
                                        [1,0,1]./sqrt(2);];  % xz
                    obj.nbValues     = size(obj.bDir,1);
                    obj.seqSpecName  = '';
                    obj.readoutTime  = 0.600e-3;
                    obj.maxGrad      = 190;
                    obj.maxSlew      = 190;
                    obj.maxDiffSlew  = 120;
                    obj.distanceFactorPercentage = 150;
                case 'se_epi2d_diff_mb'
                    obj.TE           = 64e-3;
                    obj.TR           = 130e-3;
                    obj.alpha        = 90;
                    obj.fov          = 220e-3;
                    obj.Nx           = 80;
                    obj.Ny           = 80;
                    obj.thickness    = 3e-3;
                    obj.nSlices      = 12;
                    obj.nDummy       = 5;   % totalNofDummy=nDummy*nSlices*bEncoding (without FM trigger)
                    obj.accFacPE     = 1;
                    obj.addPhaseCorrLines = 1;
                    obj.doPlayFatSat = 0;
                    obj.bFactor      = [0, 2000, 2000, 2000, 2000, 2000, 2000];
                    obj.bDir         = [0,0,0;... % b0
                                        1,0,0;... % x
                                        0,1,0;... % y
                                        0,0,1;... % z
                                        [1,1,0]./sqrt(2);... % xy
                                        [0,1,1]./sqrt(2);... % yz
                                        [1,0,1]./sqrt(2);];  % xz
                    obj.nbValues     = size(obj.bDir,1);
                    obj.seqSpecName  = '';
                    obj.readoutTime  = 0.600e-3;
                    obj.maxGrad      = 190;
                    obj.maxSlew      = 190;
                    obj.maxDiffSlew  = 120;
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
                    obj.maxGrad      = 35;
                    obj.maxSlew      = 150;
                    obj.nDummy       = 50;
                    obj.distanceFactorPercentage = 200;
                case 'spiral2d'
                    obj.fov          = 192e-3;
                    obj.Nx           = 192;
                    obj.Ny           = 16;
                    obj.alpha        = 15;
                    obj.thickness    = 3e-3;
                    obj.nSlices      = 12;
                    obj.TE           = 2.5e-3;
                    obj.TR           = 140e-3;
                    obj.readoutTime  = 8e-3;
                    obj.mode         = 'multiShot';
                    obj.maxGrad      = 40;
                    obj.maxSlew      = 150;
                    obj.distanceFactorPercentage = 150;
                    obj.nDummy       = 5;
                otherwise
                    error('SequenceParams_Cima_X: unknown sequence "%s".', seqName)
            end
        end
    end
end