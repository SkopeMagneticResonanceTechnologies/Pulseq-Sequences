classdef SequenceParams_uM_Jupiter_5T
    %SEQUENCEPARAMS_uM_Jupiter_5T  All sequence parameters for the United uM Jupiter 5T.
    %   Called by SequenceParams for imaging sequences (gre2d, epi2d, etc.).
    %   Edit this file to adjust defaults specifically for this scanner.

    methods (Static)
        function obj = apply(obj, seqName, mode) 
            switch lower(seqName)
                case 'gre2d'
                    obj.fov          = 220e-3;
                    obj.Nx           = 128;
                    obj.Ny           = 128;
                    obj.alpha        = 7;
                    obj.thickness    = 3e-3;
                    obj.nSlices      = 12;
                    obj.TE           = [5 10] * 1e-3;
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
                    obj.maxGrad      = 110;
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
                    obj.accFacPE     = 1;
                    obj.addPhaseCorrLines = 1;
                    obj.doPlayFatSat = 0;
                    obj.bFactor      = [0, 1000, 1000, 1000, 1000, 1000, 1000];
                    obj.bDir         = [0,0,0;... % b0
                                        1,0,0;... % x
                                        0,1,0;... % y
                                        0,0,1;... % z
                                        1,1,0;... % xy
                                        0,1,1;... % yz
                                        1,0,1;];  % xz
                    obj.nbValues     = size(obj.bDir,1);
                    obj.seqSpecName  = '';
                    obj.readoutTime  = 0.600e-3;
                    obj.maxGrad      = 110;
                    obj.maxSlew      = 190;
                    obj.maxDiffSlew  = 120;
                    obj.distanceFactorPercentage = 150;               
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
                    error('SequenceParams_uM_Jupiter_5T: unknown sequence "%s".', seqName)
            end
        end
    end
end