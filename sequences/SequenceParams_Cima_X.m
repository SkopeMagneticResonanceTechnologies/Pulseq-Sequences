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
                    obj.TE           = 71e-3;
                    obj.TR           = 117e-3;
                    obj.alpha        = 90;
                    obj.fov          = 220e-3;
                    obj.Nx           = 180;
                    obj.Ny           = 180;
                    obj.thickness    = 3e-3;
                    obj.nSlices      = 12;
                    obj.nDummy       = 1;   % no trigger played out. Total number of readout dummies = nDummy*nSlices.
                    obj.nRep         = 1;
                    obj.accFacPE     = 2;
                    obj.partFourierFactor = 7/8;
                    obj.addPhaseCorrLines = true;
                    obj.doPlayFatSat = true;
                    obj.bFactor = 2000 * [0, 1, 1, 1, 1, 1, 1];
                    obj.bDir         =  [0,0,0;... % b0
                                        1,0,0;... % x
                                        0,1,0;... % y
                                        0,0,1;... % z
                                        [1,1,0]./sqrt(2);... % xy
                                        [0,1,1]./sqrt(2);... % yz
                                        [1,0,1]./sqrt(2);];  % xz
                    obj.nbValues     = size(obj.bDir,1);
                    obj.seqSpecName  = '';
                    obj.readoutTime  = 0.650e-3;
                    obj.maxGrad      = 190;
                    obj.maxSlew      = 190;
                    %failed 65/180, try 40/180 and 60/160
                    obj.maxDiffSlew  = 63; %60  
                    obj.maxDiffGrad  = 172;
                    obj.distanceFactorPercentage = 150;
                    obj.gradFreeTime = 300e-6;
                case 'se_spiral_2d_diff'
                    obj.fov          = 192e-3;
                    obj.Nx           = 192;
                    obj.Ny           = 16;
                    obj.alpha        = 90;
                    obj.thickness    = 3e-3;
                    obj.nSlices      = 12;
                    obj.TE           = 2.5e-3;
                    obj.TR           = 140e-3;
                    obj.readoutTime  = 8e-3;
                    obj.mode         = 'multiShot';
                    obj.maxGrad      = 190;
                    obj.maxSlew      = 190;
                    obj.distanceFactorPercentage = 150;
                    obj.nDummy       = 1;
                    obj.accFac       = 1;
                    obj.doPlayFatSat = true;
                    obj.bFactor = 2000 * [0, 1, 1, 1, 1, 1, 1];
                    obj.bDir         =  [0,0,0;... % b0
                                        1,0,0;... % x
                                        0,1,0;... % y
                                        0,0,1;... % z
                                        [1,1,0]./sqrt(2);... % xy
                                        [0,1,1]./sqrt(2);... % yz
                                        [1,0,1]./sqrt(2);];  % xz
                    obj.nbValues     = size(obj.bDir,1);
                    obj.maxDiffSlew  = 63;
                    obj.maxDiffGrad  = 172;
                    obj.gradFreeTime = 300e-6;
                case 'gre3d'
                    obj.fov          = [0.22 0.22 0.22];
                    obj.Nx           = 112;
                    obj.Ny           = 112;
                    obj.Nz           = 112;
                    obj.alpha        = 1;
                    obj.TE           = [5.2 9.7 14.2 18.7 23.2 27.7 32.2 36.7]*1e-3;
                    obj.TR           = 55e-3;
                    obj.seqSpecName  = '8echoes';
                    % obj.TE           = [5.2 9.7]*1e-3;
                    % obj.TR           = 25e-3;
                    % obj.seqSpecName  = '2echoes';
                    obj.readoutTime  = 20e-6 * obj.Nx;
                    obj.maxGrad      = 190;
                    obj.maxSlew      = 190;
                    obj.nDummy       = 1;               %Total number of readout dummies = nDummy*Ny
                    obj.distanceFactorPercentage = 200;
                    obj.accFacPE     = 2;
                    obj.accFacZ      = 1;
                    obj.doPlayFatSat = true;
                    obj.gradFreeTime = 200e-6;
                otherwise
                    error('SequenceParams_Cima_X: unknown sequence "%s".', seqName)
            end
        end
    end
end