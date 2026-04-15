classdef skope_se_epi_2d_diff_mb < PulseqBase
% This is a spin echo EPI demo sequence with diffusion encoding gradients, 
% which includes synchronization scans for field-monitoring with a Skope Field Camera 
% and uses ramp-sampling to provide an efficient readout. T
% he member method plot() can be used to display the generated sequence.
% 
% Notes:
% - The sequence file is written into the current folder.
% - The k-space trajectory during the synchronization scans will not be
%   correctly shown by the member method plot().
% - The x-axis was flipped because of a bug in the Siemens Pulseq 
%   interpreter 1.4.0. 
% - dummy scans are played out before every bencoding volume
% - the diffusion encoding allows for:
% i. definition of the b-encoding vector of magnitude
%       e.g.: obj.bFactor=[0, 1000, 1000, 1000];
% ii. definition of the b-encoding vector of directions 
% e.g.: obj.bDir = [0,0,0; 1,0,0; 0,1,1]; % equals to b0, x and y,z cross terms
%
% Example:
%  epi = skope_se_epi_2d_diff(sequenceParams);
%  epi.plot();
%  epi.test();
%
% See also PulseqBase

% (c) 2026 Skope Magnetic Resonance Technologies AG

    properties        
    
        % A flag to quickly disable phase encoding (1/0) as needed for the delay calibration
        pe_enable = 1;             

        % Partial Fourier factor: 1: full sampling 0: start with ky=0
        partFourierFactor = 1; 

        % Add phase correction lines
        addPhaseCorrLines = true;
           
    end

    properties(SetAccess=protected, GetAccess=public)
        % Spacing of EPI echoes
        echoSpacing
    end    
  
    properties (Access=private)
        
        %delay RF90 to RF 180
        delayTE1

        %delay RF 180 to prep readout
        delayTE2 

        % Pulseq transmit object
        rf

        % Pulseq transmit object for SMS
        rfSMS

        % Pulseq fat saturation object
        rf_fs

        % Pulseq pi pulse
        rf180

        % Pulseq transmit object for SMS180
        rfSMS180

        % Pulseq ADC event
        adc

        % Pulseq diffusion gradients
        gDiff
        
        % Pulseq read prewinding gradient
        gxPre

        % Pulseq readout gradient
        gx

        % Pulseq phase prewinding gradient
        gyPre

        % Pulseq phase encoding gradient
        gy

        % Pulseq slice selection gradient
        gz

        % Pulseq slice selection gradient for SMS pulse including refocusing gradient
        gzSMS

        % Pulseq slice refocusing gradient for SMS pulse including refocusing gradient
        gzSMS180

        % z-Blip for SMS
        gzBlip

        % Pulseq spoiling gradient
        gz_fs

        % Pulseq slice selection ref gradient
        gz180

        % Pulseq blip gradient
        gy_blipup
        
        % Pulseq blip gradient
        gy_blipdown

        % Pulseq phase-encoding blip gradient
        gy_blipdownup

        % Pulseq CAIPI blip gradient
        gz_blipup
        gz_blipdown
        gz_blipdowndown
        gz_blipPre

        % Pulseq slice refocusing gradient
        gzReph

        gxSpoil
        gzSpoil

        % Echo train length
        echoTrainLength

        % Fat shift [Unit: ppm]
        sat_ppm = -3.45;

        % Play out fat saturation pulse
        doPlayFatSat = false;

        % Acceleration factor (Phase)
        accFacPE

        % Set trigger output channel
        triggerOutput

        % Slice-select gradient plateau amplitude (Hz/m) of the SMS pulse;
        % use as gzSMSAmplitude * slicePosition to compute RF frequency offsets.
        gzSMSAmplitude = 0
        gzSMSAmplitude180 = 0

        % Spoiling phase of RF and ADC events
        rf_phase = 0
        rf_inc = 0

        % RF spoiling increment (degrees)
        rfSpoilingInc = 117;       

        % Difference between center of SMS-pulse and standard pulse to end
        % of slice refocusing pulse
        fillTESMS = 0
        fillTRSMS = 0

        % Slices indices 
        chronologicalSliceSMS
        anatomicalSliceIndexChronologically
        
        % Diffusion Weighting b-encoding strength
        bFactor
        
        % Diffusion Weighting b-encoding direction (matrix)
        bDir
        
        % number of Diffusion Weighting b-encoding volumes
        nbValues
    end

    methods

        function obj = skope_se_epi_2d_diff_mb(seqParams)

            %% Check input structure
            if not(isa(seqParams,'SequenceParams'))
                error('Input need to be a SequenceParams object.');
            end

            %% Get system limits
            specs = GetMRSystemSpecs(seqParams.scannerType); 

            if not(strcmpi(specs.maxGrad_unit,'mT/m'))
                error('Expected mT/m for maximum gradient.');
            end

            if not(strcmpi(specs.maxSlew_unit,'T/m/s'))
                error('Expected T/m/s for slew rate.');
            end

            %% Check specs
            if seqParams.maxGrad > specs.maxGrad
                error('Scanner does not support requested gradient amplitude.');
            end
            if seqParams.maxSlew > specs.maxSlew
                error('Scanner does not support requested slew rate.');
            end
             if seqParams.maxDiffSlew > specs.maxSlew
                error('Scanner does not support requested Diffusion slew rate.');
            end

            % Set system limits
            obj.sys = mr.opts('MaxGrad', seqParams.maxGrad, ...
                              'GradUnit','mT/m',...
                              'MaxSlew', seqParams.maxSlew, ...
                              'SlewUnit','T/m/s',...
                              'rfRingdownTime', 30e-6, ...
                              'rfDeadtime', 100e-6,...
                              'B0', specs.B0 ...
            );      

            % set Diff grad limits
            obj.sysDiff = mr.opts('MaxGrad', seqParams.maxGrad, ...
                              'GradUnit','mT/m',...
                              'MaxSlew', seqParams.maxDiffSlew, ...
                              'SlewUnit','T/m/s',...
                              'rfRingdownTime', 30e-6, ...
                              'rfDeadtime', 100e-6,...
                              'B0', specs.B0 ...
            ); 

            % Copy all sequence parameters
            fieldNames = fields(seqParams);
            for i = 1:numel(fieldNames)
                fieldname = fieldNames{i};
                if isprop(obj,fieldname)
                    obj.(fieldname) = seqParams.(fieldname);
                end
            end

            %% Check number of repetitions
            if obj.nRep ~= 1
                error('This sequence uses the ONCE flag to mark sync and dummy scans. The number of repetitions can be set on the Sequence Special Card on the scanner.')
            end

            %% Create a new sequence object
            obj.seq = mr.Sequence(obj.sys);  
            
            %% Time for probe excitation
            obj.gradFreeTime = obj.roundUpToGRT(200e-6);

            %% Axes order
            [obj.axesOrder, obj.axesSign, readDir_SCT, phaseDir_SCT, sliceDir_SCT] ...
                = GetAxesOrderAndSign(obj.sliceOrientation,obj.phaseEncDir);

            %% Create fat-sat pulse 
            if obj.doPlayFatSat
                sat_freq = obj.sat_ppm * 1e-6 * obj.sys.B0 * obj.sys.gamma;
                obj.rf_fs = mr.makeGaussPulse(  110*pi/180, ...
                                                'system', obj.sys, ...
                                                'Duration', 8e-3, ...
                                                'dwell', 10e-6,...
                                                'bandwidth', abs(sat_freq), ...
                                                'freqOffset', sat_freq, ...
                                                'use', 'saturation');
    
                % Compensate for the frequency-offset induced phase  
                obj.rf_fs.phaseOffset = -2*pi * obj.rf_fs.freqOffset * mr.calcRfCenter(obj.rf_fs);  
    
                % Spoil up to 0.1mm
                obj.gz_fs = mr.makeTrapezoid(obj.axesOrder{3}, obj.sys, ...
                                             'delay', mr.calcDuration(obj.rf_fs), ...
                                             'Area', 0.1/1e-4); 
            end

            %% Create 90 degree slice selection pulse and gradient
            [obj.rf, obj.gz, obj.gzReph] = mr.makeSincPulse(obj.alpha*pi/180, ...
                                                'system', obj.sys, ...
                                                'Duration',2e-3,...
                                                'SliceThickness', obj.thickness, ...
                                                'apodization', 0.42, ...
                                                'timeBwProduct', 4, ...
                                                'use','excitation');

            %% Create 180 degree slice selection pulse and gradient
            [obj.rf180, obj.gz180] = mr.makeSincPulse(pi, ...
                                                'system', obj.sys, ...
                                                'Duration',4e-3,...
                                                'SliceThickness', obj.thickness, ...
                                                'apodization', 0.5, ...
                                                'timeBwProduct', 4, ...
                                                'use','refocusing');

            % Set correct axis
            obj.gz.channel =  obj.axesOrder{3};  
            obj.gzReph.channel =  obj.axesOrder{3}; 

            %% Create multi-band RF pulse
            if obj.multiBandFactor > 1
                % The center of the SMS RF pulse is not accurately determined by the Pulseq simulation
                warning('OFF', 'mr:restoreShape')

                % Input check
                if mod(obj.nSlices, obj.multiBandFactor)
                    error('Number of slices needs to be divisable by multi-band factor.')
                end
                sliceSep = obj.nSlices/obj.multiBandFactor*obj.thickness*(1+obj.distanceFactorPercentage/100);
                [obj.rfSMS, obj.gzSMS, obj.gzSMSAmplitude, t_rf_center] = CreateSMSPulse(obj.alpha, ...                 % flip angle
                                                                                        obj.thickness, ...              % slice thickness
                                                                                        4, ...                          % timeBWProduct
                                                                                        8e-3, ...                       % duration
                                                                                        obj.multiBandFactor, ...        % MSM factor (1-8)
                                                                                        sliceSep, ...                   % centre-to-centre slice separation (meters)
                                                                                        obj.sys, ...
                                                                                        'doSim', true, ...              % Plot simulated SMS slice profile
                                                                                        'type', 'st', ...               % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
                                                                                        'noRfOffset', false, ...        % don't shift slice (slab) for 3D
                                                                                        'ftype', 'ls');                 % filter design. 'ls' = least squares
                
                % Set correct axis
                obj.gzSMS.channel =  obj.axesOrder{3}; 

                [obj.rfSMS180, obj.gzSMS180, obj.gzSMSAmplitude180, t_rf_center180] = CreateSMSPulse(180, ...           % flip angle
                                                                                        obj.thickness, ...              % slice thickness
                                                                                        4, ...                          % timeBWProduct
                                                                                        28e-3, ...                       % duration
                                                                                        obj.multiBandFactor, ...        % MSM factor (1-8)
                                                                                        sliceSep, ...                   % centre-to-centre slice separation (meters)
                                                                                        obj.sys, ...
                                                                                        'doSim', true, ...              % Plot simulated SMS slice profile
                                                                                        'type', 'inv', ...               % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
                                                                                        'noRfOffset', false, ...        % don't shift slice (slab) for 3D
                                                                                        'ftype', 'ls');                 % filter design. 'ls' = least squares
                
                % Set correct axis
                obj.gzSMS180.channel =  obj.axesOrder{3};
            end

            %% Define other gradients and ADC events
            deltakx = 1/obj.fov;
            deltaky = 1/obj.fov * obj.accFacPE;

            if obj.multiBandFactor > 1
                deltakz = 1/sliceSep/obj.fovShiftFactor;
            else
                deltakz = 0;
            end

            kWidth = obj.Nx * deltakx;
            
            % Phase blip in shortest possible time
            % We round-up the duration to 2x the gradient raster time
            blip_dur = ceil(2*sqrt(deltaky/obj.sys.maxSlew)/10e-6/2)*10e-6*2; 
            blip_dur = max(blip_dur, ceil(2*sqrt(deltakz/obj.sys.maxSlew)/10e-6/2)*10e-6*2);

            % The split code below fails if this really makes a trapezoid instead of a triangle.
            % We use negative blips to save one k-space line on our way towards the k-space center
            obj.gy = mr.makeTrapezoid(obj.axesOrder{2}, obj.sys, ...
                                      'Area', -deltaky, ...
                                      'Duration', blip_dur); 
            
            % Create z-blips
            obj.gzBlip = mr.makeTrapezoid(obj.axesOrder{3}, obj.sys, ...
                                'Area', deltakz, ...
                                'Duration', blip_dur);

            % readout gradient is a truncated trapezoid with dead times at the beginning
            % and at the end each equal to a half of blip_dur
            % the area between the blips should be defined by kWidth
            % we do a two-step calculation: we first increase the area assuming maximum
            % slewrate and then scale down the amplitude to fix the area 
            extra_area = blip_dur/2 * blip_dur/2 * obj.sys.maxSlew; % check unit!;

            obj.gx = mr.makeTrapezoid(obj.axesOrder{1}, obj.sys, ...
                                  'Area', kWidth+extra_area, ...
                                  'duration', obj.readoutTime + blip_dur);

            actual_area = obj.gx.area - obj.gx.amplitude/obj.gx.riseTime * blip_dur/2 * blip_dur/2/2 ...
                        - obj.gx.amplitude/obj.gx.fallTime*blip_dur / 2 * blip_dur/2/2;

            obj.gx.amplitude = obj.gx.amplitude/actual_area*kWidth;
            obj.gx.area = obj.gx.amplitude*(obj.gx.flatTime + obj.gx.riseTime/2 + obj.gx.fallTime/2);
            obj.gx.flatArea = obj.gx.amplitude*obj.gx.flatTime;

            % Calculate ADC
            % we use ramp sampling, so we have to calculate the dwell time and the
            % number of samples, which are will be quite different from Nx and
            % readoutTime/Nx, respectively. 
            adcDwellNyquist = deltakx/obj.gx.amplitude/obj.ro_os;

            % round-down dwell time to 100 ns
            adcDwell = floor(adcDwellNyquist*1e7)*1e-7;

            % on Siemens the number of ADC samples need to be divisible by 4
            adcSamples = floor(obj.readoutTime/adcDwell/4)*4; 

            % MZ: no idea, whether ceil,round or floor is better for the adcSamples...
            obj.adc = mr.makeAdc(adcSamples, ...
                                 'Dwell', adcDwell, ...
                                 'Delay',blip_dur/2);

            % realign the ADC with respect to the gradient
            time_to_center = obj.adc.dwell*((adcSamples-1)/2+0.5);

            % we adjust the delay to align the trajectory with the gradient. We have to align the delay to 1us 
            obj.adc.delay = round((obj.gx.riseTime + obj.gx.flatTime/2-time_to_center)*1e6)*1e-6; 
            
            %% split the blip into two halves and produce a combined synthetic gradient
            gy_parts = mr.splitGradientAt(obj.gy, blip_dur/2, obj.sys);
            [obj.gy_blipup, obj.gy_blipdown,~] = mr.align('right',gy_parts(1),'left',gy_parts(2), obj.gx);
            obj.gy_blipdownup = mr.addGradients({obj.gy_blipdown, obj.gy_blipup}, obj.sys);
            
            % pe_enable support
            obj.gy_blipup.waveform = obj.gy_blipup.waveform * obj.pe_enable;
            obj.gy_blipdown.waveform = obj.gy_blipdown.waveform * obj.pe_enable;
            obj.gy_blipdownup.waveform = obj.gy_blipdownup.waveform * obj.pe_enable;
            
            gz_partsPos = mr.splitGradientAt(obj.gzBlip, blip_dur/2, obj.sys);  
            [obj.gz_blipup, obj.gz_blipdown,~] = mr.align('right',gz_partsPos(1),'left',gz_partsPos(2), obj.gx);
            obj.gz_blipup = mr.scaleGrad(obj.gz_blipup,-1.0);
            obj.gz_blipdown = mr.scaleGrad(obj.gz_blipdown,-1.0);

            gz_partsNeg = mr.splitGradientAt(mr.scaleGrad(obj.gzBlip,-1.0), blip_dur/2, obj.sys);
            [grad1, grad2,~] = mr.align('right',gz_partsNeg(1),'left',gz_partsPos(2), obj.gx);
            obj.gz_blipdowndown = mr.scaleGrad(mr.addGradients({grad1, grad2}, obj.sys),-1.0);
           
            % pe_enable support
            obj.gz_blipup.waveform = obj.gz_blipup.waveform * obj.pe_enable;
            obj.gz_blipdown.waveform = obj.gz_blipdown.waveform * obj.pe_enable;
            obj.gz_blipdowndown.waveform = obj.gz_blipdowndown.waveform * obj.pe_enable;

            %% Phase encoding and partial Fourier         
            % PE steps prior to ky=0, excluding the central line
            Ny_pre = round(obj.partFourierFactor*obj.Ny/2/obj.accFacPE-1); 
            
            % PE lines after the k-space center including the central line
            Ny_post = round(obj.Ny/2/obj.accFacPE + 1);
            obj.echoTrainLength = Ny_pre + Ny_post;
            
            % Pre-phasing gradients
            obj.gxPre = mr.makeTrapezoid(obj.axesOrder{1}, obj.sys, ...
                                         'Area',-obj.gx.area/2);
            obj.gyPre = mr.makeTrapezoid(obj.axesOrder{2}, obj.sys, 'Area', Ny_pre*deltaky);

            [obj.gxPre, obj.gyPre] = mr.align('right', obj.gxPre, ...
                                             'left', obj.gyPre);

            % relax the PE prephaser to reduce stimulation
            obj.gyPre = mr.makeTrapezoid(obj.axesOrder{2}, obj.sys, ...
                                         'Area', obj.gyPre.area, ...
                                         'Duration', mr.calcDuration(obj.gxPre,obj.gyPre,obj.gzReph));
            obj.gyPre.amplitude = obj.gyPre.amplitude*obj.pe_enable;

            obj.gz_blipPre = mr.makeTrapezoid(obj.axesOrder{3}, obj.sys, ...
                            'Area', deltakz, ...
                            'Duration', mr.calcDuration(obj.gxPre,obj.gyPre,obj.gzReph));
            obj.gz_blipPre.amplitude = obj.gz_blipPre.amplitude*obj.pe_enable;

            % gradient spoiling
            obj.gxSpoil = mr.makeTrapezoid(obj.axesOrder{1},'Area', 2*obj.Nx*deltakx,'system', obj.sys);
            obj.gzSpoil = mr.makeTrapezoid(obj.axesOrder{3},'Area',4/obj.thickness,'system', obj.sys);
            spoilTime = 2*mr.calcDuration({obj.gxSpoil,obj.gzSpoil});

            obj.gxSpoil = mr.makeTrapezoid(obj.axesOrder{1},'system', obj.sys,'Duration',spoilTime,'Area',2*obj.Nx*deltakx);
            obj.gzSpoil = mr.makeTrapezoid(obj.axesOrder{3},'system', obj.sys,'Duration',spoilTime,'Area',4/obj.thickness);

            %% Create external trigger
            obj.extTrigger = mr.makeDigitalOutputPulse(obj.triggerOutput,'duration', obj.sys.gradRasterTime);

            %% Calculate minimal TE
            sliceTimeTE = obj.gz.flatTime/2 ...
                        + obj.gz.fallTime ...
                        + mr.calcDuration(obj.gzReph);

            if obj.multiBandFactor > 1
                % Difference between standard and SMS pulse
                % Note that the gradient for the SMS pulse includes the rewinder
                obj.fillTESMS = obj.roundUpToGRT(mr.calcDuration(obj.gzSMS) - t_rf_center - sliceTimeTE);
                assert(obj.fillTESMS > 0, 'SMS pulse is supposed to be longer than normal excitation pulse.')
            else
                obj.fillTESMS = 0;
            end

            if obj.addPhaseCorrLines
                prepareTime = mr.calcDuration(obj.gxPre) + ...
                            3*mr.calcDuration(obj.gx) + ...
                            mr.calcDuration(obj.gyPre);
            else
                prepareTime = mr.calcDuration(obj.gxPre, obj.gyPre);
            end
                   
            TE1 = obj.TE/2;
            TE2 = TE1;

            if obj.multiBandFactor > 1
                obj.delayTE1 = obj.roundUpToGRT(TE1 - (sliceTimeTE ...
                      + obj.fillTESMS ... 
                      + mr.calcDuration(obj.gzSMS180)/2 ));        
                obj.delayTE2 = obj.roundUpToGRT(TE2 - (prepareTime  ...
                      + Ny_pre * mr.calcDuration(obj.gx) ...
                      + mr.calcDuration(obj.gx)/2 ...
                      + mr.calcDuration(obj.gzSMS180)/2 ...
                      + obj.gradFreeTime ));
            else
                obj.delayTE1 = obj.roundUpToGRT(TE1 - (sliceTimeTE ...
                      + obj.fillTESMS ... 
                      + mr.calcDuration(obj.gz180)/2 ));        
                obj.delayTE2 = obj.roundUpToGRT(TE2 - (prepareTime  ...
                      + Ny_pre * mr.calcDuration(obj.gx) ...
                      + mr.calcDuration(obj.gx)/2 ...
                      + mr.calcDuration(obj.gz180)/2 ...
                      + obj.gradFreeTime ));
            end
												  
            %% Calculate minimal TR
            sliceTimeTR = mr.calcDuration(obj.gz) ...
                + mr.calcDuration(obj.gzReph);

            if obj.multiBandFactor > 1
                % Difference between standard and SMS pulse
                % Note that the gradient for the SMS pulse includes the rewinder
                obj.fillTRSMS = mr.calcDuration(obj.gzSMS) - sliceTimeTR;
                assert(obj.fillTRSMS > 0, 'SMS pulse is supposed to be longer than normal excitation pulse.')
            else
                obj.fillTRSMS = 0;
            end
            
            if obj.multiBandFactor > 1
                minTR = sliceTimeTR ...
                      + obj.fillTRSMS ...
                      + obj.delayTE1 + obj.delayTE2 + obj.gradFreeTime ...
                      + mr.calcDuration(obj.gzSMS180) ...
                      + prepareTime ...
                      + obj.echoTrainLength * mr.calcDuration(obj.gx) ...
                      + mr.calcDuration(obj.gxSpoil, obj.gzSpoil); 
            else
                minTR = sliceTimeTR ...
                      + obj.fillTRSMS ...
                      + obj.delayTE1 + obj.delayTE2 + obj.gradFreeTime ...
                      + mr.calcDuration(obj.gz180) ...
                      + prepareTime ...
                      + obj.echoTrainLength * mr.calcDuration(obj.gx) ...
                      + mr.calcDuration(obj.gxSpoil, obj.gzSpoil);
            end

            if obj.doPlayFatSat
                minTR = minTR + mr.calcDuration(obj.gz_fs);
            end

            disp(['Minimal TR is ' num2str(minTR*1000) ' ms'])
             
            obj.fillTR = obj.roundUpToGRT(obj.TR - minTR);
            assert(obj.fillTR >= 0, 'Assertion for TR failed.');

            %% Preparation of diffusion gradients
            for i = 2:seqParams.nbValues %i=1 always b0
                % diffusion weighting calculation
                % delayTE2 is our window for small_delta
                % delayTE1+delayTE2-delayTE2 is our big delta
                % we anticipate that we will use the maximum gradient amplitude, so we need
                % to shorten delayTE2 by gmax/max_sr to accommodate the ramp down 
                bFactor = seqParams.bFactor(i);
                
                for ja = 1:3 %axis index
                    dir = obj.axesOrder{ja}; % 1:x, 2:y, 3:z
                    if seqParams.bDir(i,ja) > 0                                                              
                        small_delta=obj.delayTE2-ceil(obj.sysDiff.maxGrad/obj.sysDiff.maxSlew/obj.sysDiff.gradRasterTime)*obj.sysDiff.gradRasterTime;
                        if obj.multiBandFactor > 1
                            big_delta=obj.delayTE1+mr.calcDuration(obj.rf180,obj.gzSMS180);
                        else
                            big_delta=obj.delayTE1+mr.calcDuration(obj.rf180,obj.gz180);
                        end
                        % we define bFactCalc function below to eventually calculate time-optimal 
                        % gradients. For now we just abuse it with g=1 to give us the coefficient
                        g=sqrt(bFactor*1e6/bFactCalc(1,small_delta,big_delta))*obj.axesSign(ja); 
                        
                        gr=ceil(abs(g)/obj.sysDiff.maxSlew/obj.sysDiff.gradRasterTime)*obj.sysDiff.gradRasterTime;
                        
                        obj.gDiff{i,ja}=mr.makeTrapezoid(dir,'amplitude',g,'riseTime',gr,'flatTime',small_delta-gr,'system',obj.sysDiff);
                        assert(mr.calcDuration(obj.gDiff{i,ja})<=obj.delayTE1);
                        assert(mr.calcDuration(obj.gDiff{i,ja})<=obj.delayTE2);
                    else
                        obj.gDiff{i,ja}=mr.makeTrapezoid(dir,'amplitude',0,'riseTime',gr,'flatTime',small_delta-gr,'system',obj.sysDiff);
                    end
                end
            end
            
            %% Time from trigger to scanner acquisition
            if obj.addPhaseCorrLines %gxPre and gyPre are split when navigator is ON
                obj.triggerToScannerAcqDelay = mr.calcDuration(obj.gxPre) ...
                                           + mr.calcDuration(obj.gyPre) ...
                                           + obj.gradFreeTime ...
                                           + obj.adc.delay;  
            else %gxPre and gyPre are played simultaneously when navigator is OFF
                obj.triggerToScannerAcqDelay = mr.calcDuration(obj.gxPre,obj.gyPre)...
                                           + obj.gradFreeTime ...
                                           + obj.adc.delay;  
            end

            if obj.addPhaseCorrLines
                obj.triggerToScannerAcqDelay = obj.triggerToScannerAcqDelay + 3*mr.calcDuration(obj.gx);
            end
            
            %% Calculate required camera acquisition duration
            obj.cameraAcqDuration = mr.calcDuration(obj.gxPre, obj.gyPre) ...
                                  + obj.echoTrainLength * mr.calcDuration(obj.gx) ...
                                  + 1e-3; % To be safe

            if obj.addPhaseCorrLines
                obj.cameraAcqDuration = obj.cameraAcqDuration + ...
                    3*mr.calcDuration(obj.gx);
            end
            obj.cameraAcqDuration = ceil(obj.cameraAcqDuration*1000)/1000;

            %% Determine chronological order for slice positions

            % Example for 10 slices
            %  Anatomical      Chronological
            %   10              05
            %   09              10
            %   08              04
            %   07              09
            %   06              03
            %   05              08
            %   04              02
            %   03              07
            %   02              01
            %   01              06

            obj.slicePositionAnatomical = [obj.thickness*([1:obj.nSlices]-1-(obj.nSlices-1)/2)]*(1+obj.distanceFactorPercentage/100);

            if mod(obj.nSlices,2) % odd
                sliceOrder = [1:2:obj.nSlices, 2:2:obj.nSlices];
            else
                sliceOrder = [2:2:obj.nSlices, 1:2:obj.nSlices];
            end

            obj.slicePositionChronological = obj.slicePositionAnatomical(sliceOrder);

            obj.anatomicalSliceIndexChronologically = [1:obj.nSlices];
            obj.anatomicalSliceIndexChronologically = obj.anatomicalSliceIndexChronologically(sliceOrder);

            if obj.multiBandFactor > 1

                % Example for 10 slices and MB 2
                %  Anatomical      Chronological
                %   10              05
                %   09              10
                %   08              04
                %   07              09
                %   06              03
                % ----------------------------> Lower ones are measured
                %                              Chronological SMS
                %   05              08          #3
                %   04              02          #5
                %   03              07          #2
                %   02              01          #4
                %   01              06          #1

                np = obj.nSlices/obj.multiBandFactor;

                % Find the lowest nSli/multiBandFactor slices
                for i=1:obj.nSlices/obj.multiBandFactor 
                    lowestChronoSlice(i) = find(i==sliceOrder);
                end
                
                sliceOrderSMS = [1:2:np 2:2:np];
                if ~mod(np,2)
                    % for np = even, change order of last two partitions/shots
                    l = length(sliceOrderSMS);
                    sliceOrderSMS = sliceOrderSMS([1:(l-2) l l-1]);
                end
                obj.chronologicalSliceSMS = lowestChronoSlice(sliceOrderSMS);

            end

            %% Determine the echo spacing
            obj.echoSpacing = mr.calcDuration(obj.gx);
                        
            %% All LABELS / counters an flags are automatically initialized to 0 in the beginning, no need to define initial 0's  
            % so we will just increment LIN after the ADC event (e.g. during the spoiler)
                      
            % Older scanners like Trio may need this dummy delay to keep up
            % with timing
            % obj.addBlock(mr.makeDelay(1)); 
                                                        
            %% Synchronization
            if obj.nSyncDynamics > 0
                for avg = 1:obj.nSyncDynamics
                    slc = 1;
                    rep = 1;
                    obj = runKernel(obj, slc, avg, rep, 1, KernelMode.Sync);
                end
                
                %% Add pause and reset flags
                if obj.preScanPause < 4
                    warning('The pause between the synchronization and imaging scans should be equal or larger than 4 seconds. The current value is okay for simulation purposes.');
                end

                obj.addBlock(mr.makeDelay(obj.preScanPause), mr.makeLabel('SET','LIN', 0), mr.makeLabel('SET','SLC', 0), mr.makeLabel('SET','AVG', 0));
            end
            
            %% Reference scan - Single-band imaging of all slices
                if obj.multiBandFactor > 1
                    for slc = 1:obj.nSlices
                        avg = 1;
                        rep = 1;
                        obj = runKernel(obj, slc, avg, rep, 1, KernelMode.Reference);
                    end    
                end

            %% Main sequence body
            for bValue=1:seqParams.nbValues

                % Dummy scans
                for rep=1:obj.nDummy
                    for slc = 1:obj.nSlices/obj.multiBandFactor
                        avg = 1;
                        obj = runKernel(obj, slc, avg, rep, bValue, KernelMode.Dummy);
                    end
                end
    
                % Actual imaging sequence
                for rep=1:obj.nRep
                    for slc = 1:obj.nSlices/obj.multiBandFactor
                        avg = 1;
                        obj = runKernel(obj, slc, avg, rep, bValue, KernelMode.Imaging);
                    end
                end
            end

            %% Set the number of imaging triggers
            if obj.multiBandFactor == 1
                obj.nTrig = obj.nSlices*obj.nRep*seqParams.nbValues;
            else
                 obj.nTrig = obj.nSlices ...                                                % Reference scans
                            + obj.nSlices/obj.multiBandFactor*obj.nRep*seqParams.nbValues;  % Actual imaging
            end

            %% Calculate Camera Interleave TR (blank time)
            obj.CalculateInterleaveTR(obj.TR);

            %% check whether the timing of the sequence is correct
            [ok, error_report] = obj.seq.checkTiming;
            
            if (ok)
                fprintf('Timing check passed successfully\n');
            else
                fprintf('Timing check failed! Error listing follows:\n');
                fprintf([error_report{:}]);
                fprintf('\n');
            end
            
            %% Prepare sequence export
            obj.seq.setDefinition('FOV', [obj.fov obj.fov obj.thickness*obj.nSlices*(1+obj.distanceFactorPercentage/100)]);
            obj.seq.setDefinition('Name', 'epi2d');
            obj.seq.setDefinition('TE', obj.TE);
            obj.seq.setDefinition('TR', obj.TR);

            %% Parameters needed be added to the scanner data header for trajectory merging
            obj.seq.setDefinition('EchoSpacing', obj.echoSpacing);  
            obj.seq.setDefinition('EchoTrainLength', obj.echoTrainLength); 
            obj.seq.setDefinition('TriggerToScannerAcqDelay', obj.triggerToScannerAcqDelay); 

            %% Parameters to be set on the user interface of the Field Camera            
            % The number of actually acquired dynamics depends on the cameraInterleaveTR.
            obj.seq.setDefinition('CameraNrDynamics', ceil(obj.nTrig/obj.skipFactor));  
            obj.seq.setDefinition('CameraNrSyncDynamics', obj.nSyncDynamics); 
            obj.seq.setDefinition('CameraAcqDuration', obj.cameraAcqDuration);  
            obj.seq.setDefinition('CameraInterleaveTR', obj.cameraInterleaveTR); 
            obj.seq.setDefinition('CameraAqDelay', 0); 
            obj.seq.setDefinition('AdcSampleTime', obj.adc.dwell);    
            obj.seq.setDefinition('Matrix', [obj.Nx obj.Ny]);
            obj.seq.setDefinition('EncodingMatrix', [obj.adc.numSamples obj.echoTrainLength]);
            obj.seq.setDefinition('InplaneAcceleration', obj.accFacPE);
            obj.seq.setDefinition('SliceShifts', obj.slicePositionChronological); 
            obj.seq.setDefinition('readDir_SCT', readDir_SCT);
            obj.seq.setDefinition('phaseDir_SCT', phaseDir_SCT);
            obj.seq.setDefinition('sliceDir_SCT', sliceDir_SCT);
            obj.seq.setDefinition('SequenceType', 'SE');
            obj.seq.setDefinition('SliceOrdering', 'INTERLEAVED');
            if obj.multiBandFactor > 1  
                obj.seq.setDefinition('MultiBandFactor', obj.multiBandFactor);  
                obj.seq.setDefinition('FovShiftFactor', obj.fovShiftFactor); 
            end

            %% Echo spacing check to comply with scanner forbidden bands
             if isfield(specs,'forbiddenBandsEchoSpacingLimits')
                 for i=1:1:size(specs.forbiddenBandsEchoSpacingLimits,2)
                    if obj.echoSpacing>=specs.forbiddenBandsEchoSpacingLimits(i,1) && obj.echoSpacing<=specs.forbiddenBandsEchoSpacingLimits(i,2)
                        error(['Forbidden echo spacing (' num2str(obj.echoSpacing*1000,2) ' ms) for ' seqParams.scannerType ' gradient coil.'])
                    end
                 end
            else
                error('Scanner type does not provide echo spacing limits. Echo spacing might hit scanner forbidden bands. Please provide that information.')
            end

            %% Write to Pulseq file
            if not(isfolder('exports'))
                mkdir('exports')
            end

            if not(isfolder(strcat('exports/',string(seqParams.scannerType))))
                mkdir(strcat('exports/',string(seqParams.scannerType)))
            end

            filename = strcat('exports/',string(seqParams.scannerType),'/skope_se_epi_2d_diff_mb','_',string(obj.sliceOrientation),'_',string(obj.phaseEncDir));           

            if obj.doPlayFatSat == 1
                filename = strcat(filename, '_fs');
            end

            if obj.accFacPE > 1  
                filename = strcat(filename, '_R', num2str(obj.accFacPE));  
            end

            if obj.multiBandFactor > 1    
                filename = strcat(filename, '_MB', num2str(obj.multiBandFactor));
            end

            if isprop(seqParams, 'seqSpecName') && ~isempty(seqParams.seqSpecName)
                filename = strcat(filename, '_', seqParams.seqSpecName);																				
            end

            obj.seq.write(strcat(filename,'.seq')); 
            disp(['Storing sequence file "', char(filename), '.seq"'])

        end

    end

    methods (Access=private)               
        function obj = runKernel(obj, slc, avg, rep, bValue, mode)

            if not(isa(mode, 'KernelMode'))
                error('Expected a kernel mode argument')
            end

             %% Set ONCE-flag to avoid repeating sync and dummy scans
            if mode == KernelMode.Sync || mode==KernelMode.Dummy || mode==KernelMode.Reference
                % ONCE=1 marks the blocks that are only executed in the first repetition
                obj.addBlock(mr.makeLabel('SET','ONCE', 1));
            else
                % Blocks with ONCE=0 are executed on every repetition
                obj.addBlock(mr.makeLabel('SET','ONCE', 0));
            end

            % Flag the single-slice images
            if mode==KernelMode.Reference
                obj.addBlock(mr.makeLabel('SET','SMS', true));
            else
                obj.addBlock(mr.makeLabel('SET','SMS', false));
            end


            %% RF and ADC settings
            if mode==KernelMode.Dummy || mode==KernelMode.Imaging || mode==KernelMode.Reference
                if obj.doPlayFatSat
                    obj.addBlock(obj.rf_fs, obj.gz_fs);
                end

                if obj.multiBandFactor > 1 && (mode==KernelMode.Imaging || mode==KernelMode.Dummy)
                    % Play out multi-band pulse 
                    % Get the chronological slice index from the slice counter
                    sli = obj.chronologicalSliceSMS(slc);

                    % Compensate for the slice-offset induced phase
                    obj.rfSMS.freqOffset = obj.gzSMSAmplitude * obj.slicePositionChronological(sli);
                      
                    % Excitation pulse and RF spoiling
                    obj.rfSMS.phaseOffset = obj.rf_phase/180*pi - 2*pi*obj.rfSMS.freqOffset * mr.calcRfCenter(obj.rfSMS);  % align the phase for off-center slices
                    obj.adc.phaseOffset = obj.rf_phase/180*pi;
                    obj.addBlock(obj.rfSMS, obj.gzSMS, mr.makeLabel('SET','PMC',false));            

                    % For debugging
                    % disp([char(mode) ' SlC ' num2str(slc) ' SLI '  num2str(sli) ' POS ' num2str(obj.slicePositionChronological(sli)*1000) ' mm' ])
                else
                    % Play out single-band pulse 
                    % Slice counter and slice index are identical
                    sli = slc;                    
                    % Play out the standard pulse
                    obj.rf.freqOffset = obj.gz.amplitude * obj.slicePositionChronological(sli);
                    % Compensate for the slice-offset induced phase
                    obj.rf.phaseOffset = obj.rf_phase/180*pi - 2*pi*obj.rf.freqOffset * mr.calcRfCenter(obj.rf); 
                    obj.adc.phaseOffset = obj.rf_phase/180*pi;
                    obj.addBlock(obj.rf, obj.gz, mr.makeLabel('SET','PMC',false));
                    obj.addBlock(obj.gzReph);
                end
            else
                % if obj.doPlayFatSat
                %     obj.addBlock(obj.gz_fs);
                % end
                sli = slc;
                obj.addBlock(obj.gz, mr.makeLabel('SET','PMC',true));
                obj.addBlock(obj.gzReph);
            end  

            % Update RF spoiling (RR: not really used...)
            obj.rf_inc = mod(obj.rf_inc+obj.rfSpoilingInc, 360.0);
            obj.rf_phase = mod(obj.rf_phase+obj.rf_inc, 360.0);

            % The standard pulse is shorter than the SMS pulse
            % Include the difference as a delay here for the sync and
            % reference pulses
            % if obj.multiBandFactor > 1 
            %     if mode==KernelMode.Sync || mode==KernelMode.Reference					  																							
            %         obj.addBlock(mr.makeDelay(obj.fillTESMS));
            %     end
            % end
            
            if bValue<2 % b0
                obj.addBlock(mr.makeDelay(obj.delayTE1));
            else % nonzero b-encoding 
                obj.addBlock(mr.makeDelay(obj.delayTE1-mr.calcDuration(obj.gDiff{2,1}))); 
                obj.addBlock(obj.gDiff{bValue,1}, obj.gDiff{bValue,2}, obj.gDiff{bValue,3});
            end

            if mode==KernelMode.Dummy || mode==KernelMode.Imaging || mode==KernelMode.Reference   
                if obj.multiBandFactor > 1 && (mode==KernelMode.Imaging || mode==KernelMode.Dummy)
                    % Play out multi-band pulse 
                    % Get the chronological slice index from the slice counter
                    sli = obj.chronologicalSliceSMS(slc);

                    % Compensate for the slice-offset induced phase
                    obj.rfSMS180.freqOffset = obj.gzSMSAmplitude180 * obj.slicePositionChronological(sli);
                      
                    % Excitation pulse and RF spoiling
                    obj.rfSMS180.phaseOffset = obj.rf_phase/180*pi - 2*pi*obj.rfSMS180.freqOffset * mr.calcRfCenter(obj.rfSMS180);  % align the phase for off-center slices
                    obj.addBlock(obj.rfSMS180, obj.gzSMS180);       
                else
                    sli = slc;  
                    obj.rf180.freqOffset=obj.gz180.amplitude * obj.slicePositionChronological(sli);                     
                    obj.rf180.phaseOffset = obj.rf_phase/180*pi - 2*pi*obj.rf180.freqOffset * mr.calcRfCenter(obj.rf180);                 
                    obj.addBlock(obj.rf180, obj.gz180);
                end
            else
                obj.addBlock(obj.gz180);
            end
            
            if bValue<2 % b0
                obj.addBlock(mr.makeDelay(obj.delayTE2));
            else % nonzero b-encoding 
                obj.addBlock(obj.gDiff{bValue,1}, obj.gDiff{bValue,2}, obj.gDiff{bValue,3});
                obj.addBlock(mr.makeDelay(obj.delayTE2-mr.calcDuration(obj.gDiff{2,1})));                
            end

            if mode==KernelMode.Sync || mode==KernelMode.Imaging || mode==KernelMode.Reference 
                obj.addBlock(obj.extTrigger,mr.makeDelay(obj.gradFreeTime)); 												   
            else
                obj.addBlock(mr.makeDelay(obj.gradFreeTime));
            end
                      
            if obj.addPhaseCorrLines
               
                % Start with flip gx amplitude
                obj.gxPre.amplitude = -obj.gxPre.amplitude;   
                obj.addBlock(obj.gxPre); 

                % First phase correction line
                labels = { mr.makeLabel('SET','LIN', obj.echoTrainLength/2), ...
                           mr.makeLabel('SET','AVG', 0), ...
                           mr.makeLabel('SET','SEG', 1), ...
                           mr.makeLabel('SET','REP', rep-1), ...
                           mr.makeLabel('SET','SLC', sli-1), ...
                           mr.makeLabel('SET','SET', bValue-1), ...
                           mr.makeLabel('SET','NAV',true)};
                obj.gx.amplitude = -obj.gx.amplitude;
                
                if mode==KernelMode.Sync || mode==KernelMode.Imaging || mode==KernelMode.Reference
                    obj.addBlock(obj.gx, labels{:}, obj.adc);
                else
                    obj.addBlock(obj.gx);
                end

                % Second phase correction line
                labels = { mr.makeLabel('SET','LIN', obj.echoTrainLength/2), ...
                           mr.makeLabel('SET','AVG', 0), ...
                           mr.makeLabel('SET','SEG', 0), ...
                           mr.makeLabel('SET','REP', rep-1), ...
                           mr.makeLabel('SET','SLC', sli-1), ...
                           mr.makeLabel('SET','SET', bValue-1), ...
                           mr.makeLabel('SET','NAV',true)};

                obj.gx.amplitude = -obj.gx.amplitude;
                
                if mode==KernelMode.Sync || mode==KernelMode.Imaging || mode==KernelMode.Reference
                    obj.addBlock(obj.gx, labels{:}, obj.adc); 
                else
                    obj.addBlock(obj.gx); 
                end

                % Third phase correction line
                labels = { mr.makeLabel('SET','LIN', obj.echoTrainLength/2), ...
                           mr.makeLabel('SET','AVG', 1), ...
                           mr.makeLabel('SET','SEG', 1), ...
                           mr.makeLabel('SET','REP', rep-1), ...
                           mr.makeLabel('SET','SLC', sli-1), ...
                           mr.makeLabel('SET','SET', bValue-1), ...
                           mr.makeLabel('SET','NAV',true)};
                obj.gx.amplitude = -obj.gx.amplitude;

                if mode==KernelMode.Sync || mode==KernelMode.Imaging || mode==KernelMode.Reference
                    obj.addBlock(obj.gx, labels{:}, obj.adc); 
                else
                    obj.addBlock(obj.gx); 
                end

                % Restore original polarity
                obj.gx.amplitude = -obj.gx.amplitude;
                obj.gxPre.amplitude = -obj.gxPre.amplitude; 

                % Play out phase pre-winding gradient
                if obj.multiBandFactor == 1
                    obj.addBlock(obj.gyPre);
                else
                    if mode == KernelMode.Imaging || mode == KernelMode.Dummy
                        obj.addBlock(obj.gyPre, obj.gz_blipPre);
                    else
                        obj.addBlock(obj.gyPre);
                    end
                end
            else
                if obj.multiBandFactor == 1
                    obj.addBlock(obj.gxPre, obj.gyPre);
                else
                    if mode == KernelMode.Imaging || mode == KernelMode.Dummy
                        obj.addBlock(obj.gxPre, obj.gyPre, obj.gz_blipPre);
                    else
                        obj.addBlock(obj.gxPre, obj.gyPre);
                    end
                end
            end

            for lin = 1:obj.echoTrainLength
                
                if mod(lin,2) % odd line
                    segment = 0;
                    reverse = false;
                else % even line
                    segment = 1;
                    reverse = true;
                end

                % Set labels
                if lin == 1                    
                    labels = { mr.makeLabel('SET','LIN', 0), ...
                               mr.makeLabel('SET','AVG', avg-1), ...
                               mr.makeLabel('SET','REP', rep-1), ...
                               mr.makeLabel('SET','SLC', sli-1), ...
                               mr.makeLabel('SET','SET', bValue-1), ...
                               mr.makeLabel('SET','NAV', false), ...
                               mr.makeLabel('SET','SEG', segment), ...
                               mr.makeLabel('SET','REV', reverse)};
                else
                    labels = {mr.makeLabel('INC','LIN', 1), ...
                               mr.makeLabel('SET','SEG', segment), ...
                               mr.makeLabel('SET','REV', reverse)};
                end

                if lin == 1
                    % Read the first line of k-space with a single half-blip at the end
                    if mode==KernelMode.Imaging
                        if obj.multiBandFactor ==1
                            obj.addBlock(obj.gx, obj.gy_blipup, labels{:}, obj.adc);
                        else
                            obj.addBlock(obj.gx, obj.gy_blipup, obj.gz_blipup, labels{:}, obj.adc);
                        end
                    elseif mode==KernelMode.Sync || mode==KernelMode.Reference 
                        obj.addBlock(obj.gx, obj.gy_blipup, labels{:}, obj.adc);
                    else % Dummy - no ADC
                        if obj.multiBandFactor ==1
                            obj.addBlock(obj.gx, obj.gy_blipup); 
                        else
                            obj.addBlock(obj.gx, obj.gy_blipup, obj.gz_blipup); 
                        end
                    end
                    
                elseif lin==obj.echoTrainLength
                    % Read the last line of k-space with a single half-blip at the beginning
                    if mode==KernelMode.Imaging  
                        if obj.multiBandFactor == 1
                            obj.addBlock(obj.gx, obj.gy_blipdown, labels{:}, obj.adc); 
                        else
                            if mod(lin,2) % odd
                                obj.addBlock(obj.gx, obj.gy_blipdown, mr.scaleGrad(obj.gz_blipdown,-1), labels{:}, obj.adc);
                            else
                                obj.addBlock(obj.gx, obj.gy_blipdown, obj.gz_blipdown, labels{:}, obj.adc);
                            end
                        end
                    elseif mode==KernelMode.Sync || mode==KernelMode.Reference 
                        obj.addBlock(obj.gx, obj.gy_blipdown, labels{:}, obj.adc); 
                    else
                        if obj.multiBandFactor == 1
                            obj.addBlock(obj.gx, obj.gy_blipdown); 
                        else
                            if mod(lin,2) % odd
                                obj.addBlock(obj.gx, obj.gy_blipdown, mr.scaleGrad(obj.gz_blipdown,-1)); 
                            else
                                obj.addBlock(obj.gx, obj.gy_blipdown, obj.gz_blipdown); 
                            end
                        end
                    end

                else
                   % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                    if mode==KernelMode.Imaging 
                        if obj.multiBandFactor == 1
                            obj.addBlock(obj.gx, obj.gy_blipdownup, labels{:}, obj.adc); 
                        else
                            obj.addBlock(obj.gx, obj.gy_blipdownup, mr.scaleGrad(obj.gz_blipdowndown,(-1)^(mod(lin,2))), labels{:}, obj.adc); 
                        end
                    elseif mode==KernelMode.Sync ||  mode==KernelMode.Reference 
                        obj.addBlock(obj.gx, obj.gy_blipdownup, labels{:}, obj.adc); 
                    else
                        if obj.multiBandFactor == 1
                            obj.addBlock(obj.gx, obj.gy_blipdownup); 
                        else
                            obj.addBlock(obj.gx, obj.gy_blipdownup, mr.scaleGrad(obj.gz_blipdowndown,(-1)^(mod(lin,2))));
                        end
                    end
                end 
                obj.gx.amplitude = -obj.gx.amplitude;   % Reverse polarity of read gradient
            end

            %% Spoiling
            spoilBlockContents = {obj.gxSpoil, obj.gzSpoil};
            obj.addBlock(spoilBlockContents{:});

            %% TR filling
            obj.addBlock(mr.makeDelay(obj.fillTR));
        
        end
    end

end

function b=bFactCalc(g, delta, DELTA)
    % see DAVY SINNAEVE Concepts in Magnetic Resonance Part A, Vol. 40A(2) 39–65 (2012) DOI 10.1002/cmr.a
    % b = gamma^2  g^2 delta^2 sigma^2 (DELTA + 2 (kappa - lambda) delta)
    % in pulseq we don't need gamma as our gradinets are Hz/m
    % however, we do need 2pi as diffusion equations are all based on phase
    % for rect gradients: sigma=1 lambda=1/2 kappa=1/3 
    % for trapezoid gradients: TODO
    sigma=1;
    %lambda=1/2;
    %kappa=1/3;
    kappa_minus_lambda=1/3-1/2;
    b= (2*pi * g * delta * sigma)^2 * (DELTA + 2*kappa_minus_lambda*delta);
end