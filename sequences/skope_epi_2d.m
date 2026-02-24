classdef skope_epi_2d < PulseqBase
% This is a demoEPI sequence, which includes synchronization scans for
% field-monitoring with a Skope Field Camera and uses ramp-sampling to
% provide an efficient readout. The member method plot() can be used to
% display the generated sequence.
% 
% Notes:
% - The sequence file is written into the current folder.
% - The k-space trajectory during the synchronization scans will not be
%   correctly shown by the member method plot().
% - The x-axis was flipped because of a bug in the Siemens Pulseq 
%   interpreter 1.4.0. 
%
% Example:
%  epi = skope_epi_2d(sequenceParams);
%  epi.plot();
%  epi.test();
%
% See also PulseqBase

% (c) 2026 Skope Magnetic Resonance Technologies AG

    properties        
    
        % A flag to quickly disable phase encoding (1/0) as needed for the delay calibration
        pe_enable = 1             

        % Oversampling factor (in contrast to the product sequence we don't really need it)
        ro_os = 1    

        % Partial Fourier factor: 1: full sampling 0: start with ky=0
        partFourierFactor = 1 

        % Add phase correction lines
        addPhaseCorrLines = true;
           
    end

    properties(SetAccess=protected, GetAccess=public)
        % Spacing of EPI echoes
        echoSpacing
    end    
  
    properties (Access=private)

        % Pulseq transmit object
        rf

        % Pulseq transmit object for SMS
        rfSMS

        % Pulseq fat saturation object
        rf_fs

        % Pulseq ADC event
        adc
        
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

        % z-Blip for SMS
        gzBlip

        % Pulseq spoiling gradient
        gz_fs

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

        % Echo train length
        echoTrainLength

        % Fat shift [Unit: ppm]
        sat_ppm = -3.45;

        % Play out fat saturation pulse
        doPlayFatSat = false;

        % Slice distance factor percentage
        distanceFactorPercentage = 250;

        % Acceleration factor (Phase)
        accFacPE

        % Set trigger output channel
        triggerOutput

        % Frequency offset (Hz) per slice thickness; use as
        % carrier-frequency step when looping over slice groups.
        freqSMS = 0

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

    end

    methods

        function obj = skope_epi_2d(seqParams)

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

            % Set system limits
            obj.sys = mr.opts('MaxGrad', seqParams.maxGrad, ...
                              'GradUnit','mT/m',...
                              'MaxSlew', seqParams.maxSlew, ...
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

            % Set correct axis
            obj.gz.channel =  obj.axesOrder{3};  
            obj.gzReph.channel =  obj.axesOrder{3}; 

            %% Create multi-band pulse
            if obj.multiBandFactor > 1
                % The center of the SMS RF pulse is not accurately
                % determined by the Pulseq simulation
                warning('OFF', 'mr:restoreShape')

                % Input check
                if mod(obj.nSlices, obj.multiBandFactor)
                    error('Number of slices needs to be divisable by multi-band factor.')
                end
                sliceSep = obj.nSlices/obj.multiBandFactor*obj.thickness*(1+obj.distanceFactorPercentage/100);
                [obj.rfSMS, obj.gzSMS, obj.freqSMS, t_rf_center] = CreateSMSPulse(obj.alpha, ...
                                                                obj.thickness, ...
                                                                4, ... timeBwProduct
                                                                8e-3, ... 
                                                                obj.multiBandFactor, ...
                                                                sliceSep, ...
                                                                obj.sys, ...
                                                                'doSim', true, ...              % Plot simulated SMS slice profile
                                                                'type', 'st', ...               % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
                                                                'noRfOffset', false, ...        % don't shift slice (slab) for 3D
                                                                'ftype', 'ls');                 % filter design. 'ls' = least squares
    
                obj.gzSMS.waveform(1) = 0;
                
                % Set correct axis
                obj.gzSMS.channel =  obj.axesOrder{3};  
            end

            %% Define other gradients and ADC events
            deltakx = 1/obj.fov;
            deltaky = 1/obj.fov * obj.accFacPE;
            if obj.multiBandFactor > 1
                deltakz = 1/sliceSep;
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

            %% Create z-blips
            obj.gzBlip = mr.makeTrapezoid(obj.axesOrder{3}, obj.sys, ...
                                'Area', deltakz, ...
                                'Duration', blip_dur);

            %gy = mr.makeTrapezoid(obj.axesOrder{2},lims,'amplitude',deltak/blip_dur*2,'riseTime',blip_dur/2, 'flatTime', 0);
            
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

            %% split the blip into two halves and produce a combined synthetic gradient
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


            %%
            obj.gz_blipPre = mr.makeTrapezoid(obj.axesOrder{3}, obj.sys, ...
                            'Area', deltakz, ...
                            'Duration', mr.calcDuration(obj.gxPre,obj.gyPre,obj.gzReph));
            obj.gz_blipPre.amplitude = obj.gz_blipPre.amplitude*obj.pe_enable;


            %% Create external trigger
            obj.extTrigger = mr.makeDigitalOutputPulse(obj.triggerOutput,'duration', obj.sys.gradRasterTime);

            %% Calculate minimal TE
            if obj.addPhaseCorrLines
                prepareTime = mr.calcDuration(obj.gxPre) + ...
                            3*mr.calcDuration(obj.gx) + ...
                            mr.calcDuration(obj.gyPre);
            else
                prepareTime = mr.calcDuration(obj.gxPre, obj.gyPre);
            end

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

            minTE = sliceTimeTE ...
                  + obj.fillTESMS ...
                  + prepareTime  ...
                  + Ny_pre * mr.calcDuration(obj.gx) ...
                  + mr.calcDuration(obj.gx)/2;
            disp(['Minimal TE is ' num2str((minTE + obj.gradFreeTime)*1000) ' ms'])
            
            obj.fillTE = obj.roundUpToGRT(obj.TE - minTE);
            assert(obj.fillTE >= obj.gradFreeTime, 'Assertion for TE failed');

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

            minTR = sliceTimeTR ...
                  + obj.fillTRSMS ...
                  + obj.fillTE ...
                  + prepareTime ...
                  + obj.echoTrainLength * mr.calcDuration(obj.gx);  
            if obj.doPlayFatSat
                minTR = minTR + mr.calcDuration(obj.gz_fs);
            end

            disp(['Minimal TR is ' num2str(minTR*1000) ' ms'])
             
            obj.fillTR = obj.roundUpToGRT(obj.TR - minTR);
            assert(obj.fillTR >= 0, 'Assertion for TR failed.');

            %% Time from trigger to scanner acquisition
            if obj.addPhaseCorrLines %gxPre and gyPre are split when navigator is ON
                obj.triggerToScannerAcqDelay =  obj.fillTE ...
                                           + mr.calcDuration(obj.gxPre) ...
                                           + mr.calcDuration(obj.gyPre) ...
                                           + obj.adc.delay;  
            else %gxPre and gyPre are played simultaneously when navigator is OFF
                obj.triggerToScannerAcqDelay =  obj.fillTE ...
                                           + mr.calcDuration(obj.gxPre, obj.gyPre) ...
                                           + obj.adc.delay; 
            end

            if obj.addPhaseCorrLines
                obj.triggerToScannerAcqDelay = obj.triggerToScannerAcqDelay + 3*mr.calcDuration(obj.gx);
            end
            
            %% Calculate required camera acquisition duration
            obj.cameraAcqDuration = obj.fillTE ...
                                  + mr.calcDuration(obj.gxPre, obj.gyPre) ...
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
                    obj = runKernel(obj, slc, avg, rep, KernelMode.Sync);
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
                    obj = runKernel(obj, slc, avg, rep, KernelMode.Reference);
                end    
            end

            %% Dummy scans
            for rep=1:obj.nDummy
                for slc = 1:obj.nSlices/obj.multiBandFactor
                    avg = 1;
                    obj = runKernel(obj, slc, avg, rep, KernelMode.Dummy);
                end
            end

            %% Actual imaging sequence
            for rep=1:obj.nRep
                for slc = 1:obj.nSlices/obj.multiBandFactor
                    avg = 1;
                    obj = runKernel(obj, slc, avg, rep, KernelMode.Imaging);
                end
            end

            %% Set the number of imaging triggers
            obj.nTrig = obj.nSlices*obj.nRep;

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
            obj.seq.setDefinition('SliceShifts', obj.slicePositionChronological); 
            obj.seq.setDefinition('readDir_SCT', readDir_SCT);
            obj.seq.setDefinition('phaseDir_SCT', phaseDir_SCT);
            obj.seq.setDefinition('sliceDir_SCT', sliceDir_SCT);    
            obj.seq.setDefinition('SequenceType', 'GRE');

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

            filename = strcat('exports/',string(seqParams.scannerType),'/skope_epi_2d','_',string(obj.sliceOrientation),'_',string(obj.phaseEncDir));           

            if obj.doPlayFatSat == 1
                filename = strcat(filename, '_fs');
            end
      
            filename = strcat(filename, '_R', num2str(obj.accFacPE));  

            if isprop(seqParams, 'seqSpecName') && ~isempty(seqParams.seqSpecName)
                filename = strcat(filename, '_', seqParams.seqSpecName);																				
            end

            obj.seq.write(strcat(filename,'.seq')); 

        end

    end

    methods (Access=private)               
        function obj = runKernel(obj, slc, avg, rep, mode)

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
            if mode==KernelMode.Dummy || mode==KernelMode.Reference || mode==KernelMode.Imaging 
                if obj.doPlayFatSat
                    obj.addBlock(obj.rf_fs, obj.gz_fs);
                end
                
                if obj.multiBandFactor > 1 && (mode==KernelMode.Imaging || mode==KernelMode.Dummy)
                    % Play out the SMS pulse 

                    % Frequency offset (Hz) for SMS slice shift
                    obj.rfSMS.freqOffset = round((slc-1)*obj.freqSMS);

                    % Get the chronological slice index from the slice counter
                    sli = obj.chronologicalSliceSMS(slc);

                    % Excitation pulse and RF spoiling
                    obj.rfSMS.phaseOffset = obj.rf_phase/180*pi - 2*pi*obj.rfSMS.freqOffset * mr.calcRfCenter(obj.rfSMS);  % align the phase for off-center slices
                    obj.adc.phaseOffset = obj.rf_phase/180*pi;
                    obj.addBlock(obj.rfSMS, obj.gzSMS, mr.makeLabel('SET','PMC',false));                    
                else
                    % Slice counter and slice index are identical
                    sli = slc;                    
                    % Play out the standard pulse
                    obj.rf.freqOffset = obj.gz.amplitude * obj.slicePositionChronological(sli);
                     % Compensate for the slice-offset induced phase
                    obj.rf.phaseOffset = obj.rf_phase/180*pi - 2*pi*obj.rf.freqOffset * mr.calcRfCenter(obj.rf); 
                    obj.addBlock(obj.rf, obj.gz, mr.makeLabel('SET','PMC',false));
                    obj.addBlock(obj.gzReph);
                end
            else
                sli = slc; 
                obj.addBlock(obj.gz, mr.makeLabel('SET','PMC',true));
                obj.addBlock(obj.gzReph);
            end

            % Update RF spoiling
            obj.rf_inc = mod(obj.rf_inc+obj.rfSpoilingInc, 360.0);
            obj.rf_phase = mod(obj.rf_phase+obj.rf_inc, 360.0);

            % The standard pulse is shorter than the SMS pulse
            % Include the difference as a delay here for the sync and
            % reference pulses
            if obj.multiBandFactor > 1 
                if mode==KernelMode.Sync || mode==KernelMode.Reference
                    obj.addBlock(mr.makeDelay(obj.fillTESMS));
                end
            end

            if mode==KernelMode.Sync || mode==KernelMode.Imaging || mode==KernelMode.Reference
                obj.addBlock(obj.extTrigger,mr.makeDelay(obj.fillTE));
            else
               obj.addBlock(mr.makeDelay(obj.fillTE)); 
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
                               mr.makeLabel('SET','NAV', false), ...
                               mr.makeLabel('SET','SEG', segment), ...
                               mr.makeLabel('SET','REV', reverse)};
                else
                    labels = {mr.makeLabel('INC','LIN', 1), ...
                              mr.makeLabel('SET','SEG', segment) ...
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
                        if obj.multiBandFactor ==1
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
                        if obj.multiBandFactor ==1
                            obj.addBlock(obj.gx, obj.gy_blipdownup); 
                        else
                            obj.addBlock(obj.gx, obj.gy_blipdownup, mr.scaleGrad(obj.gz_blipdowndown,(-1)^(mod(lin,2))));
                        end
                    end
                end 
                obj.gx.amplitude = -obj.gx.amplitude;   % Reverse polarity of read gradient
            end

            %% TR filling
            if obj.multiBandFactor > 1 && (mode==KernelMode.Sync || mode==KernelMode.Reference)
                    obj.addBlock(mr.makeDelay(obj.fillTR + obj.fillTRSMS));
            else
                obj.addBlock(mr.makeDelay(obj.fillTR));
            end
        
        end
    end

end