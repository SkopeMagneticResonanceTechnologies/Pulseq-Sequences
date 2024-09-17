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

% (c) 2022 Skope Magnetic Resonance Technologies AG

    properties        
    
        % A flag to quickly disable phase encoding (1/0) as needed for the delay calibration
        pe_enable = 1             

        % Oversampling factor (in contrast to the product sequence we don't really need it)
        ro_os = 1    

        % Partial Fourier factor: 1: full sampling 0: start with ky=0
        partFourierFactor = 1 

        % Add phase correction lines
        addPhaseCorrLines = false;
           
    end

    properties(SetAccess=protected, GetAccess=public)
        % Spacing of EPI echoes
        echoSpacing
    end    
  
    properties (Access=private)

        % Pulseq transmit object
        rf

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

        % Pulseq spoiling gradient
        gz_fs

        % Pulseq blip gradient
        gy_blipup
        
        % Pulseq blip gradient
        gy_blipdown

        % Pulseq blip gradient
        gy_blipdownup

        % Pulseq slice refocusing gradient
        gzReph

        % Echo train length
        echoTrainLength

        % Fat shift [Unit: ppm]
        sat_ppm = -3.45;

        doPlayFatSat = false;

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
                     
            % Echo time
            obj.TE = seqParams.TE;

            % Repetition time
            obj.TR = seqParams.TR;

            % Duration of scanner readout event
            obj.readoutTime = seqParams.readoutTime;

            % Flip angle
            obj.alpha = seqParams.alpha;

            % Field of view [Unit: m]
            obj.fov = seqParams.fov;

            % Numbre of readout samples
            obj.Nx = seqParams.Nx;

            % Number of phase encoding steps
            obj.Ny = seqParams.Ny;

            % Slice thickness [Unit: m]
            obj.thickness = seqParams.thickness;

            % Number of slices
            obj.nSlices = seqParams.nSlices;

            % Bug fix for Pulseq error in version 1.4.0.
            obj.doFlipXAxis = seqParams.doFlipXAxis;

            obj.nRep = seqParams.nRep;
            obj.nAve = seqParams.nAve;

            obj.nDummy = seqParams.nDummy;
            obj.doPlayFatSat = seqParams.doPlayFatSat;

            obj.sliceOrientation = seqParams.sliceOrientation;
            obj.phaseEncDir = seqParams.phaseEncDir;

            %% Create a new sequence object
            obj.seq = mr.Sequence(obj.sys);  
            
            %% Time for probe excitation
            obj.gradFreeTime = obj.roundUpToGRT(200e-6);

            %% Axes order
            [obj.axesOrder, obj.axesSign, readDir_SCT, phaseDir_SCT, sliceDir_SCT] ...
                = GetAxesOrderAndSign(obj.sliceOrientation,obj.phaseEncDir, obj.doFlipXAxis);

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

            %% Define other gradients and ADC events
            deltak = 1/obj.fov;
            kWidth = obj.Nx * deltak;
            
            % Phase blip in shortest possible time
            % We round-up the duration to 2x the gradient raster time
            blip_dur = ceil(2*sqrt(deltak/obj.sys.maxSlew)/10e-6/2)*10e-6*2; 

            % The split code below fails if this really makes a trpezoid instead of a triangle.
            % We use negative blips to save one k-space line on our way towards the k-space center
            obj.gy = mr.makeTrapezoid(obj.axesOrder{2}, obj.sys, ...
                                      'Area', -deltak, ...
                                      'Duration', blip_dur); 
            %gy = mr.makeTrapezoid(obj.axesOrder{2},lims,'amplitude',deltak/blip_dur*2,'riseTime',blip_dur/2, 'flatTime', 0);
            
            % readout gradient is a truncated trapezoid with dead times at the beginnig
            % and at the end each equal to a half of blip_dur
            % the area between the blips should be defined by kWidth
            % we do a two-step calculation: we first increase the area assuming maximum
            % slewrate and then scale down the amlitude to fix the area 
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
            adcDwellNyquist = deltak/obj.gx.amplitude/obj.ro_os;

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
            
            % Phase encoding and partial Fourier         
            % PE steps prior to ky=0, excluding the central line
            Ny_pre = round(obj.partFourierFactor*obj.Ny/2-1); 
            
            % PE lines after the k-space center including the central line
            Ny_post = round(obj.Ny/2 + 1);
            obj.echoTrainLength = Ny_pre + Ny_post;
            
            % Pre-phasing gradients
            obj.gxPre = mr.makeTrapezoid(obj.axesOrder{1}, obj.sys, ...
                                         'Area',-obj.gx.area/2);
            obj.gyPre = mr.makeTrapezoid(obj.axesOrder{2}, obj.sys, 'Area', Ny_pre*deltak);

            [obj.gxPre, obj.gyPre] = mr.align('right', obj.gxPre, ...
                                             'left', obj.gyPre);

            % relax the PE prepahser to reduce stimulation
            obj.gyPre = mr.makeTrapezoid(obj.axesOrder{2}, obj.sys, ...
                                         'Area', obj.gyPre.area, ...
                                         'Duration', mr.calcDuration(obj.gxPre,obj.gyPre,obj.gzReph));
            obj.gyPre.amplitude = obj.gyPre.amplitude*obj.pe_enable;

  
            %% Create external trigger
            obj.extTrigger = mr.makeDigitalOutputPulse('ext1','duration', obj.sys.gradRasterTime);

            %% Calculate minimal TE
            if obj.addPhaseCorrLines
                prepareTime = mr.calcDuration(obj.gxPre) + ...
                            3*mr.calcDuration(obj.gx) + ...
                            mr.calcDuration(obj.gyPre);
            else
                prepareTime = mr.calcDuration(obj.gxPre, obj.gyPre);
            end

            minTE = obj.gz.flatTime/2 ...
                  + obj.gz.fallTime ...
                  + mr.calcDuration(obj.gzReph) ...
                  + prepareTime  ...
                  + Ny_pre * mr.calcDuration(obj.gx) ...
                  + mr.calcDuration(obj.gx)/2;
            disp(['Minimal TE is ' num2str((minTE + obj.gradFreeTime)*1000) ' ms'])
            
            obj.fillTE = obj.roundUpToGRT(obj.TE - minTE);
            assert(obj.fillTE >= obj.gradFreeTime, 'Assertion for TE failed');

            %% Calculate minimal TR
            minTR = mr.calcDuration(obj.gz) ...
                  + mr.calcDuration(obj.gzReph) ...
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
            obj.triggerToScannerAcqDelay =  obj.fillTE ...
                                       + mr.calcDuration(obj.gxPre, obj.gyPre) ...
                                       + obj.adc.delay;  

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

            %% Dummy scans
            for rep=1:obj.nDummy
                for slc = 1:obj.nSlices
                    avg = 1;
                    obj = runKernel(obj, slc, avg, rep, KernelMode.Dummy);
                end
            end

            %% Actual imaging sequence
            for rep=1:obj.nRep
                for slc = 1:obj.nSlices
                    avg = 1;
                    obj = runKernel(obj, slc, avg, rep, KernelMode.Imaging);
                end
            end

            %% Set the number of imaging triggers
            obj.nTrig = obj.nSlices;

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
            obj.seq.setDefinition('FOV', [obj.fov obj.fov obj.thickness*obj.nSlices]);
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
            obj.seq.setDefinition('SliceShifts', [obj.thickness*([1:obj.nSlices]-1-(obj.nSlices-1)/2)]); 
            obj.seq.setDefinition('readDir_SCT', readDir_SCT);
            obj.seq.setDefinition('phaseDir_SCT', phaseDir_SCT);
            obj.seq.setDefinition('sliceDir_SCT', sliceDir_SCT);

            %% Write to pulseq file
            if not(isfolder('exports'))
                mkdir('exports')
            end
            obj.seq.write(strcat('exports/skope_epi_2d','_',string(obj.sliceOrientation),'_',string(obj.phaseEncDir),'.seq'));       
            
        end

    end

    methods (Access=private)               
        function obj = runKernel(obj, slc, avg, rep, mode)

            if not(isa(mode, 'KernelMode'))
                error('Expected a kernel mode argument')
            end

            %% RF and ADC settings
            if mode==KernelMode.Dummy || mode==KernelMode.Imaging
                if obj.doPlayFatSat
                    obj.addBlock(obj.rf_fs, obj.gz_fs);
                end
                obj.rf.freqOffset = obj.gz.amplitude * obj.thickness*(slc-1-(obj.nSlices-1)/2);
                 % Compensate for the slice-offset induced phase
                obj.rf.phaseOffset = -2*pi*obj.rf.freqOffset * mr.calcRfCenter(obj.rf); 
                obj.addBlock(obj.rf, obj.gz, mr.makeLabel('SET','PMC',false));
            else
                obj.addBlock(obj.gz, mr.makeLabel('SET','PMC',true));
            end
            
            obj.addBlock(obj.gzReph);

            if mode==KernelMode.Sync || mode==KernelMode.Imaging
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
                           mr.makeLabel('SET','REP', rep-1), ...
                           mr.makeLabel('SET','SLC', slc-1), ...
                           mr.makeLabel('SET','NAV',true)};
                obj.gx.amplitude = -obj.gx.amplitude;
                
                if mode==KernelMode.Sync || mode==KernelMode.Imaging
                    obj.addBlock(obj.gx, labels{:}, obj.adc);
                else
                    obj.addBlock(obj.gx);
                end

                % Second phase correction line
                labels = { mr.makeLabel('SET','LIN', obj.echoTrainLength/2), ...
                           mr.makeLabel('SET','AVG', 1), ...
                           mr.makeLabel('SET','REP', rep-1), ...
                           mr.makeLabel('SET','SLC', slc-1), ...
                           mr.makeLabel('SET','NAV',true)};

                obj.gx.amplitude = -obj.gx.amplitude;
                
                if mode==KernelMode.Sync || mode==KernelMode.Imaging
                    obj.addBlock(obj.gx, labels{:}, obj.adc); 
                else
                    obj.addBlock(obj.gx); 
                end

                % Third phase correction line
                labels = { mr.makeLabel('SET','LIN', obj.echoTrainLength/2), ...
                           mr.makeLabel('SET','AVG', 2), ...
                           mr.makeLabel('SET','REP', rep-1), ...
                           mr.makeLabel('SET','SLC', slc-1), ...
                           mr.makeLabel('SET','NAV',true)};
                obj.gx.amplitude = -obj.gx.amplitude;

                if mode==KernelMode.Sync || mode==KernelMode.Imaging
                    obj.addBlock(obj.gx, labels{:}, obj.adc); 
                else
                    obj.addBlock(obj.gx); 
                end

                % Restore original polarity
                obj.gx.amplitude = -obj.gx.amplitude;
                obj.gxPre.amplitude = -obj.gxPre.amplitude; 

                % Play out phase pre-winding gradient
                obj.addBlock(obj.gyPre);
            else
                obj.addBlock(obj.gxPre, obj.gyPre);
            end

            for lin = 1:obj.echoTrainLength

                % Set labels
                if lin == 1                    
                    labels = { mr.makeLabel('SET','LIN', 0), ...
                               mr.makeLabel('SET','AVG', avg-1), ...
                               mr.makeLabel('SET','REP', rep-1), ...
                               mr.makeLabel('SET','SLC', slc-1), ...
                               mr.makeLabel('SET','NAV', false)};
                else
                    labels = {mr.makeLabel('INC','LIN', 1)};
                end
              
                if lin == 1
                    % Read the first line of k-space with a single half-blip at the end
                    if mode==KernelMode.Sync || mode==KernelMode.Imaging
                        obj.addBlock(obj.gx, obj.gy_blipup, labels{:}, obj.adc); 
                    else
                        obj.addBlock(obj.gx, obj.gy_blipup); 
                    end
                elseif lin==obj.echoTrainLength
                    % Read the last line of k-space with a single half-blip at the beginning
                    if mode==KernelMode.Sync || mode==KernelMode.Imaging
                        obj.addBlock(obj.gx, obj.gy_blipdown, labels{:}, obj.adc); 
                    else
                        obj.addBlock(obj.gx, obj.gy_blipdown); 
                    end
                else
                    % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                    if mode==KernelMode.Sync || mode==KernelMode.Imaging
                        obj.addBlock(obj.gx, obj.gy_blipdownup, labels{:}, obj.adc); 
                    else
                        obj.addBlock(obj.gx, obj.gy_blipdownup); 
                    end
                end 
                obj.gx.amplitude = -obj.gx.amplitude;   % Reverse polarity of read gradient
            end

            %% TR filling
            obj.addBlock(mr.makeDelay(obj.fillTR));
        
        end
    end

end