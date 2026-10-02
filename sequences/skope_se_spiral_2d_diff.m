classdef skope_se_spiral_2d_diff < PulseqBase
% This is a spin echo spiral demo sequence with diffusion encoding gradients, 
% which includes synchronization scans for field-monitoring with a Skope Field Camera 
% . The member method plot() can be used to display the generated sequence.
% 
% Notes:
% - The sequence is not adaptive, i.e. the parameters can't be changed, as
%   the spiral waveform is taken as argument
% - The k-space trajectory during the synchronization scans will not be
%   correctly shown by the member method plot().
% - The x-axis was flipped because of a bug in the Siemens Pulseq 
%   interpreter 1.4.0. 
% - dummy scans are played out for no b-encoding
% - the diffusion encoding allows for:
% i. definition of the b-encoding vector of magnitude
%       e.g.: obj.bFactor=[0, 1000, 1000, 1000];
% ii. definition of the b-encoding vector of directions 
% e.g.: obj.bDir = [0,0,0; 1,0,0; 0,1,1]; % equals to b0, x and y,z cross terms
%
% Example:
%  spiral = skope_se_spiral_2d_diff(sequenceParams, waveform);
%  spiral.plot();
%  spiral.test();
%
% See also PulseqBase

% (c) 2026 Skope Magnetic Resonance Technologies AG

    properties        

        % Trigger output channel
        triggerOutput;
           
    end  
  
    properties (Access=private)

        %delay RF90 to RF 180
        delayTE1

        %delay RF 180 to prep readout
        delayTE2 

        % Pulseq transmit object
        rf

        % Pulseq fat saturation object
        rf_fs

        % Pulseq pi pulse
        rf180

        % Pulseq ADC event
        adc

        % Spiral gradient (loaded)
        gspiral

        % Pulseq diffusion gradients
        gDiff

        % Pulseq spiral readout gradient (x)
        gx

        % Pulseq spiral readout gradient (y)
        gy

        % Pulseq slice selection gradient
        gz

        % Pulseq spoiling gradient
        gz_fs

        % Pulseq slice selection ref gradient
        gz180

        % ADC dwelltime
        adcDwelltime

        % Pulseq slice refocusing gradient
        gzReph

        % Refocussing gradients for spiral
        gxRefoc
        gyRefoc
        
        % Fat suppression gradients
        gx_fs
        gy_fs
        gx_fs_pre
        gy_fs_pre
        gz_fs_pre

        % Fat shift [Unit: ppm]
        sat_ppm = -3.45;

        % Play out fat saturation pulse
        doPlayFatSat = false; 

        bFactor

        bDir

        nbValues

        mode = 'default'
    end

    methods

        function obj = skope_se_spiral_2d_diff(seqParams, waveform)

            warning('OFF', 'mr:restoreShape')

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
            obj.sys = mr.opts( ...
                            'MaxGrad', seqParams.maxGrad, ...
                            'GradUnit','mT/m', ...
                            'MaxSlew', seqParams.maxSlew, ...
                            'SlewUnit','T/m/s', ...
                            'rfRingdownTime',30e-6, ...
                            'rfDeadtime',100e-6, ...
                            'adcDeadTime',10e-6, ...
                            'adcSamplesLimit',specs.adcSamplesLimit, ...
                            'adcSamplesDivisor',specs.adcSamplesDivisor, ...
                            'B0',specs.B0);

            % set Diff grad limits
            obj.sysDiff = mr.opts('MaxGrad', seqParams.maxDiffGrad, ...
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

            % ADC dwelltime
            obj.adcDwelltime = 2e-6;

            %% Check number of repetitions
            if obj.nRep ~= 1
                error('This sequence uses the ONCE flag to mark sync and dummy scans. The number of repetitions can be set on the Sequence Special Card on the scanner.')
            end

            %% Create a new sequence object
            obj.seq = mr.Sequence(obj.sys);  
            
            %% Time for probe excitation
            obj.gradFreeTime = obj.roundUpToGRT(300e-6);

            %% Axes order
            [obj.axesOrder, obj.axesSign, readDir_SCT, phaseDir_SCT, sliceDir_SCT] ...
                = GetAxesOrderAndSign(obj.sliceOrientation,obj.phaseEncDir);

            %% Create fat-sat pulse 
            if obj.doPlayFatSat
                sat_freq = obj.sat_ppm * 1e-6 * obj.sys.B0 * obj.sys.gamma;
                obj.rf_fs = mr.makeGaussPulse(  110*pi/180, ...
                                                'system', obj.sys, ...
                                                'Duration', 5e-3, ... % 8
                                                'dwell', 10e-6,...
                                                'bandwidth', abs(sat_freq), ...
                                                'freqOffset', sat_freq, ...
                                                'use', 'saturation');
    
                % Compensate for the frequency-offset induced phase  
                obj.rf_fs.phaseOffset = -2*pi * obj.rf_fs.freqOffset * mr.calcRfCenter(obj.rf_fs);  
    
                % Spoil up to 0.1mm
                % with a limited amplitude 8mT/m
                g_fs = mr.convert(8,'mT/m','Hz/m','gamma',42576000);

                obj.gz_fs_pre = mr.makeTrapezoid(obj.axesOrder{3}, obj.sys, ...
                                             'Area', -0.1/1e-4, ...
                                             'maxGrad', g_fs); 
                obj.gx_fs_pre = mr.makeTrapezoid(obj.axesOrder{1}, obj.sys, ...
                                             'Area', -0.1/1e-4, ...
                                             'maxGrad', g_fs); 
                obj.gy_fs_pre = mr.makeTrapezoid(obj.axesOrder{2}, obj.sys, ...
                                             'Area', -0.1/1e-4, ...
                                             'maxGrad', g_fs);  

                obj.gz_fs = mr.makeTrapezoid(obj.axesOrder{3}, obj.sys, ...
                                             'delay', mr.calcDuration(obj.rf_fs), ...
                                             'Area', 0.1/1e-4, ...
                                             'maxGrad', g_fs); 
                obj.gx_fs = mr.makeTrapezoid(obj.axesOrder{1}, obj.sys, ...
                                             'delay', mr.calcDuration(obj.rf_fs), ...
                                             'Area', 0.1/1e-4, ...
                                             'maxGrad', g_fs); 
                obj.gy_fs = mr.makeTrapezoid(obj.axesOrder{2}, obj.sys, ...
                                             'delay', mr.calcDuration(obj.rf_fs), ...
                                             'Area', 0.1/1e-4, ...
                                             'maxGrad', g_fs); 

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

            %% Determine resolution
            % Convert waveform [mT/m] to gradient in [Hz/m] and integrate
            gspiral_tmp = waveform.' * obj.sys.gamma / 1000; % [Hz/m]
            k_tmp = cumsum(gspiral_tmp, 2) * obj.sys.gradRasterTime; % k-space trajectory [1/m]
            kmax = max(sqrt(k_tmp(1,:).^2 + k_tmp(2,:).^2)); % max k-space radius [1/m]
            resolution = 1 / (2 * kmax); % spatial resolution [m]
            n = round(obj.fov/resolution);
            n = n - mod(n,2);
            
            %% Define other gradients and ADC events
            % Create spiral waveform
            obj.gspiral = waveform.' * obj.sys.gamma / 1000;
            obj.gx = mr.makeArbitraryGrad(obj.axesOrder{1},obj.gspiral(1,:), ...
                'system', obj.sys, ...
                'first', 0, ...
                'last', 0);
            obj.gy = mr.makeArbitraryGrad(obj.axesOrder{2},obj.gspiral(2,:), ...
                'system', obj.sys, ...
                'first', 0, ...
                'last', 0);
            durADC = mr.calcDuration(obj.gx);

            % Let's make the number of samples divisible by 8 and 10
            nSamplesADC = round(durADC/obj.adcDwelltime/80)*80;
            adcSamplesPerSegment = nSamplesADC;

            if nSamplesADC > obj.sys.adcSamplesLimit
                [adcSegments,adcSamplesPerSegment] = mr.calcAdcSeg(nSamplesADC, ...
                                                    obj.adcDwelltime, ...
                                                    obj.sys, ...
                                                    'shorten');
                nSamplesADC = adcSegments*adcSamplesPerSegment;
            end
            
            % Update duration
            durADC = nSamplesADC*obj.adcDwelltime;

            obj.adc = mr.makeAdc(nSamplesADC, ...
                                'Duration', durADC, ...
                                'system', obj.sys);

            % Compute moment
            mx = -trapz(obj.gspiral(1,:))*obj.sys.gradRasterTime;
            my = -trapz(obj.gspiral(2,:))*obj.sys.gradRasterTime;
           
            % gradient spoiling and refocussing
            obj.gxRefoc = mr.makeTrapezoid(obj.axesOrder{1},'Area',mx,'system', obj.sys);
            obj.gyRefoc = mr.makeTrapezoid(obj.axesOrder{2},'Area',my,'system', obj.sys);            
            refocTime = max([mr.calcDuration(obj.gxRefoc),mr.calcDuration(obj.gyRefoc)]);

            % Stretch gradients
            obj.gxRefoc = mr.makeTrapezoid(obj.axesOrder{1},'Area',mx,'system', obj.sys, 'Duration',refocTime);
            obj.gyRefoc = mr.makeTrapezoid(obj.axesOrder{2},'Area',my,'system', obj.sys, 'Duration',refocTime);

            %% Create external trigger
            obj.extTrigger = mr.makeDigitalOutputPulse(obj.triggerOutput,'duration', obj.sys.gradRasterTime);

            %% Calculate minimal TE
            TE1 = obj.TE/2;
            TE2 = TE1;
            obj.delayTE1 = obj.roundUpToGRT(TE1 - (obj.gz.flatTime/2 ...
                  + obj.gz.fallTime ...
                  + mr.calcDuration(obj.gzReph) ...
                  + mr.calcDuration(obj.gz180)/2 ));         

            obj.delayTE2 = obj.roundUpToGRT(TE2 - (...
                  + mr.calcDuration(obj.gz180)/2 ...
                  + obj.gradFreeTime ));

            %% Calculate minimal TR
            minTR = mr.calcDuration(obj.gz) ...
                  + mr.calcDuration(obj.gzReph) ...
                  + obj.delayTE1 + obj.delayTE2 + obj.gradFreeTime ...
                  + mr.calcDuration(obj.gz180) ...
                  + mr.calcDuration(obj.gx) ...
                  + refocTime;
            if obj.doPlayFatSat
                minTR = minTR + mr.calcDuration(obj.gz_fs) + mr.calcDuration(obj.gz_fs_pre);
            end

            disp(['Minimal TR is ' num2str(minTR*1000) ' ms'])
             
            obj.fillTR = obj.roundUpToGRT(obj.TR - minTR);
            assert(obj.fillTR >= 0, 'Assertion for TR failed.');

            %% Preparation of diffusion gradients            
            %        90° RF                           180° RF
            %          |                                |
            %          |        ______       _          |       ______
            %                  /      \      ¦          |      /      \
            %          |      /        \     ¦ G        |     /        
            %          |     /          \    ¦          |    /          \
            %          |____/            \___¦__________|___/            \____
            %               <tau> 
            %               <  delta  >
            %               <             Delta            >    

            for i = 1:seqParams.nbValues %i=1 always b0
                bFactor = seqParams.bFactor(i);
                
                if i<2 %run this once
                    % estimate timing for diffusion gradients when bFactor is max
                    bFactor_max = max(seqParams.bFactor);
                    % for the max bvalue applied to one axis, calculate the min duration of the gradient
                    tau = round(obj.sysDiff.maxGrad / obj.sysDiff.maxSlew,4); %ramp duration for diff gradients (shortest possible from specs) 
                    %Delta: distance between the two diffusion gradients
                    Delta = obj.delayTE1 + mr.calcDuration(obj.gz180); %function of TE and pulse duration
                    %delta: duration of diffusion gradient, considering its trapezoidal shape and inverting the b formula (RS)
                    c = bFactor_max*1e6/((2*pi)^2*obj.sysDiff.maxGrad^2);
                    p = [1, -3*Delta, tau^2/2, -tau^3/10 + 3*c];
                    delta_all = roots(p);
                    delta = min(real(delta_all( ...
                        abs(imag(delta_all)) < 1e-10 & ...
                        real(delta_all) > 0 )));
                    delta = obj.roundUpToGRT(delta);
                    bRef=bFactCalc(obj.sysDiff.maxGrad/obj.sys.gamma, tau, delta, Delta, obj);
                    disp(['Reference b-value for delta: ' num2str(bRef) 's/mm2'])                  
                    gDiff_flattime = delta-tau;
                    %is the gradient fitting the delays we have
                    assert(delta+tau<=obj.delayTE2, 'Assertion for delayTE2 failed.'); 
                    assert(delta+tau<=obj.delayTE1, 'Assertion for delayTE1 failed.'); 
                end
                
                if bFactor > 0
                    bDirNorm = obj.bDir(i,:) / norm(obj.bDir(i,:)); % make sure directions are normalized
                else 
                    bDirNorm = [0 0 0];
                end 

                g = obj.sysDiff.maxGrad * sqrt(bFactor/bRef);

                bAxis = zeros(1,3);
                gAxisMHz = zeros(1,3);

                for ja = 1:3 %axis index
                    dir = obj.axesOrder{ja}; % 1:x, 2:y, 3:z
                                                                                      
                    gAxis = g * bDirNorm(ja) * obj.axesSign(ja);
                    assert(abs(gAxis)<=obj.sysDiff.maxGrad, 'Assertion for Diff max Grad failed.'); %gAxis shall be smaller than max allowed G. 
                    
                    obj.gDiff{i,ja}=mr.makeTrapezoid(dir,'amplitude',gAxis,'riseTime',tau,'flatTime',gDiff_flattime,'system',obj.sysDiff);
                    bAxis(ja)=bFactCalc(gAxis/obj.sys.gamma, tau, delta, Delta, obj);
                    gAxisMHz(ja) = gAxis / 1e6;
                end
                fprintf(['b-enc = %d, b = (%.1f, %.1f, %.1f) s/mm2, g = (%.1f, %.1f, %.1f) MHz/m\n'], i, bAxis, gAxisMHz);
            end

            %% Time from trigger to scanner acquisition
            obj.triggerToScannerAcqDelay = obj.gradFreeTime + obj.adc.delay; 
            
            %% Calculate required camera acquisition duration
            obj.cameraAcqDuration = mr.calcDuration(obj.gx) ...
                                  + 1e-3; % To be safe 
            obj.cameraAcqDuration = ceil(obj.cameraAcqDuration*1000)/1000;

            %% Determine chronological order for slice positions

            %  Example for 10 slices
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
                    lin = 1;
                    obj = runKernel(obj, lin, slc, avg, rep, 1, KernelMode.Sync);
                end
                
                %% Add pause and reset flags
                if obj.preScanPause < 4
                    warning('The pause between the synchronization and imaging scans should be equal or larger than 4 seconds. The current value is okay for simulation purposes.');
                end

                obj.addBlock(mr.makeDelay(obj.preScanPause), mr.makeLabel('SET','LIN', 0), mr.makeLabel('SET','SLC', 0), mr.makeLabel('SET','AVG', 0));
    
            end
            
            %% Main sequence body
            
            % Dummy scans: repeated nDummy-times for no b-encoding
            for rep = 1:obj.nDummy %number of dummy volumes (i.e., they loop through the whole slice volume)
                for lin = 1:obj.Ny
                    for slc = 1:obj.nSlices
                        avg = 1;
                        obj = runKernel(obj, lin, slc, avg, rep, 1, KernelMode.Dummy); %bValue = 1 has no b-encoding
                    end
                end
            end

            for bValue=1:seqParams.nbValues               
                % Actual imaging sequence
                for rep=1:obj.nRep
                    for lin = 1:obj.Ny
                        for slc = 1:obj.nSlices
                            avg = 1;
                            obj = runKernel(obj, lin, slc, avg, rep, bValue, KernelMode.Imaging);
                        end
                    end
                end
            end

            %% Set the number of imaging triggers
            obj.nTrig = obj.nSlices*obj.nRep*obj.Ny*seqParams.nbValues;

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
            obj.seq.setDefinition('Name', 'spir2d');
            obj.seq.setDefinition('TE', obj.TE);
            obj.seq.setDefinition('TR', obj.TR);
            obj.seq.setDefinition('TRvolume', obj.TR * obj.nSlices * obj.Ny);

            %% Parameters needed be added to the scanner data header for trajectory merging
            obj.seq.setDefinition('TriggerToScannerAcqDelay', obj.triggerToScannerAcqDelay); 

            %% Parameters to be set on the user interface of the Field Camera            
            % The number of actually acquired dynamics depends on the CameraTrigIgnore.
            obj.seq.setDefinition('CameraNrSyncDynamics', obj.nSyncDynamics); 
            obj.seq.setDefinition('CameraNrDynamics', ceil(obj.nTrig/obj.skipFactor));  
            obj.seq.setDefinition('CameraAcqDuration', obj.cameraAcqDuration); 
            obj.seq.setDefinition('CameraAqDelay', 0);
            obj.seq.setDefinition('CameraTrigIgnore', obj.cameraInterleaveTR); 
            obj.seq.setDefinition('AdcSampleTime', obj.adc.dwell);             
            obj.seq.setDefinition('Matrix', [n n]);
            obj.seq.setDefinition('EncodingMatrix', [obj.adc.numSamples obj.Ny]);
            obj.seq.setDefinition('NrSpiralInterleaves', obj.Ny);
            obj.seq.setDefinition('SpiralSamplesPerInterleave', obj.adc.numSamples);
            obj.seq.setDefinition('InplaneAcceleration', 1);
            obj.seq.setDefinition('SliceShifts', obj.slicePositionChronological); 
            obj.seq.setDefinition('readDir_SCT', readDir_SCT);
            obj.seq.setDefinition('phaseDir_SCT', phaseDir_SCT);
            obj.seq.setDefinition('sliceDir_SCT', sliceDir_SCT);
            obj.seq.setDefinition('SequenceType', 'SE');
            obj.seq.setDefinition('SliceOrdering', 'INTERLEAVED');
            % this is important for making the sequence run automatically
            % on siemens scanners without further parameter tweaking
            obj.seq.setDefinition('MaxAdcSegmentLength', adcSamplesPerSegment); 

            %% Write to pulseq file
            if not(isfolder('exports'))
                mkdir('exports')
            end

            if not(isfolder(strcat('exports/',string(seqParams.scannerType))))
                mkdir(strcat('exports/',string(seqParams.scannerType)))
            end

            filename = strcat('exports/',string(seqParams.scannerType),'/skope_se_spiral_2d_diff','_',string(obj.sliceOrientation),'_',string(obj.phaseEncDir),'_',string(obj.mode));           

            if obj.doPlayFatSat == 1
                filename = strcat(filename, '_fs');
            end

            if isprop(seqParams, 'seqSpecName') && ~isempty(seqParams.seqSpecName)
                filename = strcat(filename, '_', seqParams.seqSpecName);																				
            end

            obj.seq.write(strcat(filename,'.seq')); 
            

            %% Safety checks
            ascfile = fullfile('dependencies','asc', specs.HWfilename);
            if isfile(ascfile)
                % PNS
                fprintf('PNS and CNS computation: using hardware file %s \n', ascfile);
                [ok, pns_norm, pns_comp, t_axis] = obj.seq.calcPNS(ascfile);
                maxPNS = max(pns_norm(1,:));
                maxCNS = max(pns_norm(2,:));
    
                if ok(1)
                    fprintf('PNS check passed successfully (%.1f%% < 100%%)\n', ...
                        100*maxPNS);
                else
                    fprintf('PNS check failed (%.1f %% >= 100%%)\n', ...
                        100*maxPNS);
                end
                
                if ok(2)
                    fprintf('CNS check passed successfully (%.1f%% < 100%%)\n', ...
                        100*maxCNS);
                else
                    fprintf('CNS check failed (%.1f %% >= 100%%)\n', ...
                        100*maxCNS);
                end
    
                % Gradient spectrum check
                [R, Rax, F] = gradSpectrum_latest(obj.seq,ascfile); % waiting for pulseq 1.5.2
            else
                warning('Hardware .asc file not found in dependencies/asc. PNS/CNS and Gradient Spectrum checks were skipped.')
            end

        end

    end

    methods (Access=private)  

        function [Gx_rot, Gy_rot] = rotate_spiralArm(obj, Gx, Gy, phi)
            Gx_rot = cos(phi) * Gx - sin(phi) * Gy;
            Gy_rot = sin(phi) * Gx + cos(phi) * Gy;    
        end

        function obj = runKernel(obj, lin, slc, avg, rep, bValue, mode)

            if not(isa(mode, 'KernelMode'))
                error('Expected a kernel mode argument')
            end

             %% Set ONCE-flag to avoid repeating sync and dummy scans
            if mode == KernelMode.Sync || mode==KernelMode.Dummy
                % ONCE=1 marks the blocks that are only executed in the first repetition
                obj.addBlock(mr.makeLabel('SET','ONCE', 1));
            else
                % Blocks with ONCE=0 are executed on every repetition
                obj.addBlock(mr.makeLabel('SET','ONCE', 0));
            end

            %% RF and ADC settings
            if mode==KernelMode.Dummy || mode==KernelMode.Imaging
                if obj.doPlayFatSat
                    obj.addBlock(obj.gz_fs_pre, obj.gx_fs_pre, obj.gy_fs_pre);
                    obj.addBlock(obj.rf_fs, obj.gz_fs, obj.gx_fs, obj.gy_fs);
                end
                obj.rf.freqOffset = obj.gz.amplitude * obj.slicePositionChronological(slc);                															
                obj.rf.phaseOffset = -2*pi*obj.rf.freqOffset * mr.calcRfCenter(obj.rf); % Compensate for the slice-offset induced phase
                obj.addBlock(obj.rf, obj.gz, mr.makeLabel('SET','PMC',false));
            else
                if obj.doPlayFatSat
                    obj.addBlock(obj.gz_fs_pre, obj.gx_fs_pre, obj.gy_fs_pre);
                    obj.addBlock(obj.gz_fs, obj.gx_fs, obj.gy_fs);
                end
                obj.addBlock(obj.gz, mr.makeLabel('SET','PMC',true));
            end  
            obj.addBlock(obj.gzReph);
                
            %% Diffusion gradients
            % if bValue<2 % b0
            %     obj.addBlock(mr.makeDelay(obj.delayTE1));
            % else % nonzero b-encoding               
            obj.addBlock(obj.gDiff{bValue,1}, obj.gDiff{bValue,2}, obj.gDiff{bValue,3});
            obj.addBlock(mr.makeDelay(obj.delayTE1-mr.calcDuration(obj.gDiff{2,1}))); 
            % end

            if mode==KernelMode.Dummy || mode==KernelMode.Imaging              
                obj.rf180.freqOffset=obj.gz180.amplitude * obj.slicePositionChronological(slc); 
                obj.rf180.phaseOffset=-2*pi*obj.rf180.freqOffset * mr.calcRfCenter(obj.rf180); % compensate for the slice-offset induced phase
                obj.addBlock(obj.rf180, obj.gz180);
            else
                obj.addBlock(obj.gz180);
            end
            
            % if bValue<2 % b0
            %     obj.addBlock(mr.makeDelay(obj.delayTE2));
            % else % nonzero b-encoding 
            obj.addBlock(obj.gDiff{bValue,1}, obj.gDiff{bValue,2}, obj.gDiff{bValue,3});
            obj.addBlock(mr.makeDelay(obj.delayTE2-mr.calcDuration(obj.gDiff{2,1})));                
            % end

            if mode==KernelMode.Sync || mode==KernelMode.Imaging
                obj.addBlock(obj.extTrigger,mr.makeDelay(obj.gradFreeTime)); 												   
            else
                obj.addBlock(mr.makeDelay(obj.gradFreeTime));
            end
          
            %% Set labels
            labels = [  {mr.makeLabel('SET','LIN', lin-1)}, ...
                        {mr.makeLabel('SET','SLC', slc-1)}, ...        
                        {mr.makeLabel('SET','AVG', avg-1)}, ...
                        {mr.makeLabel('SET','REP', rep-1)}, ...
                        {mr.makeLabel('SET','SET', bValue-1)}, ...
                        ];

            %% Rotate and add readout
            phi = (2 * pi * lin / obj.Ny);
            [gx, gy] = obj.rotate_spiralArm(obj.gspiral(1,:), obj.gspiral(2,:), phi);
            obj.gx = mr.makeArbitraryGrad(obj.axesOrder{1},gx, ...
                'system', obj.sys, ...
                'first', 0, ...
                'last', 0);
            obj.gy = mr.makeArbitraryGrad(obj.axesOrder{2},gy, ...
                'system', obj.sys, ...
                'first', 0, ...
                'last', 0);

            if mode==KernelMode.Sync || mode==KernelMode.Imaging
                obj.addBlock(obj.gx, obj.gy, obj.adc, mr.makeLabel('SET','ECO', 0), labels{:});
            else
                obj.addBlock(obj.gx, obj.gy, mr.makeLabel('SET','ECO', 0), labels{:});
            end
        
            %% Refocusing
            [gxRefoc, gyRefoc] = mr.rotate(obj.axesOrder{3},phi,{obj.gxRefoc, obj.gyRefoc});
            obj.addBlock(gxRefoc, gyRefoc);

            %% Add delay
            obj.addBlock(mr.makeDelay(obj.fillTR));
        
        end
    end
end

function b=bFactCalc(G, tau, delta, Delta, obj)
    gamma = obj.sys.gamma*2*pi;
    b = (gamma)^2 * G^2 * ( delta^2 * (Delta-delta/3) + tau^3/30 - tau^2*delta/6 ) *1e-6;
end