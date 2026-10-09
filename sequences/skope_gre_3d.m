classdef skope_gre_3d < PulseqBase
% This is a demo monopolar multi-echo gradient-echo sequence, which includes
% synchronization scans for field-monitoring with a Skope Field Camera and
% uses the LABEL extension. The member method plot() can be used to display
% the generated sequence. 
% 
% Notes:
% - The sequence file is written into the current folder.
% - The TR refers to the excitation repetition time in this example and not
%   the slice TR.
% - The k-space trajectory during the synchronization scans will not be
%   correctly shown by the member method plot().
% - The x-axis is flipped because of a bug in the Siemens Pulseq 
%   interpreter 1.4.0. 
%
% Example:
%  gre = skope_gre_3d(sequenceParams);
%  gre.plot();
%  gre.test();
%
% See also PulseqBase

% (c) 2026 Skope Magnetic Resonance Technologies AG

    properties (Access=private)

        % Pulseq transmit object
        rf

        % Pulseq ADC event
        adc
        
        % Pulseq prewinding gradient
        gxPre

        % Pulseq readout gradient
        gx

        % Pulseq partition selection gradient
        gz

        % Pulseq partition refocusing gradient
        gzReph

        % Pulseq read rewinding gradient
        gxFlyBack

        % Pulseq read spoiling gradient
        gxSpoil

        % Pulseq partition spoiling gradient
        gzSpoil

        % Phase encoding moments
        phaseAreaY

        % Partition encoding moments
        phaseAreaZ

        % RF phase
        rf_phase

        % RF phase increment
        rf_inc

        % Fat suppression gradients
        rf_fs

        gx_fs
        gy_fs
        gz_fs
        gx_fs_pre
        gy_fs_pre
        gz_fs_pre

         % Fat shift [Unit: ppm]
        sat_ppm = -3.45;

        % Play out fat saturation pulse
        doPlayFatSat = false; 

        % Phase increment for RF spoiling
        rfSpoilingInc = 117 

        % Trigger output channel
        triggerOutput;

        % Acceleration factor (Phase)
        accFacPE

    end

    methods

        function obj = skope_gre_3d(seqParams)

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

            %% Create object
            obj.sys = mr.opts(  'MaxGrad', seqParams.maxGrad, ...
                                'GradUnit', 'mT/m', ...
                                'MaxSlew', seqParams.maxSlew, ...
                                'SlewUnit', 'T/m/s', ...
                                'rfRingdownTime', 20e-6, ...
                                'rfDeadTime', 100e-6, ...
                                'adcDeadTime', 10e-6);  

            % Copy all sequence parameters
            fieldNames = fields(seqParams);
            for i = 1:numel(fieldNames)
                fieldname = fieldNames{i};
                if isprop(obj,fieldname)
                    obj.(fieldname) = seqParams.(fieldname);
                end
            end
            
            % Number of phase encoding steps with undersampling
            obj.Ny = round(obj.Nx/obj.accFacPE); 

            % Number of phase encoding steps
            obj.Nz = obj.Nx; 

            %% Axes order
            [obj.axesOrder, obj.axesSign, readDir_SCT, phaseDir_SCT, sliceDir_SCT] ...
                = GetAxesOrderAndSign(obj.sliceOrientation,obj.phaseEncDir);

            % Readout time
            Tpre = obj.readoutTime;

            %% Create a new sequence object
            obj.seq = mr.Sequence(obj.sys);  

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

            %% Create non-selective pulse
            obj.rf = mr.makeBlockPulse(obj.alpha*pi/180,...
                'system',obj.sys,...
                'Duration',0.2e-3,...
                'use','excitation');
            obj.rf_phase = 0;
            obj.rf_inc = 0;

            %% Time for probe excitation
            obj.gradFreeTime = obj.roundUpToGRT(200e-6);
          
            %% Define other gradients and ADC events (Not that X gradient has been flipped here)
            deltak = 1./obj.fov;
            deltak(2) = deltak(2) * obj.accFacPE;
            obj.gx = mr.makeTrapezoid(  obj.axesOrder{1}, ...
                                        'FlatArea', obj.Nx*deltak(1), ...
                                        'FlatTime', obj.readoutTime, ...
                                        'system',obj.sys);
            obj.adc = mr.makeAdc(obj.Nx, ...
                                'Duration', obj.gx.flatTime, ...
                                'Delay', obj.gx.riseTime, ...
                                'system', obj.sys);
            obj.gxPre = mr.makeTrapezoid(obj.axesOrder{1},obj.sys,'Area',-obj.gx.area/2,'Duration',Tpre);
            % Create flyback gradients
            nEchoes = length(obj.TE);
            obj.gxFlyBack = repmat( ...
                mr.makeTrapezoid(obj.axesOrder{1}, 'Area', -obj.gx.area, 'system', obj.sys), ...
                1, nEchoes-1);
            obj.gxSpoil = mr.makeTrapezoid(obj.axesOrder{1},obj.sys,'Area',obj.gx.area,'Duration',Tpre*2);
            obj.phaseAreaY = ([(obj.Ny-1):-1:0]-obj.Ny/2)*deltak(2);
            obj.phaseAreaZ = ([(obj.Nz-1):-1:0]-obj.Nz/2)*deltak(3);

            %% Calculate minimal TEs
            minTE = zeros(size(obj.TE));
            obj.fillTE = zeros(size(obj.TE));
            
            % First echo
            minTE(1) = mr.calcDuration(obj.rf) ...
                     - mr.calcRfCenter(obj.rf) ...
                     - obj.rf.delay ...
                     + obj.gradFreeTime ...
                     + mr.calcDuration(obj.gxPre) ...
                     + mr.calcDuration(obj.gx)/2;
            
            disp(['Minimal TE1 is ' num2str(minTE(1)*1000) ' ms']);
            
            obj.fillTE(1) = obj.roundUpToGRT(obj.TE(1) - minTE(1));
            
            assert(obj.fillTE(1) >= 0, 'Assertion for TE1 failed.');

            %% Subsequent echoes
            for i = 2:nEchoes
                minTE(i) = obj.TE(i-1) ...
                         + mr.calcDuration(obj.gx)/2 ...
                         + mr.calcDuration(obj.gxFlyBack(i-1)) ...
                         + mr.calcDuration(obj.gx)/2;
            
                disp(['Minimal TE' num2str(i) ' is ' num2str(minTE(i)*1000) ' ms']);
            
                obj.fillTE(i) = obj.roundUpToGRT(obj.TE(i) - minTE(i));
                assert(obj.fillTE(i) >= 0, 'Assertion for TE(%d) failed.', i);
            
                % Absorb spacing into preceding flyback gradient
                if obj.fillTE(i) > 0
                    obj.gxFlyBack(i-1) = mr.makeTrapezoid( ...
                        obj.axesOrder{1}, ...
                        'Area', -obj.gx.area, ...
                        'system', obj.sys, ...
                        'Duration', mr.calcDuration(obj.gxFlyBack(i-1)) + obj.fillTE(i));
                    obj.fillTE(i) = 0;
                end
            end

            %% Total flyback duration
            flyBackDuration = 0;
            for i = 1:nEchoes-1
                flyBackDuration = flyBackDuration + mr.calcDuration(obj.gxFlyBack(i));
            end
            
            %% Calculate minimal TR
            minTR = mr.calcDuration(obj.rf) ...
                  + obj.gradFreeTime ...
                  + obj.fillTE(1) ...
                  + mr.calcDuration(obj.gxPre) ...
                  + nEchoes * mr.calcDuration(obj.gx) ...
                  + flyBackDuration ...      
                  + mr.calcDuration(obj.gxSpoil);
            if obj.doPlayFatSat
                minTR = minTR + mr.calcDuration(obj.gz_fs) + mr.calcDuration(obj.gz_fs_pre);
            end
            
            disp(['Minimal TR is ' num2str(minTR*1000) ' ms']);
            
            obj.fillTR = obj.roundUpToGRT(obj.TR - minTR);
            assert(obj.fillTR >= 0, 'Assertion for TR failed.');

            %% Time from trigger to scanner acquisition
            obj.triggerToScannerAcqDelay = obj.gradFreeTime ...
                                           + mr.calcDuration(obj.gxPre) ...
                                           + obj.adc.delay;

            %% Prepare trigger
            obj.extTrigger = mr.makeDigitalOutputPulse(obj.triggerOutput,'duration', obj.sys.gradRasterTime);

                        
            %% Calculate required camera acquisition duration
            obj.cameraAcqDuration = obj.gradFreeTime ...
                                  + mr.calcDuration(obj.gxPre) ...
                                  + nEchoes * mr.calcDuration(obj.gx) ...
                                  + flyBackDuration ...      
                                  + 1e-3; % To be safe 
            
            %% Synchronization
            if obj.nSyncDynamics > 0
                for avg = 1:obj.nSyncDynamics
                    par = 1;
                    lin = 1;
                    obj = runKernel(obj, lin, par, avg, KernelMode.Sync);
                end
            
                %% Add pause and reset flags
                if obj.preScanPause < 4
                    warning('The pause between the synchronization and imaging scans should be equal or larger than 4 seconds. The current value is okay for simulation purposes.');
                end
            
                obj.addBlock(mr.makeDelay(obj.preScanPause), ...
                    mr.makeLabel('SET','LIN', 0), ...
                    mr.makeLabel('SET','PAR', 0), ...
                    mr.makeLabel('SET','AVG', 0), ...
                    mr.makeLabel('SET','ECO', 0));
            end

            %% Drive magnetization to steady state with Ny dummies
            for rep = 1:obj.nDummy
                for i = 1:obj.Ny
                    lin = i;
                    par = 1;
                    obj = runKernel(obj, lin, par, 1, KernelMode.Dummy);
                end
            end
            %% Actual imaging sequence
            % loop over phase encodes and define sequence blocks
            for par = 1:obj.Nz
                for lin = 1:obj.Ny 
                    % loop over partitions
                    avg = 1;
                    obj = runKernel(obj, lin, par, avg, KernelMode.Imaging);   
                end
            end

            %% Set the number of imaging triggers
            obj.nTrig = obj.Ny * obj.Nz;

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
            obj.seq.setDefinition(' Units of time - seconds', '');
            obj.seq.setDefinition(' Units of length - meters', '');
            obj.seq.setDefinition('Name', 'gre3d');
            obj.seq.setDefinition('FOV', obj.fov);
            obj.seq.setDefinition('TR', obj.TR);
            obj.seq.setDefinition('TE', obj.TE);            

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
            obj.seq.setDefinition('Matrix', [obj.Nx obj.Ny*obj.accFacPE obj.Nz]); 
            obj.seq.setDefinition('Encoding', [obj.Nx obj.Ny obj.Nz]);
            obj.seq.setDefinition('InplaneAcceleration', obj.accFacPE);
            obj.seq.setDefinition('readDir_SCT', readDir_SCT);
            obj.seq.setDefinition('phaseDir_SCT', phaseDir_SCT);
            obj.seq.setDefinition('sliceDir_SCT', sliceDir_SCT);  
            obj.seq.setDefinition('SequenceType', 'GRE');

            %% Write to Pulseq file
            if not(isfolder('exports'))
                mkdir('exports')
            end
            if not(isfolder(strcat('exports/',string(seqParams.scannerType))))
                mkdir(strcat('exports/',string(seqParams.scannerType)))
            end
            
            filename = strcat('exports/',string(seqParams.scannerType),'/skope_gre_3d','_',string(obj.sliceOrientation),'_',string(obj.phaseEncDir));  
            
            if obj.doPlayFatSat == 1
                filename = strcat(filename, '_fs');
            end

            if obj.accFacPE > 1  
                filename = strcat(filename, '_R', num2str(obj.accFacPE));  
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
                [R, Rax, F] = gradSpectrum_latest(obj.seq,ascfile);
            else
                warning('Hardware .asc file not found in dependencies/asc. PNS/CNS and Gradient Spectrum checks were skipped.')
            end
            
        end    
    end

    methods (Access=private)               
        function obj = runKernel(obj, lin, par, avg, mode)

            if not(isa(mode, 'KernelMode'))
                error('Expected a kernel mode argument')
            end
        
            %% RF and ADC settings
            if mode==KernelMode.Dummy || mode==KernelMode.Imaging
                if obj.doPlayFatSat
                    obj.addBlock(obj.gz_fs_pre, obj.gx_fs_pre, obj.gy_fs_pre);
                    obj.addBlock(obj.rf_fs, obj.gz_fs, obj.gx_fs, obj.gy_fs);
                end
                obj.rf.phaseOffset = obj.rf_phase * pi/180;
                obj.adc.phaseOffset = obj.rf.phaseOffset;
                obj.addBlock(obj.rf, mr.makeDelay(obj.fillTE(1) + mr.calcDuration(obj.rf)));
                obj.rf_inc = mod(obj.rf_inc + obj.rfSpoilingInc, 360);
                obj.rf_phase = mod(obj.rf_phase + obj.rf_inc, 360);
            else
                if obj.doPlayFatSat
                    obj.addBlock(obj.gz_fs_pre, obj.gx_fs_pre, obj.gy_fs_pre);
                    obj.addBlock(obj.gz_fs, obj.gx_fs, obj.gy_fs);
                end
                obj.rf.phaseOffset = 0;
                obj.adc.phaseOffset = 0;
                obj.addBlock(mr.makeDelay(mr.calcDuration(obj.rf) + obj.fillTE(1)));
            end                      
        
            %% External trigger and gradient-free interval
            if mode==KernelMode.Sync || mode==KernelMode.Imaging
                obj.addBlock(obj.extTrigger,mr.makeDelay(obj.gradFreeTime)); 												   
            else
                obj.addBlock(mr.makeDelay(obj.gradFreeTime));
            end

            %% Read-prewinding and phase encoding gradients
            gyPre = mr.makeTrapezoid(obj.axesOrder{2}, ...
                    'Area', obj.phaseAreaY(lin), ...
                    'Duration', mr.calcDuration(obj.gxPre), ...
                    'system',obj.sys);
            gzPre = mr.makeTrapezoid(obj.axesOrder{3}, ...
                    'Area', obj.phaseAreaZ(par), ...
                    'Duration', mr.calcDuration(obj.gxPre), ...
                    'system',obj.sys);
            obj.addBlock(obj.gxPre,gyPre,gzPre);

            %% All LABELS / counters an flags are automatically initialized to 0 in the beginning, no need to define initial 0's  
            % so we will just increment LIN after the ADC event (e.g. during the spoiler)         
            %seq.addBlock(mr.makeDelay(1)); % older scanners like Trio may need this
            % dummy delay to keep up with timing

            %% Readout gradients
            if mode==KernelMode.Sync || mode==KernelMode.Imaging
                % Set labels
                labels = [  {mr.makeLabel('SET','LIN', lin-1)}, ...
                         {mr.makeLabel('SET','PAR', par-1)}, ...        
                         {mr.makeLabel('SET','AVG', avg-1)}];

                % First echo
                obj.addBlock(obj.gx, obj.adc, mr.makeLabel('SET','ECO', 0), labels{:});
                % Remaining echoes 
                for i = 2:length(obj.TE)
                    obj.addBlock(obj.gxFlyBack(i-1));
                    obj.addBlock(obj.gx, obj.adc, mr.makeLabel('SET','ECO', i-1), labels{:});
                end
            else 
                % First echo
                obj.addBlock(obj.gx);
                % Remaining echoes 
                for i = 2:length(obj.TE)
                    obj.addBlock(obj.gxFlyBack(i-1));
                    obj.addBlock(obj.gx);
                end
            end
          
            %% Negative Phase encoding
            gyPre.amplitude = -gyPre.amplitude;
            gzPre.amplitude = -gzPre.amplitude;
        
            %% Spoiling
            spoilBlockContents = {obj.gxSpoil, gyPre, gzPre};
            obj.addBlock(spoilBlockContents{:});

            %% Add delay
            obj.addBlock(mr.makeDelay(obj.fillTR));
        
        end
    end
end