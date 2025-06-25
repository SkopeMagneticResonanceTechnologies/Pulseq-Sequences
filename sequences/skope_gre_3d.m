classdef skope_gre_3d < PulseqBase
% This is a demo monopolar dual-echo gradient-echo sequence, which includes
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
%
% Example:
%  gre = skope_gre_3d(sequenceParams);
%  gre.plot();
%  gre.test();
%
% See also PulseqBase

% (c) 2024 Skope Magnetic Resonance Technologies AG

    properties (Access=private)

        % Pulseq transmit object
        rf

        % Pulseq ADC event
        adc
        
        % Pulseq prewinding gradient
        gxPre

        % Pulseq readout gradient
        gx

        % Pulseq slice selection gradient
        gz

        % Pulseq slice refocusing gradient
        gzReph

        % Pulseq read rewinding gradient
        gxFlyBack

        % Pulseq read spoiling gradient
        gxSpoil

        % Pulseq slice spoiling gradient
        gzSpoil

        % Phase encoding moments
        phaseAreaY

        % Partition encoding moments
        phaseAreaZ

        % RF phase
        rf_phase

        % RF phase increment
        rf_inc

        % Phase increment for RF spoiling
        rfSpoilingInc = 117 

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

            %% Used gradient amplitude and slew rate by this sequence
            
           
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

            % Field of view [Unit: m]
            obj.fov = seqParams.fov; 
            
            % % Number of readout samples
            obj.Nx = seqParams.Nx; 
            
            % Number of phase encoding steps
            obj.Ny = seqParams.Ny; 

            % Number of phase encoding steps
            obj.Nz = seqParams.Nz; 
            
            % Flip angle [Unit: deg]
            obj.alpha = seqParams.alpha;     
                       
            % Echo times [Unit: s] + one millisecond for phase estimation
            obj.TE = seqParams.TE;
            
            % Excitation repetition time [Unit: s]
            obj.TR = seqParams.TR;                       
            
            % ADC duration [Unit: s]
            obj.readoutTime = seqParams.readoutTime;  

            obj.sliceOrientation = seqParams.sliceOrientation;
            obj.phaseEncDir = seqParams.phaseEncDir;

            % Number of dummy shots for steady state
            obj.nDummy = seqParams.nDummy;

            %% Axes order
            [obj.axesOrder, obj.axesSign, readDir_SCT, phaseDir_SCT, sliceDir_SCT] ...
                = GetAxesOrderAndSign(obj.sliceOrientation,obj.phaseEncDir);
            
            %% Sync scans
            obj.nSyncDynamics = 5;

            Tpre = obj.readoutTime;

            %% Create a new sequence object
            obj.seq = mr.Sequence(obj.sys);  

            %% Create non-selective pulse
            obj.rf = mr.makeBlockPulse(obj.alpha*pi/180,obj.sys,'Duration',0.2e-3);
            obj.rf_phase = 0;
            obj.rf_inc = 0;

            %% Time for probe excitation
            obj.gradFreeTime = obj.roundUpToGRT(200e-6);
          
            %% Define other gradients and ADC events (Not that X gradient has been flipped here)
            deltak = 1./obj.fov;
            obj.gx = mr.makeTrapezoid(  obj.axesOrder{1}, ...
                                        'FlatArea', obj.Nx*deltak(1), ...
                                        'FlatTime', obj.roundUpToGRT(obj.readoutTime), ...
                                        'system',obj.sys);

            % Center ADC on gradient
            adcShift = (obj.roundUpToGRT(obj.readoutTime)-obj.readoutTime)/2;
            adcShift = ceil(adcShift/1e-6)*1e-6;

            obj.adc = mr.makeAdc(obj.Nx, ...
                                'Duration', obj.readoutTime, ...
                                'Delay', obj.gx.riseTime + adcShift, ...
                                'system', obj.sys);
            obj.gxPre = mr.makeTrapezoid(obj.axesOrder{1},obj.sys,'Area',-obj.gx.area/2,'Duration',obj.roundUpToGRT(Tpre));
            obj.gxFlyBack = mr.makeTrapezoid(obj.axesOrder{1},'Area',-obj.gx.area,'system',obj.sys);  
            obj.gxSpoil = mr.makeTrapezoid(obj.axesOrder{1},obj.sys,'Area',obj.gx.area,'Duration',obj.roundUpToGRT(Tpre)*2);
            obj.phaseAreaY = ([(obj.Ny-1):-1:0]-obj.Ny/2)*deltak(2);
            obj.phaseAreaZ = ([(obj.Nz-1):-1:0]-obj.Nz/2)*deltak(3);

            %% Calculate minimal TEs
            % First echo
            minTE1 = mr.calcDuration(obj.rf) - mr.calcRfCenter(obj.rf) - obj.rf.delay ...
                     + obj.gradFreeTime ...
                     + mr.calcDuration(obj.gxPre) ...
                     + mr.calcDuration(obj.gx)/2;

            disp(['Minimal TE 1 is ' num2str(minTE1*1000) ' ms'])

            obj.fillTE(1) = obj.roundUpToGRT(obj.TE(1) - minTE1);
            assert(obj.fillTE(1) >= 0, 'Assertion for TE1 failed');

            % Calculate minTE for all subsequent echoes
            previousMinTE = minTE1;
            for i=2:numel(obj.TE)
                    minTEx = previousMinTE + mr.calcDuration(obj.gx)/2 ...
                    + mr.calcDuration(obj.gxFlyBack) ...
                    + mr.calcDuration(obj.gx)/2 ...
                    + obj.fillTE(i-1);
                disp(['Minimal TE ' num2str(i) ' is ' num2str(minTEx*1000) ' ms'])

                obj.fillTE(i) = obj.roundUpToGRT(obj.TE(i) - minTEx);
                assert(obj.fillTE(i) >= 0, ['Assertion for TE' num2str(i) ' failed.']);

                previousMinTE = minTEx;
            end
            
            %% Increase duration of fly-back gradient by TE2 fill time
            minTE = min(obj.fillTE(2:end));
            if minTE > 0 
                obj.gxFlyBack = mr.makeTrapezoid(obj.axesOrder{1},'Area',-obj.gx.area,'system',obj.sys,'Duration',mr.calcDuration(obj.gxFlyBack) + minTE); 
                obj.fillTE(2:end) = obj.fillTE(2:end)-minTE;
            end

            %% Calculate minimal TR
            minTR = mr.calcDuration(obj.rf) ...
                  + mr.calcDuration(obj.gradFreeTime) ...
                  + obj.fillTE(1) ...
                  + mr.calcDuration(obj.gxPre) ...
                  + mr.calcDuration(obj.gx);
            
            for i=2:numel(obj.TE)
                minTR = minTR ...
                  + mr.calcDuration(obj.gxFlyBack) ...
                  + obj.fillTE(i) ...
                  + mr.calcDuration(obj.gx);
            end
            minTR = minTR + mr.calcDuration(obj.gxSpoil);
            disp(['Minimal TR is ' num2str(minTR*1000) ' ms'])
            
            obj.fillTR = obj.roundUpToGRT(obj.TR - minTR);
            assert(obj.fillTR >= 0, 'Assertion for TR failed.');

            %% Time from trigger to scanner acquisition
            obj.triggerToScannerAcqDelay = obj.fillTE(1) + mr.calcDuration(obj.gradFreeTime) ...
                                           + mr.calcDuration(obj.gxPre) ...
                                           + obj.adc.delay;

            %% Prepare trigger
            obj.extTrigger = mr.makeDigitalOutputPulse('ext1','duration', obj.sys.gradRasterTime);
                                        
            %% Calculate required camera acquisition duration
            obj.cameraAcqDuration = obj.fillTE(1) + mr.calcDuration(obj.gradFreeTime) ...
                                  + mr.calcDuration(obj.gxPre) ...
                                  + mr.calcDuration(obj.gx) ...
                                  + 1e-3; % To be safe ;
            
            for i=2:numel(obj.TE)
                obj.cameraAcqDuration = obj.cameraAcqDuration  ...
                                      + mr.calcDuration(obj.gxFlyBack) ...
                                      + obj.fillTE(i) ...
                                      + mr.calcDuration(obj.gx);
            end           
            
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

                obj.addBlock(mr.makeDelay(obj.preScanPause), mr.makeLabel('SET','LIN', 0), mr.makeLabel('SET','SLC', 0), mr.makeLabel('SET','AVG', 0));
            end

            %% Drive magnetization to steady state
            for i=1:obj.nDummy
                runKernel(obj, floor(obj.Ny/2), floor(obj.Nz/2), 1, KernelMode.Dummy);
            end

            %% Actual imaging sequence
            % loop over phase encodes and define sequence blocks
            f = waitbar(0,'Creating partitions for 3D sequence. Please wait...');
            for par = 1:obj.Nz
                if isvalid(f)
                    waitbar(par/obj.Nz,f);
                end
                for lin = 1:obj.Ny      
                    % loop over slices
                    avg = 1;
                    obj = runKernel(obj, lin, par, avg, KernelMode.Imaging);   
                end
            end
            if isvalid(f)
                close(f)
            end

            % Set number of expected external triggers
            obj.nTrig = obj.nDummy + obj.Ny * obj.Nz;
            
            %% check whether the timing of the sequence is correct
            fprintf('Checking timing...\n');
            [ok, error_report] = obj.seq.checkTiming;
            
            if (ok)
                fprintf('Timing check passed successfully\n');
            else
                fprintf('Timing check failed! Error listing follows:\n');
                fprintf([error_report{:}]);
                fprintf('\n');
            end
            
            %% Calculate Camera Interleave TR (blank time)
            obj.CalculateInterleaveTR(obj.TR);	        

            %% Prepare sequence export
            obj.seq.setDefinition('Name', 'gre3d');
            obj.seq.setDefinition('FOV', obj.fov);
            obj.seq.setDefinition('TR', obj.TR);
            obj.seq.setDefinition('TE', obj.TE);            

            %% Parameters needed be added to the scanner data header for trajectory merging
            obj.seq.setDefinition('TriggerToScannerAcqDelay', obj.triggerToScannerAcqDelay);   
           
            %% Parameters to be set on the user interface of the Field Camera            
            % The number of actually acquired dynamics depends on the cameraInterleaveTR.
            obj.seq.setDefinition('CameraNrDynamics', obj.nTrig);  
            obj.seq.setDefinition('CameraNrSyncDynamics', obj.nSyncDynamics); 
            obj.seq.setDefinition('CameraAcqDuration', obj.cameraAcqDuration);  
            obj.seq.setDefinition('CameraInterleaveTR', obj.cameraInterleaveTR); 
            obj.seq.setDefinition('CameraAqDelay', 0); 
            obj.seq.setDefinition('AdcSampleTime', obj.adc.dwell); 
            obj.seq.setDefinition('Matrix', [obj.Nx obj.Ny obj.Nz]); 
            obj.seq.setDefinition('readDir_SCT', readDir_SCT);
            obj.seq.setDefinition('phaseDir_SCT', phaseDir_SCT);
            obj.seq.setDefinition('sliceDir_SCT', sliceDir_SCT);
            
            %% Write to Pulseq file
            if not(isfolder('exports'))
                mkdir('exports')
            end 
            obj.seq.write(strcat('exports/skope_gre_3d','_',string(obj.sliceOrientation),'_',string(obj.phaseEncDir),'.seq'));  
            
            
        end    
    end

    methods (Access=private)               
        function obj = runKernel(obj, lin, par, avg, mode)

            if not(isa(mode, 'KernelMode'))
                error('Expected a kernel mode argument')
            end
        
            %% RF and ADC settings
            if mode==KernelMode.Dummy || mode==KernelMode.Imaging
                obj.rf.phaseOffset = mod(117*(lin^2+lin+2)*pi/180,2*pi);
                obj.adc.phaseOffset = obj.rf.phaseOffset;
                obj.addBlock(obj.rf, mr.makeLabel('SET','PMC',false), mr.makeLabel('SET','AVG',avg-1));
            else
                obj.rf.phaseOffset = 0;
                obj.adc.phaseOffset = 0;
                obj.addBlock(mr.makeDelay(mr.calcDuration(obj.rf)),mr.makeLabel('SET','PMC',true), mr.makeLabel('SET','AVG',avg-1));
            end       
        
            %% External trigger and gradient-free interval
            % We send the trigger here always for the dummies to get a
            % steady state field probe signal
            obj.addBlock(obj.extTrigger, mr.makeDelay(obj.gradFreeTime + obj.fillTE(1)));

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

            %% Set labels
            labels = [  {mr.makeLabel('SET','LIN', lin-1)}, ...
                         {mr.makeLabel('SET','PAR', par-1)}, ...        
                         {mr.makeLabel('SET','AVG', avg-1)}];

            %% First readout gradient
            if mode==KernelMode.Sync || mode==KernelMode.Imaging
                obj.addBlock(obj.gx, obj.adc, mr.makeLabel('SET','ECO', 0), labels{:});
            else
                obj.addBlock(obj.gx, mr.makeLabel('SET','ECO', 0), labels{:});
            end
                       
            for i=2:numel(obj.TE)
                %% Fly back
                obj.addBlock(obj.gxFlyBack); % Fill time has been absorbed in gradient duration

                %% Fill time
                if obj.fillTE(i) > 0
                    obj.addBlock(mr.makeDelay(obj.fillTE(i)));
                end
            
                %% Following readout gradient   
                if mode==KernelMode.Sync || mode==KernelMode.Imaging
                    obj.addBlock(obj.gx, obj.adc, mr.makeLabel('SET','ECO', i-1), labels{:});
                else
                    obj.addBlock(obj.gx, mr.makeLabel('SET','ECO', i-1), labels{:});
                end

            end

            %% Negative Phase encoding
            gyPre.amplitude = -gyPre.amplitude;
            gzPre.amplitude = -gzPre.amplitude;
        
            %% Spoiling
            spoilBlockContents = {obj.gxSpoil, gyPre, gzPre};
            obj.addBlock(spoilBlockContents{:});

            %% Add delay
            if mode==KernelMode.Sync 
                % Add one second between each sync scan
                obj.addBlock(mr.makeDelay(obj.fillTR+1));
            else
                obj.addBlock(mr.makeDelay(obj.fillTR));
            end
        end
    end
end