classdef skope_gre_2d < PulseqBase
% This is a demo mono-polar dual-echo gradient-echo sequence, which includes
% synchronization scans for field-monitoring with a Skope Field Camera and
% uses the LABEL extension. The member method plot() can be used to display
% the generated sequence. 
% 
% Notes:
% - The TR refers to the excitation repetition time in this example and not
%   the slice TR.
% - The k-space trajectory during the synchronization scans will not be
%   correctly shown by the member method plot().
%
% Example:
%  gre = skope_gre_2d(sequenceParams);
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
        phaseAreas

        % RF phase
        rf_phase

        % RF phase increment
        rf_inc

        % Phase increment for RF spoiling
        rfSpoilingInc = 117 

        % Play out trigger earlier to compensate for field measurement latency
        triggerLatency = 150e-6;

        doMonitoringDuringRF;

        triggerOutput

    end

    methods

        function obj = skope_gre_2d(seqParams)

            %% Check input structure
            if not(isa(seqParams,'SequenceParams'))
                error('Input needs to be a SequenceParams object.');
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

            % Field of view [Unit: m]
            obj.fov = seqParams.fov; 

            % Number of readout samples
            obj.Nx = seqParams.Nx; 

            % Number of phase encoding steps
            obj.Ny = seqParams.Ny; 
                       
            % Flip angle [Unit: deg]
            obj.alpha = seqParams.alpha;     
            
            % Slice thickness [Unit: m]
            obj.thickness = seqParams.thickness; 

            % Number of slices
            obj.nSlices = seqParams.nSlices;
            
            % Echo times [Unit: s]
            obj.TE = seqParams.TE;

            % Excitation repetition time [Unit: s]
            obj.TR = seqParams.TR;                       
            
            % ADC duration [Unit: s]
            obj.readoutTime = seqParams.readoutTime;  

            % Slice orientation and phase-encoding direction
            obj.sliceOrientation = seqParams.sliceOrientation;
            obj.phaseEncDir = seqParams.phaseEncDir;

            % Number of dummies to reach steady-state
            obj.nDummy = seqParams.nDummy;

            % Perform monitoring during RF
            obj.doMonitoringDuringRF = seqParams.doMonitoringDuringRF;

            % Set trigger output channel
            obj.triggerOutput = seqParams.triggerOutput;

            % Slice distance factor in percent
            obj.distanceFactorPercentage = seqParams.distanceFactorPercentage;

            %% Axes order
            [obj.axesOrder, obj.axesSign, readDir_SCT, phaseDir_SCT, sliceDir_SCT] ...
                = GetAxesOrderAndSign(obj.sliceOrientation,obj.phaseEncDir);

            %% Create a new sequence object
            obj.seq = mr.Sequence(obj.sys);  
            
            %% Time for probe excitation
            if obj.doMonitoringDuringRF
                obj.gradFreeTime = 0;
            else
                obj.gradFreeTime = obj.roundUpToGRT(200e-6);
            end                      

            %% Create alpha-degree slice selection pulse and gradient
            [obj.rf, obj.gz] = mr.makeSincPulse( obj.alpha*pi/180, ...
                                                'Duration', 2e-3, ...
                                                'SliceThickness', obj.thickness, ...
                                                'apodization', 0.5, ...
                                                'timeBwProduct', 4, ...
                                                'system',obj.sys);
            % Set correct axis
            obj.gz.channel =  obj.axesOrder{3};   

            %% Define other gradients and ADC events (Not that X gradient has been flipped here)
            deltak = 1/obj.fov;
            obj.gx = mr.makeTrapezoid(  obj.axesOrder{1}, ...
                                        'FlatArea', obj.Nx*deltak, ...
                                        'FlatTime', obj.readoutTime, ...
                                        'system',obj.sys);
            obj.adc = mr.makeAdc(obj.Nx*obj.ro_os, ...
                                'Duration', obj.gx.flatTime, ...
                                'Delay', obj.gx.riseTime, ...
                                'system', obj.sys);
            obj.gxPre = mr.makeTrapezoid(obj.axesOrder{1},'Area',-obj.gx.area/2,'Duration',1e-3,'system',obj.sys);
            obj.gzReph = mr.makeTrapezoid(obj.axesOrder{3},'Area',-obj.gz.area/2,'Duration',2e-3,'system',obj.sys);
            obj.gxFlyBack = mr.makeTrapezoid(obj.axesOrder{1},'Area',-obj.gx.area,'system',obj.sys);
            obj.phaseAreas = -((0:obj.Ny-1)-obj.Ny/2)*deltak; % phase area should be Kmax for clin=0 and -Kmax for clin=Ny... strange
            
            % gradient spoiling
            obj.gxSpoil = mr.makeTrapezoid(obj.axesOrder{1},'Area', 2*obj.Nx*deltak,'system', obj.sys);
            obj.gzSpoil = mr.makeTrapezoid(obj.axesOrder{3},'Area',4/obj.thickness,'system', obj.sys);
            
            %% Calculate minimal TEs
            % First echo
            minTE1 = obj.gz.flatTime/2 + obj.gz.fallTime +  mr.calcDuration(obj.gzReph) ...
                    + mr.calcDuration(obj.gxPre)  ...
                    + mr.calcDuration(obj.gx)/2;
            disp(['Minimal TE1 is ' num2str((minTE1 + obj.gradFreeTime)*1000) ' ms'])
            
            obj.fillTE(1) = obj.roundUpToGRT(obj.TE(1) - minTE1);
            assert(obj.fillTE(1) >= obj.gradFreeTime, 'Assertion for TE1 failed');
            
            % Second echo
            minTE2 = minTE1 + mr.calcDuration(obj.gx)/2 ...
                    + mr.calcDuration(obj.gxFlyBack) ...
                    + mr.calcDuration(obj.gx)/2 ...
                    + obj.fillTE(1);
            disp(['Minimal TE2 is ' num2str(minTE2*1000) ' ms'])
            
            obj.fillTE(2) = obj.roundUpToGRT(obj.TE(2) - minTE2);
            assert(obj.fillTE(2) >= 0, 'Assertion for TE2 failed.');
            
            %% Time from trigger to scanner acquisition
            obj.triggerToScannerAcqDelay =  obj.fillTE(1) ...
                                           + mr.calcDuration(obj.gxPre) ...
                                           + obj.adc.delay;
            if obj.doMonitoringDuringRF
                obj.triggerToScannerAcqDelay = obj.triggerToScannerAcqDelay + ...
                                                obj.triggerLatency + ...
                                                obj.gz.flatTime/2 + ...
                                                obj.gz.fallTime + ...
                                                mr.calcDuration(obj.gzReph);
            end

            %% Absorb delayTE2 in gradient
            obj.gxFlyBack = mr.makeTrapezoid(obj.axesOrder{1},'Area',-obj.gx.area, ...
                        'Duration', mr.calcDuration(obj.gxFlyBack) + obj.fillTE(2), ...
                        'system', obj.sys);
            obj.fillTE(2) = 0;
            
            %% Prepare trigger
            if obj.doMonitoringDuringRF
                % Play out trigger event in the middle of the slice
                % selection gradient
                triggerDelay = obj.gz.flatTime/2 + obj.gz.riseTime - obj.triggerLatency; 
                assert(triggerDelay>=0,'Trigger delay must be larger than zero.');
                obj.extTrigger = mr.makeDigitalOutputPulse( obj.triggerOutput, ...
                                                            'duration', obj.sys.gradRasterTime, ...
                                                            'delay',triggerDelay);
            else
                obj.extTrigger = mr.makeDigitalOutputPulse(obj.triggerOutput,'duration', obj.sys.gradRasterTime);
            end

            %% Calculate minimal TR
            minTR = mr.calcDuration(obj.gz) ...
                  + mr.calcDuration(obj.gzReph) ...
                  + obj.fillTE(1) ...
                  + mr.calcDuration(obj.gxPre) ...
                  + mr.calcDuration(obj.gx) ...
                  + mr.calcDuration(obj.gxFlyBack) ...
                  + mr.calcDuration(obj.gx) ...
                  + mr.calcDuration(obj.gxSpoil, obj.gzSpoil);
            disp(['Minimal TR is ' num2str(minTR*1000) ' ms'])
            
            obj.fillTR = obj.roundUpToGRT(obj.TR - minTR);
            assert(obj.fillTR >= 0, 'Assertion for TR failed.');
            
            %% Calculate required camera acquisition duration
            obj.cameraAcqDuration = obj.fillTE(1) ...
                                  + mr.calcDuration(obj.gxPre) ...
                                  + mr.calcDuration(obj.gx) ...
                                  + mr.calcDuration(obj.gxFlyBack) ...
                                  + mr.calcDuration(obj.gx) ...
                                  + 1e-3; % To be safe 

            if obj.doMonitoringDuringRF
                obj.cameraAcqDuration = obj.cameraAcqDuration + ...
                                        obj.triggerLatency + ...
                                        obj.gz.flatTime/2 + ...
                                        obj.gz.fallTime + ...
                                        mr.calcDuration(obj.gzReph);
            end

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

            %% Phase settings
            obj.rf_phase = 0;
            obj.rf_inc = 0;
            
            %% Synchronization
            if obj.nSyncDynamics > 0
                for avg = 1:obj.nSyncDynamics
                    slc = 1;
                    lin = 1;
                    obj = runKernel(obj, lin, slc, avg, KernelMode.Sync);
                end
                
                %% Add pause and reset flags
                if obj.preScanPause < 4
                    warning('The pause between the synchronization and imaging scans should be equal or larger than 4 seconds. The current value is okay for simulation purposes.');
                end
    
                obj.addBlock(mr.makeDelay(obj.preScanPause), mr.makeLabel('SET','LIN', 0), mr.makeLabel('SET','SLC', 0), mr.makeLabel('SET','AVG', 0));
            end

            %% Dummies
            for lin = mod(1:obj.nDummy,obj.Ny)+1      
                % loop over slices
                for slc = 1:obj.nSlices
                    avg = 1;
                    obj = runKernel(obj, lin, slc, avg, KernelMode.Dummy);   
                end
            end

            %% Actual imaging sequence
            % loop over phase encodes and define sequence blocks
            for lin = 1:obj.Ny      
                % loop over slices
                for slc = 1:obj.nSlices
                    avg = 1;
                    obj = runKernel(obj, lin, slc, avg, KernelMode.Imaging);   
                end
            end

            % Set number of expected external triggers
            obj.nTrig = obj.Ny * obj.nSlices;
            
            %% check whether the timing of the sequence is correct
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
            obj.seq.setDefinition('Name', 'gre2d');
            obj.seq.setDefinition('FOV', [obj.fov obj.fov obj.thickness*obj.nSlices*(1+obj.distanceFactorPercentage/100)]);
            obj.seq.setDefinition('TR', obj.TR);
            obj.seq.setDefinition('TE', obj.TE);            

            %% Parameters needed be added to the scanner data header for trajectory merging
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
            obj.seq.setDefinition('EncodingMatrix', [obj.adc.numSamples obj.Ny]);
            obj.seq.setDefinition('InplaneAcceleration', 1);
            obj.seq.setDefinition('SliceShifts', obj.slicePositionChronological); 
            obj.seq.setDefinition('readDir_SCT', readDir_SCT);
            obj.seq.setDefinition('phaseDir_SCT', phaseDir_SCT);
            obj.seq.setDefinition('sliceDir_SCT', sliceDir_SCT);
            obj.seq.setDefinition('MonitoringDuringRF', obj.doMonitoringDuringRF);
            obj.seq.setDefinition('SequenceType', 'GRE');
            obj.seq.setDefinition('SliceOrdering', 'INTERLEAVED');
                                     
            %% Write to Pulseq file
            if not(isfolder('exports'))
                mkdir('exports')
            end
            if not(isfolder(strcat('exports/',string(seqParams.scannerType))))
                mkdir(strcat('exports/',string(seqParams.scannerType)))
            end
            
            filename = strcat('exports/',string(seqParams.scannerType),'/skope_gre_2d','_',string(obj.sliceOrientation),'_',string(obj.phaseEncDir));  

            if isprop(seqParams, 'seqSpecName') && ~isempty(seqParams.seqSpecName)
                filename = strcat(filename, '_', seqParams.seqSpecName);																				
            end

            obj.seq.write(strcat(filename,'.seq')); 
            disp(['Storing sequence file "', char(filename), '.seq"'])
            
        end    
    end

    methods (Access=private)               
        function obj = runKernel(obj, lin, slc, avg, mode)

            if not(isa(mode, 'KernelMode'))
                error('Expected a kernel mode argument')
            end
        
           %% RF and ADC settings
            if mode==KernelMode.Dummy
                obj.rf.freqOffset = obj.gz.amplitude  * obj.slicePositionChronological(slc);
                obj.rf.phaseOffset = obj.rf_phase/180*pi - 2*pi*obj.rf.freqOffset * mr.calcRfCenter(obj.rf);  % compensate for slice-offset induced phase
                obj.adc.phaseOffset = obj.rf_phase/180*pi;
                obj.addBlock(obj.rf, obj.gz, mr.makeLabel('SET','PMC',false), mr.makeLabel('SET','AVG',avg-1));
            elseif mode==KernelMode.Imaging
                obj.rf.freqOffset = obj.gz.amplitude * obj.slicePositionChronological(slc);
                obj.rf.phaseOffset = obj.rf_phase/180*pi - 2*pi*obj.rf.freqOffset * mr.calcRfCenter(obj.rf);  % compensate for slice-offset induced phase
                obj.adc.phaseOffset = obj.rf_phase/180*pi;
                if obj.doMonitoringDuringRF
                    obj.addBlock(obj.rf, obj.extTrigger, obj.gz, mr.makeLabel('SET','PMC',false), mr.makeLabel('SET','AVG',avg-1));
                else
                    obj.addBlock(obj.rf, obj.gz, mr.makeLabel('SET','PMC',false), mr.makeLabel('SET','AVG',avg-1));
                end
            elseif mode==KernelMode.Sync
                if obj.doMonitoringDuringRF
                    obj.addBlock(obj.gz, obj.extTrigger, mr.makeLabel('SET','PMC',true), mr.makeLabel('SET','AVG',avg-1));
                else
                    obj.addBlock(obj.gz, mr.makeLabel('SET','PMC',true), mr.makeLabel('SET','AVG',avg-1));
                end
            else
                error('Unknown kernel mode');
            end

            % Update RF spoiling
            obj.rf_inc = mod(obj.rf_inc + obj.rfSpoilingInc, 360.0);
            obj.rf_phase = mod(obj.rf_phase + obj.rf_inc, 360.0);
            
            %% Slice refocusing gradient
            obj.addBlock(obj.gzReph);
        
             %% External trigger and gradient-free interval
            if mode==KernelMode.Sync || mode==KernelMode.Imaging
                if obj.doMonitoringDuringRF
                    obj.addBlock(mr.makeDelay(obj.fillTE(1)));
                else
                    obj.addBlock(obj.extTrigger, mr.makeDelay(obj.fillTE(1)));
                end
            else
                obj.addBlock(mr.makeDelay(obj.fillTE(1)));
            end          

            %% Read-prewinding and phase encoding gradients
            gyPre = mr.makeTrapezoid(obj.axesOrder{2}, ...
                    'Area', obj.phaseAreas(lin), ...
                    'Duration', mr.calcDuration(obj.gxPre), ...
                    'system',obj.sys);
            obj.addBlock(obj.gxPre,gyPre);

            %% All LABELS / counters an flags are automatically initialized to 0 in the beginning, no need to define initial 0's  
            % so we will just increment LIN after the ADC event (e.g. during the spoiler)         
            %seq.addBlock(mr.makeDelay(1)); % older scanners like Trio may need this
            % dummy delay to keep up with timing

            %% Set labels
            labels = [  {mr.makeLabel('SET','LIN', lin-1)}, ...
                        {mr.makeLabel('SET','SLC', slc-1)}, ...        
                        {mr.makeLabel('SET','AVG', avg-1)}];

            %% First readout gradient
            if mode==KernelMode.Sync || mode==KernelMode.Imaging
                obj.addBlock(obj.gx, obj.adc, mr.makeLabel('SET','ECO', 0), labels{:});
            else
                obj.addBlock(obj.gx, mr.makeLabel('SET','ECO', 0), labels{:});
            end
        
            %% Fly back
            obj.addBlock(obj.gxFlyBack);
        
            %% Second readout gradient   
            if mode==KernelMode.Sync || mode==KernelMode.Imaging
                obj.addBlock(obj.gx, obj.adc, mr.makeLabel('SET','ECO', 1));
            else
                obj.addBlock(obj.gx, mr.makeLabel('SET','ECO', 1));
            end
        
            %% Negative Phase encoding
            gyPre.amplitude = -gyPre.amplitude;
        
            %% Spoiling
            spoilBlockContents = {obj.gxSpoil, gyPre, obj.gzSpoil};
            obj.addBlock(spoilBlockContents{:});

            %% Add delay
            obj.addBlock(mr.makeDelay(obj.fillTR));
        
        end
    end
end