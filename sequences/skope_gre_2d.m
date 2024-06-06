classdef skope_gre_2d < PulseqBase
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
% - The x-axis is flipped because of a bug in the Siemens Pulseq 
%   interpreter 1.4.0. 
%
% Example:
%  gre = skope_gre_2d(scannerType);
%  gre.plot();
%  gre.test();
%
% See also PulseqBase

% (c) 2022 Skope Magnetic Resonance Technologies AG

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

    end

    methods

        function obj = skope_gre_2d(scannerType)

            %% Get system limits
            specs = GetMRSystemSpecs(scannerType); 

            if not(strcmpi(specs.maxGrad_unit,'mT/m'))
                error('Expected mT/m for maximum gradient.');
            end

            if not(strcmpi(specs.maxSlew_unit,'T/m/s'))
                error('Expected T/m/s for slew rate.');
            end

            %% Used gradient amplitude and slew rate by this sequence
            maxGrad = 28;
            maxSlew = 150;
           
            %% Check specs
            if maxGrad > specs.maxGrad
                error('Scanner does not support requested gradient amplitude.');
            end
            if maxSlew > specs.maxSlew
                error('Scanner does not support requested slew rate.');
            end

            %% Create object
            obj.sys = mr.opts(  'MaxGrad', maxGrad, ...
                                'GradUnit', 'mT/m', ...
                                'MaxSlew', maxSlew, ...
                                'SlewUnit', 'T/m/s', ...
                                'rfRingdownTime', 20e-6, ...
                                'rfDeadTime', 100e-6, ...
                                'adcDeadTime', 10e-6);  

            % Field of view [Unit: m]
            obj.fov = 200e-3; 

            % Number of readout samples
            obj.Nx = 128; 

            % Number of phase encoding steps
            obj.Ny = obj.Nx; 
                       
            % Flip angle [Unit: deg]
            obj.alpha = 7;     
            
            % Slice thickness [Unit: m]
            obj.thickness = 3e-3; 

            % Number of slices
            obj.nSlices = 5;
            
            % Echo times [Unit: s]
            obj.TE = [6 12] * 1e-3;

            % Excitation repetition time [Unit: s]
            obj.TR = 25e-3;                       
            
            % ADC duration [Unit: s]
            obj.readoutTime = 3.2e-3;             

            %% Create a new sequence object
            obj.seq = mr.Sequence(obj.sys);  
            
            %% Time for probe excitation
            obj.gradFreeTime = obj.roundUpToGRT(200e-6);

            %% Create alpha-degree slice selection pulse and gradient
            [obj.rf, obj.gz] = mr.makeSincPulse( obj.alpha*pi/180, ...
                                                'Duration', 1.5e-3, ...
                                                'SliceThickness', obj.thickness, ...
                                                'apodization', 0.5, ...
                                                'timeBwProduct', 4, ...
                                                'system',obj.sys);
            
            %% Define other gradients and ADC events (Not that X gradient has been flipped here)
            deltak = 1/obj.fov;
            obj.gx = mr.makeTrapezoid(  'x', ...
                                        'FlatArea', obj.signFlip*obj.Nx*deltak, ...
                                        'FlatTime', obj.readoutTime, ...
                                        'system',obj.sys);
            obj.adc = mr.makeAdc(obj.Nx, ...
                                'Duration', obj.gx.flatTime, ...
                                'Delay', obj.gx.riseTime, ...
                                'system', obj.sys);
            obj.gxPre = mr.makeTrapezoid('x','Area',-obj.gx.area/2,'Duration',1e-3,'system',obj.sys);
            obj.gzReph = mr.makeTrapezoid('z','Area',-obj.gz.area/2,'Duration',1e-3,'system',obj.sys);
            obj.gxFlyBack = mr.makeTrapezoid('x','Area',-obj.gx.area,'system',obj.sys);
            obj.phaseAreas = -((0:obj.Ny-1)-obj.Ny/2)*deltak; % phase area should be Kmax for clin=0 and -Kmax for clin=Ny... strange
            
            % gradient spoiling
            obj.gxSpoil = mr.makeTrapezoid('x','Area', obj.signFlip*2*obj.Nx*deltak,'system', obj.sys);
            obj.gzSpoil = mr.makeTrapezoid('z','Area',4/obj.thickness,'system', obj.sys);
            
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

            %% Absorb delayTE2 in gradient
            obj.gxFlyBack = mr.makeTrapezoid('x','Area',-obj.gx.area, ...
                        'Duration', mr.calcDuration(obj.gxFlyBack) + obj.fillTE(2), ...
                        'system', obj.sys);
            obj.fillTE(2) = 0;
            
            %% Prepare trigger
            obj.extTrigger = mr.makeDigitalOutputPulse('ext1','duration', obj.sys.gradRasterTime);

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

            %% Phase settings
            obj.rf_phase = 0;
            obj.rf_inc = 0;
            
            %% Synchronization
            if obj.nSyncDynamics > 0
                for avg = 1:obj.nSyncDynamics
                    slc = 1;
                    lin = 1;
                    obj = runKernel(obj, lin, slc, avg, false);
                end
                
                %% Add pause and reset flags
                if obj.preScanPause < 4
                    warning('The pause between the synchronization and imaging scans should be equal or larger than 4 seconds. The current value is okay for simulation purposes.');
                end
    
                obj.seq.addBlock({mr.makeDelay(obj.preScanPause), mr.makeLabel('SET','LIN', 0), mr.makeLabel('SET','SLC', 0), mr.makeLabel('SET','AVG', 0)});
            end

            %% Actual imaging sequence
            % loop over phase encodes and define sequence blocks
            for lin = 1:obj.Ny      
                % loop over slices
                for slc = 1:obj.nSlices
                    avg = 1;
                    obj = runKernel(obj, lin, slc, avg);   
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
            obj.seq.setDefinition('FOV', [obj.fov obj.fov obj.thickness*obj.nSlices]);
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
            obj.seq.setDefinition('Matrix', [obj.Nx obj.Ny]); 
            obj.seq.setDefinition('SliceShifts', [obj.thickness*([1:obj.nSlices]-1-(obj.nSlices-1)/2)]); 
                        
            %% Write to Pulseq file
            if not(isfolder('exports'))
                mkdir('exports')
            end
            obj.seq.write('exports/skope_gre_2d.seq')       
            
        end    
    end

    methods (Access=private)               
        function obj = runKernel(obj, lin, slc, avg, doPlayRF)

            %% Input check
            if not(exist('doPlayRF','var'))
                doPlayRF = true;
            end
        
            %% RF and ADC settings
            if doPlayRF
                obj.rf.freqOffset = obj.gz.amplitude * obj.thickness * (slc-1-(obj.nSlices-1)/2);
                obj.rf.phaseOffset = obj.rf_phase/180*pi;
                obj.adc.phaseOffset = obj.rf_phase/180*pi;
                obj.rf_inc = mod(obj.rf_inc + obj.rfSpoilingInc, 360.0);
                obj.rf_phase = mod(obj.rf_phase + obj.rf_inc, 360.0);
                obj.seq.addBlock(obj.rf, obj.gz, mr.makeLabel('SET','NAV',false), mr.makeLabel('SET','AVG',avg-1));
            else
                obj.rf.freqOffset = 0;
                obj.rf.phaseOffset = 0;
                obj.seq.addBlock(obj.gz, mr.makeLabel('SET','NAV',true), mr.makeLabel('SET','AVG',avg-1));
            end
            
            %% Slice refocusing gradient
            obj.seq.addBlock(obj.gzReph);
        
            %% External trigger and gradient-free interval a
            obj.seq.addBlock(obj.extTrigger, mr.makeDelay(obj.fillTE(1)));
        
            %% Read-prewinding and phase encoding gradients
            gyPre = mr.makeTrapezoid('y', ...
                    'Area', obj.phaseAreas(lin), ...
                    'Duration', mr.calcDuration(obj.gxPre), ...
                    'system',obj.sys);
            obj.seq.addBlock(obj.gxPre,gyPre);

            %% All LABELS / counters an flags are automatically initialized to 0 in the beginning, no need to define initial 0's  
            % so we will just increment LIN after the ADC event (e.g. during the spoiler)         
            %seq.addBlock(mr.makeDelay(1)); % older scanners like Trio may need this
            % dummy delay to keep up with timing

            %% Set labels
            labels = [  {mr.makeLabel('SET','LIN', lin-1)}, ...
                        {mr.makeLabel('SET','SLC', slc-1)}, ...        
                        {mr.makeLabel('SET','AVG', avg-1)}];

            %% First readout gradient
            obj.seq.addBlock(obj.gx, obj.adc, mr.makeLabel('SET','ECO', 0), labels{:});
        
            %% Fly back
            obj.seq.addBlock(obj.gxFlyBack);
        
            %% Second readout gradient   
            obj.seq.addBlock(obj.gx, obj.adc, mr.makeLabel('SET','ECO', 1));     
        
            %% Negative Phase encoding
            gyPre.amplitude = -gyPre.amplitude;
        
            %% Spoiling
            spoilBlockContents = {obj.gxSpoil, gyPre, obj.gzSpoil};
            obj.seq.addBlock(spoilBlockContents{:});

            %% Add delay
            obj.seq.addBlock(mr.makeDelay(obj.fillTR));
        
        end
    end
end