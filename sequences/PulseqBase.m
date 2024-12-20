classdef (Abstract) PulseqBase < handle
% Base class for Pulseq sequences that support field monitoring

% (c) 2022 Skope Magnetic Resonance Technologies AG

    properties
       
        % System properties
        sys

        % Array of echo times (1x2) [Unit: s]
        TE

        % Excitation repetition time [Unit: s]
        TR

        % Duration of scanner readout event [Unit: s]
        readoutTime

        % Field of view (1x2) [Unit: m]
        fov

        % Number of readout samples
        Nx

        % Number of phase encoding lines
        Ny

        % Number of partition encoding steps
        Nz

        % Slice thickness
        thickness

        % Number of slices
        nSlices

        % Number of synchronization triggers
        nSyncDynamics = 10

        % Flip angle
        alpha = 90

        % Number of repetitions
        nRep = 1;

        % Number of averages
        nAve = 1;

        % Slice orientation
        sliceOrientation = SliceOrientation.TRA;

        % Phase encoding direction
        phaseEncDir = PhaseEncodingDirection.AP;

        % Order of axes according to slice orientation and phase encoding
        % direction
        axesOrder;

        % Sign flip of axis
        axesSign;

        % Number of dummy pulses to reach steady state
        nDummy = 0;
        
    end

     properties(SetAccess=protected, GetAccess=public)

        % Pause between synchronization and actual imaging scans (default 4s)
        preScanPause = 4

        % Minimal camera TR [Unit: s]
        minCameraTR = 115e-3      

        % Number of imaging triggers
        nTrig = 0

        % Maximal duty cycle of Acquisition system
        maxDutyCycleAQSys = 0.25
      
        % Time from rising edge of trigger to first scanner ADC sample [Unit: s]
        triggerToScannerAcqDelay

        % Duration of camera acquisition [Unit: s]
        cameraAcqDuration

        % Number of trigger to skip by AQ system
        skipFactor = 1;

        % Camera interleave TR or blank time [Unit: s]
        % Additional external triggers from the scanner are ignored by the
        % Camera Acquisition System during this period
        cameraInterleaveTR = 0

        % Pulseq object
        seq
    end

    properties (Access=protected)

        % Gradient-free time for probe excitation [Unit: s]
        gradFreeTime = 200e-6

        % Trigger trigger object
        extTrigger

        % Fill time for echo times [Unit: s]
        fillTE

        % Fill time for repetition time [Unit: s]
        fillTR
    end

    methods

        function obj = plot(obj, timeRange)
        % Plot label information and k-space trajectory data

            if not(exist('timeRange','var'))
                timeRange = [0 25]*obj.TR;
            end

            %% plot sequence and k-space diagrams
            % Check if LIN is a used label
            labelLINPresent = false;
            for b=1:length(obj.seq.blockEvents)
                block = obj.seq.getBlock(b);
                if isfield(block,'label') && any(cellfun(@(x)strcmp(x,'LIN'),{block.label.label}))
                    labelLINPresent = true;
                    break;
                end
            end
            if labelLINPresent
                obj.seq.plot('timeRange', timeRange, 'TimeDisp', 'ms', 'label', 'lin');
            else
                obj.seq.plot('timeRange', timeRange, 'TimeDisp', 'ms');
            end

            if false                
                obj.seq.plot('timeRange', timeRange, 'TimeDisp', 'ms', 'label', 'eco');
                obj.seq.plot('timeRange', timeRange, 'TimeDisp', 'ms', 'label', 'slc');   
                obj.seq.plot('timeRange', timeRange, 'TimeDisp', 'ms', 'label', 'nav');
                obj.seq.plot('timeRange', timeRange, 'TimeDisp', 'ms', 'label', 'avg');
            end

             if not(isempty(obj.seq.labelincLibrary.data)) || not(isempty(obj.seq.labelsetLibrary.data))                
            
                % k-space trajectory calculation
                [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = obj.seq.calculateKspacePP();
                
                % plot k-spaces
                figure; 
                subplot(2,1,1);
                plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
                hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
                title('k-space components as functions of time');
                drawnow
                subplot(2,1,2);
                plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
                axis('equal'); % enforce aspect ratio for the correct trajectory display
                hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
                title('2D k-space');
                drawnow
            end
        end

        function obj = test(obj)
        % Very optional slow step but useful for testing during
        % development e.g. for the real TE, TR or for staying within
        % slew rate limits
            rep = obj.seq.testReport;
            fprintf([rep{:}]);

        end
                
        function obj = CalculateInterleaveTR(obj, triggerTR)
        % Calculate blank time for Camera Acquisition System. The field
        % probes can perform an acquisition circa every 110 ms. If the
        % triggering is faster, several triggers need to be skipped. The
        % interleaveTR/blank time allows to set the period in which triggers
        % should be ignored by the Camera Acquisition System.
	        
            % Default values
            eps = 2e-3;            

            %% Calculate minimal interleave TR
	        if triggerTR < obj.minCameraTR + eps
                obj.skipFactor = max(2, ceil(obj.minCameraTR / triggerTR));
		        obj.cameraInterleaveTR = obj.skipFactor * triggerTR - eps;	        
            else 
                obj.skipFactor = 1;
		        obj.cameraInterleaveTR = triggerTR - eps;
            end

	        %% Calculate duty cycle
            dutyCycle = obj.cameraAcqDuration /  obj.cameraInterleaveTR;
            if dutyCycle > obj.maxDutyCycleAQSys

		        % Minimal TR that Skope AQ system will allow
		        minTR = obj.cameraAcqDuration / obj.maxDutyCycleAQSys;

		        % Time to add to reach minTR
		        diff = minTR - obj.cameraInterleaveTR;

		        % Number of additional triggers to skip
		        addTrig2Skip = ceil(diff / triggerTR);

		       % Add a multiple of trigger TR
               obj.cameraInterleaveTR =  obj.cameraInterleaveTR + addTrig2Skip * triggerTR;
               obj.skipFactor = obj.skipFactor + addTrig2Skip;
            end
        end

        function t = roundUpToGRT(obj,t)
        % Round time to gradient raster time

            % First round to microseconds to avoid errors due to number representation
            t_us = round(t*1e6);
            gdt_us = obj.sys.gradRasterTime*1e6;

            % Round up to GRT
            t = ceil(t_us/gdt_us)*obj.sys.gradRasterTime;
        end

        function obj = addBlock(obj,varargin)
            for i=1:numel(varargin)
                switch varargin{i}.type
                    case {'trap','grad'}
                        ind = find(cellfun(@(x)strcmp(x,varargin{i}.channel),obj.axesOrder));                     
                        if isfield(varargin{i},'amplitude')
                            varargin{i}.amplitude = obj.axesSign(ind)*varargin{i}.amplitude;
                        elseif isfield(varargin{i},'waveform')                         
                            varargin{i}.waveform = obj.axesSign(ind)*varargin{i}.waveform;
                            varargin{i}.first = varargin{i}.waveform(1);
                            varargin{i}.last = varargin{i}.waveform(end);
                        end
                    otherwise
                        % Nothing to do
                end
                
            end
            obj.seq.addBlock(varargin);
        end
    end

    methods (Static)
        function t = CalculateDuration(set)
        % Calculate the total duration of all event objects contained in a set
            t = 0;            
            for i = set
                t = t + mr.calcDuration(i);
            end      
        end 

        function dummy = makeDummy()
            % Prepare a gradient with 0 amplitude
            % to create a [SHAPES] paragraph required by PulSeq IDEA
            dummyWaveForm = [0 0];
            dummy = mr.makeArbitraryGrad('x', dummyWaveForm);      
        end
    end
end