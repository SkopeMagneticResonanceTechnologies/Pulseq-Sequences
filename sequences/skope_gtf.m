classdef skope_gtf < PulseqBase
% Sequence for impulse response measurements using gradient blips

% (c) 2024 Skope Magnetic Resonance Technologies AG

    properties        
        axis = ['x', 'y', 'z']      % Axis with blips            
        relax = 0.9                 % Relax system limits
        nBlipsPerAxis = 19          % Number of blips per axis (one 
                                    % measurement without a gradient blib 
                                    % will be added per axis)
    end

    methods

        function obj = skope_gtf(seqParams)

            %% Check input structure
            if not(isa(seqParams,'SequenceParams'))
                error('Input need to be a SequenceParams object.');
            end

            if not(strcmpi(seqParams.mode,'linearityCheck')) && ...
                 not(strcmpi(seqParams.mode,'default'))
                error('Unknown sequence mode.')
            end

            %% Get system limits
            specs = GetMRSystemSpecs(seqParams.scannerType); 

            if not(strcmpi(specs.maxGrad_unit,'mT/m'))
                error('Expected mT/m for maximum gradient.');
            end

            if not(strcmpi(specs.maxSlew_unit,'T/m/s'))
                error('Expected T/m/s for slew rate.');
            end

            obj.TR = seqParams.TR;   
            obj.gradFreeTime = 0.5e-3;   % Delay between trigger and blip-train [Unit: s]
            obj.nAve = seqParams.nAve;
            obj.nRep = seqParams.nRep;
            
            % Bug fix for Pulseq error in version 1.4.0.
            obj.signFlip = seqParams.signFlip;

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

            %% Create a new sequence object
            obj.seq = mr.Sequence(obj.sys);  
            
            %% Set parameters            
            % Compliance to gradient raster time
            dur_base_ = 2e-5; % [s]; blip duration (factor must be even)
            dur_incr_ = 2e-5; % [s]; increment to adjust blip duration (factor must be even)
            dur_base = obj.roundUpToGRT(dur_base_);
            dur_inc = obj.roundUpToGRT(dur_incr_);
            assert(~mod((dur_base / 2), obj.sys.gradRasterTime));
            assert(~mod((dur_inc / 2), obj.sys.gradRasterTime));
        
            %% Prepare event objects
            % Trigger
            mr_trig = mr.makeDigitalOutputPulse('ext1','duration', obj.sys.gradRasterTime);
            
            % Delay
            mr_gradFreeTime = mr.makeDelay(obj.gradFreeTime);

            % Blips
            mr_blips = {};
            for ax = 1:3
                index = 0;
                for count=[0:(obj.nBlipsPerAxis-1)]
                    
                    sign = 1;
                    dur = dur_base + count * dur_inc;
        
                    % use relaxed slew-rate
                    mr_blip_ = obj.makeBlip(obj.axis(ax), obj.sys.maxSlew * obj.relax, dur, sign);
                    assert(dur - mr.calcDuration(mr_blip_) < 1e-15); % sanity check
                    
                    if strcmpi(seqParams.mode,'linearityCheck')                        
                        % Play out same blip with half the amplitude
                        mr_blips{ax,index+1} = mr.scaleGrad(mr_blip_,0.5);                    
                        mr_blips{ax,index+2} = mr_blip_;
                        index = index+2;
                    else
                        mr_blips{ax,index+1} = mr_blip_;
                        index = index+1;
                    end
                end
            end
            
            % Required by PulSeq IDEA
            mr_dummy = PulseqBase.makeDummy(); 

            %% Combines event objects to event blocks 
            if strcmpi(seqParams.mode,'linearityCheck') 
                % Play out the 'blip' with zero amplitude twice as well
                startVal = -1;
            else
                startVal = 0;
            end

            for r = 1:obj.nRep
                for ax = 1:3
                    % blip-train
                    for n=startVal:size(mr_blips(),2)
                        for i=1:obj.nAve
                            % trigger block
                            obj.seq.addBlock(mr_trig);
                            T_tot = mr_trig.duration;
                            
                            if n>0
                                % delay
                                obj.seq.addBlock(mr_gradFreeTime);
                                T_tot = T_tot + mr_gradFreeTime.delay;            
                            
                                obj.seq.addBlock(mr_blips{ax,n});
                                T_tot = T_tot + mr.calcDuration(mr_blips{ax,n});
                            end
        
                            % Delay
                            assert(obj.TR - T_tot > 0)
                            mr_delay = mr.makeDelay(obj.TR - T_tot);
                    
                            % delay
                            obj.seq.addBlock(mr_delay);
                        end
                    end
                end
            end
                    
            obj.seq.addBlock(mr_dummy);

            %% Number of triggers
            obj.nTrig = obj.nAve*obj.nRep*3*(obj.nBlipsPerAxis+1);
            if strcmpi(seqParams.mode,'linearityCheck')
                obj.nTrig = obj.nTrig * 2;
            end

            %% Prepare sequence export
            if strcmpi(seqParams.mode,'linearityCheck')
                obj.seq.setDefinition('Name', 'gtf_linCheck');
            else
                obj.seq.setDefinition('Name', 'gtf');
            end
            obj.seq.setDefinition('CameraNrDynamics', obj.nTrig);  
            obj.seq.setDefinition('CameraNrSyncDynamics', 0); 
            obj.seq.setDefinition('CameraAcqDuration', 0.040);  
            obj.seq.setDefinition('CameraInterleaveTR', 0.400); 
            obj.seq.setDefinition('CameraAqDelay', 0); 
            
            %% Write to Pulseq file
            if not(isfolder('exports'))
                mkdir('exports')
            end
            if strcmpi(seqParams.mode,'linearityCheck')
                obj.seq.write('exports/skope_gtf_linearityCheck.seq')  
            else
                obj.seq.write('exports/skope_gtf.seq')  
            end
            
                     
        end
    end

    methods (Access = private)
        function blip = makeBlip(obj, axis, slew, duration, sign)
        % Prepare a blip gradient
                  
            rise_time = duration / 2.;
            grad = slew * rise_time;
            
            if (abs(grad) >= obj.sys.maxGrad * obj.relax)
                err_msg = ['Gradient strength exceeds maximum.' , ...
                          ' Try:' , ...
                          ' (1) decreasing number of blips per axis "nBlipsPerAxis", or ', ...
                          ' (2) decreasing duration increment "dur_incr_".'];
                error(err_msg);
            end
                
            if (axis == 'x')
                % Siemens Pulseq interpreter 1.4.0 flips x-axis when 
                % transforming from physical to logical coordinate system
                sign = sign * obj.signFlip;
            end
            
            % gradient moment
            area = sign * (grad * rise_time);
            area = round(area, 12, 'decimal'); % to prevent pulseq from strange rounding
                                               % in makeTrapezoid.m line:113
            
            blip = mr.makeTrapezoid( axis, ...,
                                    'maxSlew', slew, ...
                                    'maxGrad', obj.sys.maxGrad * obj.relax, ...
                                    'Area', area, ...
                                    'system', obj.sys);   
        end 
    end
end