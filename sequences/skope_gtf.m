classdef skope_gtf < PulseqBase
% Sequence for local eddy current calibration

% (c) 2022 Skope Magnetic Resonance Technologies AG

    properties        
        axis = ['x', 'y', 'z']      % Axis with blips        
        N_rep = 1                   % Number of repetitions        
        relax = 0.9                 % Relax system limits
        nBlipsPerAxis = 17          % Number of blips per axis
    end

    methods

        function obj = skope_gtf(scannerType)

            %% Set base class properties
            obj.TR = 5;              % Repetition time [Unit: s]
            obj.gradFreeTime = 0.5e-3;   % Delay between trigger and blip-train [Unit: s]

            %% Get system limits
            specs = GetMRSystemSpecs(scannerType); 

            if not(strcmpi(specs.maxGrad_unit,'mT/m'))
                error('Expected mT/m for maximum gradient.');
            end

            if not(strcmpi(specs.maxSlew_unit,'T/m/s'))
                error('Expected T/m/s for slew rate.');
            end

            %% Used gradient amplitude and slew rate by this sequence
            maxGrad = 40;
            maxSlew = 200;

            %% Check specs
            if maxGrad > specs.maxGrad
                error('Scanner does not support requested gradient amplitude.');
            end
            if maxSlew > specs.maxSlew
                error('Scanner does not support requested slew rate.');
            end

            % Set system limits
            obj.sys = mr.opts('MaxGrad', maxGrad, ...
                              'GradUnit','mT/m',...
                              'MaxSlew', maxSlew, ...
                              'SlewUnit','T/m/s',...
                              'rfRingdownTime', 30e-6, ...
                              'rfDeadtime', 100e-6,...
                              'B0', specs.B0 ...
            ); 

            %% Create a new sequence object
            obj.seq = mr.Sequence(obj.sys);  
            
            %% Set parameters            
            % Compliance to gradient raster time
            dur_base_ = 6e-5; % [s]; blip duration (factor must be even)
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
                for count=[0:(obj.nBlipsPerAxis-1)]
                    
                    sign = 1;
                    dur = dur_base + count * dur_inc;
        
                    % use relaxed slew-rate
                    mr_blip_ = obj.makeBlip(obj.axis(ax), obj.sys.maxSlew * obj.relax, dur, sign);
                    assert(dur - mr.calcDuration(mr_blip_) < 1e-15); % sanity check
                    
                    mr_blips{ax,count+1} = mr_blip_;
                end
            end
            
            % Required by PulSeq IDEA
            mr_dummy = PulseqBase.makeDummy(); 

            %% Combines event objects to eventblocks   
            for i=[1:obj.N_rep]
                for ax = 1:3
                    % blip-train
                    for n=0:size(mr_blips(),2)
    
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
                    
            obj.seq.addBlock(mr_dummy);

            %% Number of triggers
            obj.nTrig = obj.N_rep*3*(obj.nBlipsPerAxis+1);

            %% Prepare sequence export
            obj.seq.setDefinition('Name', 'gtf');
            obj.seq.setDefinition('CameraNrDynamics', obj.nTrig);  
            obj.seq.setDefinition('CameraNrSyncDynamics', 0); 
            obj.seq.setDefinition('CameraAcqDuration', 0.100);  
            obj.seq.setDefinition('CameraInterleaveTR', 0.400); 
            obj.seq.setDefinition('CameraAqDelay', 0); 
            
            %% Write to Pulseq file
            if not(isfolder('exports'))
                mkdir('exports')
            end
            obj.seq.write('exports/skope_gtf.seq')  
                     
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