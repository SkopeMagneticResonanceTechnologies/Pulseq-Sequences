classdef skope_localEddyCalib < PulseqBase
% Sequence for local eddy current calibration

% (c) 2024 Skope Magnetic Resonance Technologies AG

    properties        
        axis = ['x', 'y', 'z']      % Axis with blips        
        N_rep = 10                  % Number of repetitions        
        relax = 0.8                 % Relax system limits
        nBlipsPerAxis = 32          % Number of blips per axis
    end

    methods

        function obj = skope_localEddyCalib(seqParams)

            %% Check input structure
            if not(isa(seqParams,'SequenceParams'))
                error('Input need to be a SequenceParams object.');
            end

            %% Set base class properties
            obj.TR = seqParams.TR;       % Repetition time [Unit: s]
            obj.gradFreeTime = 0.5e-3;   % Delay between trigger and blip-train [Unit: s]

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
        
            T_blips_limit = 40e-3; % time limit of blip-train per excitation
        
            if (mod(obj.nBlipsPerAxis, 2))
                error("Choose even number of blips per axis 'nBlipsPerAxis'!"); % balancing condition
            end

            %% Prepare event objects
            % Trigger
            mr_trig = mr.makeDigitalOutputPulse('ext1','duration', obj.sys.gradRasterTime);
            
            % Delay
            mr_gradFreeTime = mr.makeDelay(obj.gradFreeTime);

            % Blips
            mr_blips = {};
            for ax = obj.axis
        
                count = 0;
                for i=[0:(obj.nBlipsPerAxis-1)]
                    
                    sign = (-1)^i;
                    dur = dur_base + count * dur_inc;
                               
                    if (mod(i,2)) 
                        count = count + 1;
                    end
        
                    % use relaxed slew-rate
                    mr_blip_ = obj.makeBlip(ax, obj.sys.maxSlew * obj.relax, dur, sign);
                    assert(dur - mr.calcDuration(mr_blip_) < 1e-15); % sanity check
                    
                    mr_blips = [mr_blips(:)', {mr_blip_}];
                end
            end

            T_blips = PulseqBase.CalculateDuration(mr_blips);
            
            % Delay
            mr_delay = mr.makeDelay(obj.TR - (T_blips + mr_trig.duration));
            
        
            assert(T_blips_limit - (T_blips + mr_trig.duration) > 0);
            assert(abs(obj.TR - (mr_trig.duration + T_blips + mr_delay.delay)) < 1e-8);
            
            % Required by PulSeq IDEA
            mr_dummy = PulseqBase.makeDummy(); 

            %% Combines event objects to event blocks
            T_tot = 0; % counter: total time
        
            for i=[1:obj.N_rep]
        
                % trigger block
                obj.seq.addBlock(mr_trig);
                T_tot = T_tot + mr_trig.duration;
                
                % delay
                obj.seq.addBlock(mr_gradFreeTime);
                T_tot = T_tot + mr_gradFreeTime.delay;
        
                % blip-train
                for mr_blip_=mr_blips
                    obj.seq.addBlock(mr_blip_);
                    T_tot = T_tot + mr.calcDuration(mr_blip_);
                end
        
                % delay
                obj.seq.addBlock(mr_delay);
                T_tot = T_tot + mr_delay.delay;
            end
        
            assert(obj.N_rep * obj.TR - T_tot < 1e-15);
            
            obj.seq.addBlock(mr_dummy);

            %% Number of triggers
            obj.nTrig = obj.N_rep;

            %% Prepare sequence export
            obj.seq.setDefinition('Name', 'lec');
            obj.seq.setDefinition('CameraNrDynamics', obj.nTrig);  
            obj.seq.setDefinition('CameraNrSyncDynamics', 0); 
            obj.seq.setDefinition('CameraAcqDuration', 0.021);  
            obj.seq.setDefinition('CameraInterleaveTR', 0.15); 
            obj.seq.setDefinition('CameraAqDelay', 0); 
            
            %% Write to Pulseq file
            if not(isfolder('exports'))
                mkdir('exports')
            end
            obj.seq.write('exports/skope_localEddyCalib.seq')  
                     
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