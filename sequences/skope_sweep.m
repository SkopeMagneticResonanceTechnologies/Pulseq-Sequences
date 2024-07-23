classdef skope_sweep < PulseqBase
% Sequence for off-resonance and position calibration

% (c) 2024 Skope Magnetic Resonance Technologies AG

    properties        
        axis = ['x', 'y', 'z'];     % Axis with blips      
        flattopTime = 1100e-3;      % [s]
        gradientAmplitude = 2.5;    % [mT/m]
    end

    methods

        function obj = skope_sweep(seqParams, varargin)

            %% Check input structure
            if not(isa(seqParams,'SequenceParams'))
                error('Input need to be a SequenceParams object.');
            end

            %% Set base class properties
            obj.TR = seqParams.TR;       % Repetition time [Unit: s]
            obj.gradFreeTime = 0.5e-3;   % Delay between trigger and blip-train [Unit: s]

            T_trig_delay = 1e-3; % trigger delay [s]

            % Bug fix for Pulseq error in version 1.4.0.
            obj.signFlip = seqParams.signFlip;

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
            load('inSweep.mat')
            inSweep(end+1) = 0; 
            inSweep = inSweep(403:end)*100;
            amp = max(inSweep);
            grad_amp_Hzm = mr.convert(amp, 'mT/m', 'Hz/m', 'gamma', obj.sys.gamma);            
            T_inter = 2000e-3; % delay between event-blocks [s]
            
            %% Prepare event objects and combine to eventblocks
            mr_trig = mr.makeDigitalOutputPulse('ext1','duration', obj.sys.gradRasterTime, 'delay', T_trig_delay);
            mr_inter = mr.makeDelay(T_inter); 

            nAvg = 50;
            for ax = obj.axis 
                for avg = 1:nAvg                               
                    if(ax == 'x')
                        g = mr.makeArbitraryGrad('x',inSweep*grad_amp_Hzm);
                    elseif (ax == 'y')
                        g = mr.makeArbitraryGrad('y',inSweep*grad_amp_Hzm);
                    else
                        g = mr.makeArbitraryGrad('z',inSweep*grad_amp_Hzm);                    
                    end                    
                    obj.seq.addBlock(mr_trig);                
                    obj.seq.addBlock(g);
                    obj.seq.addBlock(mr_inter);
                end
            end
            
            % Required by PulSeq IDEA
            obj.seq.addBlock(PulseqBase.makeDummy);

            % Number of external triggers
            obj.nTrig = 3*nAvg;

            %% Prepare sequence export
            obj.seq.setDefinition('Name', 'sweep');
            obj.seq.setDefinition('CameraNrDynamics', obj.nTrig);  
            obj.seq.setDefinition('CameraNrSyncDynamics', 0); 
            obj.seq.setDefinition('CameraAcqDuration', 0.1);  
            obj.seq.setDefinition('CameraInterleaveTR', 0.4); 
            obj.seq.setDefinition('CameraAqDelay', 0); 
            
            %% Write to Pulseq file
            if not(isfolder('exports'))
                mkdir('exports')
            end
            obj.seq.write('exports/skope_sweep.seq')       
                      
        end
    end 
end