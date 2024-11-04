classdef skope_sweep < PulseqBase
% Sequence for off-resonance and position calibration

% (c) 2024 Skope Magnetic Resonance Technologies AG

    properties        
        axis = ['x', 'y', 'z'];     % Axis with blips      
        flattopTime = 1100e-3;      % [s]
        gradientAmplitude = 2.5;    % [mT/m]
    end

    methods

        function obj = skope_sweep(seqParams, sweepWaveform)

            %% Check input structure
            if not(isa(seqParams,'SequenceParams'))
                error('Input need to be a SequenceParams object.');
            end

            %% Set base class properties
            obj.gradFreeTime = 0.5e-3;   % Delay between trigger and blip-train [Unit: s]
            obj.nAve = seqParams.nAve;
            obj.TR = seqParams.TR;

            T_trig_delay = 1e-3; % trigger delay [s]

            %% Axes order
            [obj.axesOrder, obj.axesSign, readDir_SCT, phaseDir_SCT, sliceDir_SCT] ...
                = GetAxesOrderAndSign(obj.sliceOrientation,obj.phaseEncDir);

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
            amp = max(sweepWaveform);
            grad_amp_Hzm = mr.convert(amp, 'mT/m', 'Hz/m', 'gamma', obj.sys.gamma);            
            
            %% Prepare event objects and combine to eventblocks
            mr_trig = mr.makeDigitalOutputPulse('ext1','duration', obj.sys.gradRasterTime, 'delay', T_trig_delay);
            
            mr_gradfree = mr.makeDelay(2e-3);
            
            for ax = obj.axis 
                for avg = 1:obj.nAve   

                    % Add trigger
                    obj.addBlock(mr_trig); 

                    % Add delay after trigger
                    obj.addBlock(mr_gradfree);

                    % Play out waveform
                    g = mr.makeArbitraryGrad(ax,sweepWaveform*grad_amp_Hzm);            
                    obj.addBlock(g);

                    totalTime = mr.calcDuration(mr_trig) + ...
                                mr.calcDuration(mr_gradfree) + ...
                                mr.calcDuration(g);

                    mr_inter = mr.makeDelay(obj.TR - totalTime); 

                    % Wait for next waveform
                    obj.addBlock(mr_inter);
                end
            end
            
            % Required by PulSeq IDEA
            obj.addBlock(PulseqBase.makeDummy);

            % Number of external triggers
            obj.nTrig = 3*obj.nAve;

            %% Prepare sequence export
            obj.seq.setDefinition('Name', 'sweep');
            obj.seq.setDefinition('CameraNrDynamics', obj.nTrig);  
            obj.seq.setDefinition('CameraNrSyncDynamics', 0); 
            obj.seq.setDefinition('CameraAcqDuration', 0.070);  
            obj.seq.setDefinition('CameraInterleaveTR', 0.4); 
            obj.seq.setDefinition('CameraAqDelay', 0); 
            
            %% Write to Pulseq file
            if not(isfolder('exports'))
                mkdir('exports')
            end
            if obj.nAve == 1
                obj.seq.write('exports/skope_sweep_one_ave.seq')  
            else
                obj.seq.write('exports/skope_sweep.seq')  
            end
                                       
        end
    end 
end