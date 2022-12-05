classdef skope_offresAndPosCalib < PulseqBase
% Sequence for off-resonance and position calibration

% (c) 2022 Skope Magnetic Resonance Technologies AG

    properties        
        axis = ['x', 'y', 'z'];     % Axis with blips      
        flattopTime = 1100e-3;      % [s]
        gradientAmplitude = 2.5;    % [mT/m]
    end

    methods

        function obj = skope_offresAndPosCalib(scannerType, varargin)

            %% Set base class properties
            obj.TR = 200e-3;             % Repetition time [Unit: s]
            obj.gradFreeTime = 0.5e-3;   % Delay between trigger and blip-train [Unit: s]

            T_trig_delay = 990e-3; % trigger delay [s]

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
            grad_amp_Hzm = mr.convert(obj.gradientAmplitude, 'mT/m', 'Hz/m', 'gamma', obj.sys.gamma);
            
            T_inter = 100e-3; % delay between event-blocks [s]
            
            %% Prepare event objects and combine to eventblocks
            area_flattop = grad_amp_Hzm * obj.flattopTime;         
            mr_trig = mr.makeDigitalOutputPulse('ext1','duration', obj.sys.gradRasterTime, 'delay', T_trig_delay);
            mr_inter = mr.makeDelay(T_inter);
            
            % no-grad
            mr_G_ = mr.makeTrapezoid('x', 'FlatTime', obj.flattopTime, 'FlatArea', 0); % could be replaced by delay
            obj.seq.addBlock(mr_trig, mr_G_);
            obj.seq.addBlock(mr_inter);
        
            for ax = obj.axis                
                if(ax == 'x')
                    % Hint by Bonn-Group: On their scanner, PulSeq flips x-axis when 
                    % transforming from physical to logical coordinate system
                    area_flattop_ = area_flattop * (-1);
                else
                    area_flattop_ = area_flattop;
                end
                
                mr_G_ = mr.makeTrapezoid(ax, 'FlatTime', obj.flattopTime, 'FlatArea', area_flattop_);
                obj.seq.addBlock(mr_trig, mr_G_);
                obj.seq.addBlock(mr_inter)
            end
            
            % Required by PulSeq IDEA
            obj.seq.addBlock(PulseqBase.makeDummy);

            % Number of external triggers
            obj.nTrig = 4;

            %% Prepare sequence export
            obj.seq.setDefinition('Name', 'opc');
            obj.seq.setDefinition('CameraNrDynamics', obj.nTrig);  
            obj.seq.setDefinition('CameraNrSyncDynamics', 0); 
            obj.seq.setDefinition('CameraAcqDuration', 0.1);  
            obj.seq.setDefinition('CameraInterleaveTR', 0.4); 
            obj.seq.setDefinition('CameraAqDelay', 0); 
            
            %% Write to Pulseq file
            if not(isfolder('exports'))
                mkdir('exports')
            end
            obj.seq.write('exports/skope_offresAndPosCalib.seq')       
                      
        end
    end 
end