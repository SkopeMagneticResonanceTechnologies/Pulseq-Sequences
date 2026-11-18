classdef skope_grad_heat < PulseqBase
% Sequence that plays 10 gradient lobes at 20 mT/m with 5 ms flat-top,
% alternating sign on x/y/z, then waits 10 s, triggers, and finally waits 1 s.

% (c) 2026 Skope Magnetic Resonance Technologies AG

    properties
        nLobes = 480;
        amplitude = 20; % [mT/m]
        flatTopTime = 5e-3; % [s]
        slewFraction = 0.75;
        waitBeforeTrigger = 10; % [s]
        triggerOutput = 'ext1';
        triggerDuration = 300e-6; % [s]
        postTriggerDelay = 1; % [s]
    end

    methods

        function obj = skope_grad_heat(seqParams)

            %% Check input structure
            if not(isa(seqParams,'SequenceParams'))
                error('Input needs to be a SequenceParams object.');
            end

            %% Set base class properties
            obj.gradFreeTime = 0.5e-3;
            obj.triggerOutput = seqParams.triggerOutput;
            obj.TR = seqParams.TR;

            %% Initialize axis metadata required by PulseqBase.addBlock
            obj.axesOrder = {'x','y','z'};
            obj.axesSign = [1, 1, 1];

            %% Get system limits
            specs = GetMRSystemSpecs(seqParams.scannerType);

            if not(strcmpi(specs.maxGrad_unit,'mT/m'))
                error('Expected mT/m for maximum gradient.');
            end

            if not(strcmpi(specs.maxSlew_unit,'T/m/s'))
                error('Expected T/m/s for slew rate.');
            end

            %% Check requested limits
            if seqParams.maxGrad > specs.maxGrad
                error('Scanner does not support requested gradient amplitude.');
            end
            if seqParams.maxSlew > specs.maxSlew
                error('Scanner does not support requested slew rate.');
            end

            %% Set system limits
            obj.sys = mr.opts('MaxGrad', seqParams.maxGrad, ...
                              'GradUnit','mT/m', ...
                              'MaxSlew', seqParams.maxSlew, ...
                              'SlewUnit','T/m/s', ...
                              'rfRingdownTime', 30e-6, ...
                              'rfDeadtime', 100e-6, ...
                              'B0', specs.B0 ...
            );

            %% Create a new sequence object
            obj.seq = mr.Sequence(obj.sys);

            %% Gradient lobe timing parameters
            % Pulseq gradient amplitudes are specified in Hz/m, not mT/m.
            gradAmplitudeHz = mr.convert(abs(obj.amplitude), 'mT/m', 'Hz/m', 'gamma', obj.sys.gamma);
            riseTime = gradAmplitudeHz / (obj.slewFraction * obj.sys.maxSlew);
            riseTime = max(riseTime, obj.sys.gradRasterTime);
            riseTime = obj.sys.gradRasterTime * ceil(riseTime / obj.sys.gradRasterTime);

            assert(obj.amplitude <= seqParams.maxGrad, 'Requested gradient amplitude exceeds scanner limits.');

            refPulse = mr.makeBlockPulse(60*pi/180, obj.sys, ...
                'Duration', 500e-6, ...
                'use', 'excitation');
            trigger = mr.makeDigitalOutputPulse(obj.triggerOutput, 'duration', obj.triggerDuration);
            adc = mr.makeAdc(1024, obj.sys, 'Dwell', 100e-6);

            %% Play 10 gradient lobes with alternating sign on x, y, and z
            axes = {'x', 'y', 'z'};
            for n = 1:obj.nLobes
                signVal = (-1)^(n+1);
                for axIndex = 1:numel(axes)
                    ax = axes{axIndex};
                    g = mr.makeTrapezoid(ax, ...
                        'amplitude', signVal * gradAmplitudeHz, ...
                        'riseTime', riseTime, ...
                        'flatTime', obj.flatTopTime, ...
                        'system', obj.sys);
                    obj.addBlock(g);
                end
                obj.addBlock(mr.makeDelay(5e-3));
            end

            %% Wait, apply ref pulse, trigger, acquire, and finish
            obj.addBlock(mr.makeDelay(obj.waitBeforeTrigger));
            obj.addBlock(refPulse);
            obj.addBlock(trigger);
            obj.addBlock(adc);
            obj.addBlock(mr.makeDelay(obj.postTriggerDelay));

            % Required by Pulseq IDEA
            obj.addBlock(PulseqBase.makeDummy);

            %% Number of external triggers
            obj.nTrig = 1;

            %% Prepare sequence export
            obj.seq.setDefinition('Name', 'grad_heat');
            obj.seq.setDefinition('CameraNrDynamics', obj.nTrig);
            obj.seq.setDefinition('CameraNrSyncDynamics', 0);
            obj.seq.setDefinition('CameraAcqDuration', 0.1);
            obj.seq.setDefinition('CameraInterleaveTR', 0.4);
            obj.seq.setDefinition('CameraAqDelay', 0);

            %% Write to Pulseq file
            if not(isfolder('exports'))
                mkdir('exports')
            end

            if not(isfolder(strcat('exports/', string(seqParams.scannerType))))
                mkdir(strcat('exports/', string(seqParams.scannerType)))
            end

            filename = strcat('exports/', string(seqParams.scannerType), '/skope_grad_heat');
            obj.seq.write(strcat(filename, '.seq'));
        end
    end
end
