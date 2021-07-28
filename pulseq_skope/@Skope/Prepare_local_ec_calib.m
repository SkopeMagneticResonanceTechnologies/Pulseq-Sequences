% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

function Prepare_local_ec_calib(this)
% Prepare local eddy-current calibration

    %% Set parameters
    this.seq_params.axis = ['x', 'y', 'z']; % axis with blips
    this.seq_params.TR = 200e-3; % [s]    
    this.seq_params.N_rep = 10;
    this.seq_params.delay_trigTrain = 0.5e-3; % [s] delay between trigger and blip-train
    
    this.seq_params.relax = 0.8; % relax system limits
    
    % Compliance to gradient raster time
    dur_base_ = 6e-5; % [s]; blip duration (factor must be even)
    dur_incr_ = 2e-5; % [s]; increment to adjust blip duration (factor must be even)
    this.seq_params.dur_base = this.Rasterize(dur_base_, "grad");
    this.seq_params.dur_incr = this.Rasterize(dur_incr_, "grad");
    assert(~mod((this.seq_params.dur_base / 2), this.sys.gradRasterTime));
    assert(~mod((this.seq_params.dur_incr / 2), this.sys.gradRasterTime));

    T_blips_limit = 40e-3; % time limit of blip-train per excitation

    if (mod(this.seq_params.N_blips_axis, 2) ~= 0)
        error("Choose even number of blips per axis 'N_blips_axis'!"); % balancing condition
    end

    %% Prepare event objects
    % Trigger
    mr_trig = mr.makeDigitalOutputPulse('ext1','duration', this.sys.gradRasterTime);
    
    % Delay
    mr_delay_trigTrain = mr.makeDelay(this.seq_params.delay_trigTrain);

    % Blips
    mr_blips = {};
    for ax = this.seq_params.axis

        count = 0;
        for i=[0:(this.seq_params.N_blips_axis-1)]
            
            sign = (-1)^i;
            dur = this.seq_params.dur_base + count * this.seq_params.dur_incr;
                       
            if (mod(i,2)) 
                count = count + 1;
            end

            % use relaxed slew-rate
            mr_blip_ = this.Make_blip(ax, this.scanner.maxSlew * this.seq_params.relax, this.scanner.maxSlew_unit, dur, sign);
            assert(dur - mr.calcDuration(mr_blip_) < 1e-15); % sanity check
            
            mr_blips = [mr_blips(:)', {mr_blip_}];
        end
    end

    T_blips = this.Calc_duration_of_set(mr_blips);
    
    % Delay
    mr_delay = mr.makeDelay(this.seq_params.TR - (T_blips + mr_trig.duration));
    

    assert(T_blips_limit - (T_blips + mr_trig.duration) > 0);
    assert(abs(this.seq_params.TR - (mr_trig.duration + T_blips + mr_delay.delay)) < 1e-8);
    
    % Required by PulSeq IDEA
    mr_dummy = this.Make_dummy(); 

    %% Combines event objects to eventblocks
    T_tot = 0; % counter: total time

    for i=[1:this.seq_params.N_rep]

        % trigger block
        this.seq.addBlock(mr_trig);
        T_tot = T_tot + mr_trig.duration;
        
        % delay
        this.seq.addBlock(mr_delay_trigTrain);
        T_tot = T_tot + mr_delay_trigTrain.delay;

        % blip-train
        for mr_blip_=mr_blips
            this.seq.addBlock(mr_blip_);
            T_tot = T_tot + mr.calcDuration(mr_blip_);
        end

        % delay
        this.seq.addBlock(mr_delay);
        T_tot = T_tot + mr_delay.delay;
    end

    assert(this.seq_params.N_rep * this.seq_params.TR - T_tot < 1e-15);
    
    this.seq.addBlock(mr_dummy)

end 