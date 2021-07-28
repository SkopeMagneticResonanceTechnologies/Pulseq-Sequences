% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

function Prepare_local_offres_calib(this)
% Prepare local off-resonance and position calibration

    %% Set parameters
    this.seq_params.axis = ['x', 'y', 'z'];
    this.seq_params.T_flattop = 1100e-3; % [s]
    grad_amp_ = 2.5; % [mT/m]
    this.seq_params.grad_amp = mr.convert(grad_amp_, 'mT/m', 'Hz/m', 'gamma', this.sys.gamma);
    this.seq_params.T_trig_delay = 990e-3; % trigger delay [s]
    this.seq_params.T_inter = 100e-3; % delay between event-blocks [s]
    
    %% Prepare event objects and combine to eventblocks
    area_flattop = this.seq_params.grad_amp * this.seq_params.T_flattop; 

    mr_trig = mr.makeDigitalOutputPulse('ext1','duration', this.sys.gradRasterTime, 'delay', this.seq_params.T_trig_delay);
    mr_inter = mr.makeDelay(this.seq_params.T_inter);
    
    % no-grad
    mr_G_ = mr.makeTrapezoid('x', 'FlatTime', this.seq_params.T_flattop, 'FlatArea', 0); % could be replaced by delay
    this.seq.addBlock(mr_trig, mr_G_);
    this.seq.addBlock(mr_inter);

    for ax = this.seq_params.axis
        
        if(ax == 'x')
            % Hint by Bonn-Group: On their scanner, PulSeq flips x-axis when 
            % transforming from physical to logical coordinate system
            area_flattop_ = area_flattop * (-1);
        else
            area_flattop_ = area_flattop;
        end
        
        mr_G_ = mr.makeTrapezoid(ax, 'FlatTime', this.seq_params.T_flattop, 'FlatArea', area_flattop_);
        this.seq.addBlock(mr_trig, mr_G_);
        this.seq.addBlock(mr_inter)
    end
    
    % Required by PulSeq IDEA
    this.seq.addBlock(this.Make_dummy);
end 