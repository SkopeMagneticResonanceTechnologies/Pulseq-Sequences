% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

function blip = Make_blip(this, axis, slew, slew_unit,  duration, sign)
% Prepare a blip gradient
     
    % convert unit for slew
    slew_ = mr.convert(slew, slew_unit, 'Hz/m', 'gamma', this.sys.gamma);
 
    rise_time = duration / 2.;
    grad = slew_ * rise_time;
    
    if (abs(grad) >= this.sys.maxGrad * this.seq_params.relax)
        err_msg = ['Gradient strength exceeds maximum.' , ...
                  ' Try:' , ...
                  ' (1) decreasing number of blips per axis "N_blips_axis", or ', ...
                  ' (2) decreasing duration increment "dur_incr_".'];
        error(err_msg);
    end
        
    if (axis == 'x')
        % Hint by Bonn-Group: On their scanner, PulSeq flips x-axis when 
        % transforming from physical to logical coordinate system
        sign = sign * (-1);
    end
    
    % gradient moment
    area = sign * (grad * rise_time);
    area = round(area, 12, 'decimal'); % to prevent pulseq from strange rounding
                                       % in makeTrapezoid.m line:113
    
    blip = mr.makeTrapezoid(axis, ...,
                            'maxSlew', slew_, ...
                            'maxGrad', this.sys.maxGrad * this.seq_params.relax, ...
                            'Area', area, ...
                            'system', this.sys);   
end 