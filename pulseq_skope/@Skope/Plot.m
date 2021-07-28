% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

function Plot(this, range)
% Default PulSeq plot
    
   this.seq.plot('timeRange', range);
end 