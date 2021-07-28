% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

function Prepare(this, varargin)
% Prepare function wrapper

   if nargin > 1
       this.seq_params = varargin{1};
   end
   
   switch(this.seq_type)
       case("Local Eddy-current Calibration")
           this.Prepare_local_ec_calib();
    
       case("Local Off-res Calibration")
           this.Prepare_local_offres_calib();
            
       otherwise
           error("Unknown seq_type");  
   end
   
   this.seq_params.total_time = sum(this.seq.blockDurations);
end 