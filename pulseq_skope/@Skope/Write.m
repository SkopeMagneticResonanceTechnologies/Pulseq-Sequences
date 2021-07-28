% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

function Write(this, varargin)
% Write to file
   
   if (nargin == 1)
       name = varargin{1};
   else
       name = this.seq_name; % default name
   end
   
   this.seq.setDefinition('Name', name);
   this.seq.write(join([name, '.seq'],''));   
end 