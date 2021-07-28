%% PulSeq interface for Skope
% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

classdef Skope < handle
%%%%%%%%%%%%%%
% Basic class to use PulSeq with Skope® systems
%%%%%%%%%%%%%%
    properties
        seq_name;
        seq_type;
        scanner;
        seq_params;
        
        % PulSeq objects
        sys;
        seq;
    end
    
    methods 
        % Constructor
        function obj = Skope(seq_name, seq_type, scanner_type)
            obj.seq_name = seq_name;
            obj.seq_type = seq_type;
            [obj.sys, obj.scanner] = obj.Set_scanner(scanner_type);
            obj.seq = mr.Sequence();
        end
        
    end
end
