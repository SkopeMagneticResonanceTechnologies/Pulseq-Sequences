%% Set scanner parameters
% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

function [sys, scanner] = Set_scanner(this, scanner_type)
    scanner = [];
    switch scanner_type
        case "Siemens Terra 7T SC72CD"
            scanner.type = scanner_type;
            
            scanner.maxGrad = 40;
            scanner.maxGrad_unit = 'mT/m';

            scanner.maxSlew = 200;
            scanner.maxSlew_unit = 'T/m/s';
            
        case "Canon Galan 3T"
            % TODO

        case "Philips XX 3T"
            % TODO

        case "Philips XX 7T"
            % TODO
    end
    
    sys = mr.opts(  'MaxGrad',  scanner.maxGrad,  'GradUnit', scanner.maxGrad_unit, ...
                    'MaxSlew',  scanner.maxSlew,  'SlewUnit', scanner.maxSlew_unit);
end