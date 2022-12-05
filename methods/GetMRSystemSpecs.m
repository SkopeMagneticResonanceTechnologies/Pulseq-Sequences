function specs = GetMRSystemSpecs(scannerType)
% Define scanner specs

% (c) 2022 Skope Magnetic Resonance Technologies AG  

    specs = [];
    specs.maxGrad_unit = 'mT/m';
    specs.maxSlew_unit = 'T/m/s';

    switch scannerType
        case 'Siemens Terra 7T SC72CD'
            specs.type = scannerType;            
            specs.maxGrad = 40; 
            specs.maxSlew = 200; 
            specs.B0 = 6.98;
        case 'Siemens 3T'
            specs.type = scannerType;            
            specs.maxGrad = 40; 
            specs.maxSlew = 200; 
            specs.B0 = 2.98;
        otherwise
            error('Unknown scanner type');
    end

end