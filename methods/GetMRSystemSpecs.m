function specs = GetMRSystemSpecs(scannerType)
% Define scanner specs

% (c) 2024 Skope Magnetic Resonance Technologies AG  

    specs = [];
    specs.maxGrad_unit = 'mT/m';
    specs.maxSlew_unit = 'T/m/s';

    switch scannerType
        case 'Siemens 3T'
            specs.type = scannerType;            
            specs.maxGrad = 40; 
            specs.maxSlew = 200; 
            specs.B0 = 2.98;
        case 'Siemens Terra 7T SC72CD'
            specs.type = scannerType;            
            specs.maxGrad = 40; 
            specs.maxSlew = 200; 
            specs.B0 = 6.98;
        case 'Siemens 9.4T SC72CD'
            specs.type = scannerType;            
            specs.maxGrad = 40; 
            specs.maxSlew = 200; 
            specs.B0 = 9.385;      
        otherwise
            error('Unknown scanner type');
    end

end