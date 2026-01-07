function specs = GetMRSystemSpecs(scannerType)
% Define scanner specs

% (c) 2026 Skope Magnetic Resonance Technologies AG  

    specs = [];
    specs.maxGrad_unit = 'mT/m';
    specs.maxSlew_unit = 'T/m/s';

    switch scannerType
        % 3T
        case 'Siemens 3T'
            specs.type = scannerType;            
            specs.maxGrad = 40; 
            specs.maxSlew = 200; 
            specs.B0 = 2.98;
        case 'Siemens 3T Cima.X'
            specs.type = scannerType;            
            specs.maxGrad = 200; 
            specs.maxSlew = 200; 
            specs.B0 = 2.98;
        % 7T
        case 'Siemens 7T Terra SC72CD'
            specs.type = scannerType;            
            specs.maxGrad = 40; 
            specs.maxSlew = 200; 
            specs.B0 = 6.98;
        case 'Siemens 7T Terra.X T60'
            specs.type = scannerType;            
            specs.maxGrad = 80; 
            specs.maxSlew = 200; 
            specs.B0 = 6.98;
        % 9.4T
        case 'Siemens 9.4T SC72CD'
            specs.type = scannerType;            
            specs.maxGrad = 40; 
            specs.maxSlew = 200; 
            specs.B0 = 9.385;      
        otherwise
            error('Unknown scanner type');
    end

end