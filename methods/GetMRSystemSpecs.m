function specs = GetMRSystemSpecs(scannerType)
% Define scanner specs

% (c) 2026 Skope Magnetic Resonance Technologies AG  

    specs = [];
    specs.maxGrad_unit = 'mT/m';
    specs.maxSlew_unit = 'T/m/s';
    specs.forbiddenBandsEchoSpacingUnits = 's';

    switch scannerType
        % 3T       
        case 'Siemens 3T Cima.X'
            specs.type = scannerType;            
            specs.maxGrad = 200; 
            specs.maxSlew = 200; 
            specs.B0 = 2.89;
            specs.forbiddenBandsEchoSpacingLimits = [0.63e-3, 0.74e-3;...
                                                    1.20e-3, 1.47e-3 ]; % placeholder
        % 7T
        case 'Siemens 7T Terra SC72CD'
            specs.type = scannerType;            
            specs.maxGrad = 40; 
            specs.maxSlew = 200; 
            specs.B0 = 6.98;
            specs.forbiddenBandsEchoSpacingLimits = [0.63e-3, 0.74e-3;...
                                                     1.20e-3, 1.47e-3 ];
        % 9.4T
        case 'Siemens 9.4T SC72CD'
            specs.type = scannerType;            
            specs.maxGrad = 40; 
            specs.maxSlew = 200; 
            specs.B0 = 9.385;    
            specs.forbiddenBandsEchoSpacingLimits = [0.63e-3, 0.74e-3;...
                                                     1.20e-3, 1.47e-3 ];
        otherwise
            error('Unknown scanner type');
    end

end