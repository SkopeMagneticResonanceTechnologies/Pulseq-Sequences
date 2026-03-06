function specs = GetMRSystemSpecs(scannerType)
% Define scanner specs

% (c) 2026 Skope Magnetic Resonance Technologies AG  

    specs = [];
    specs.maxGrad_unit = 'mT/m';
    specs.maxSlew_unit = 'T/m/s';
    specs.forbiddenBandsEchoSpacingUnits = 's';
    specs.adcSamplesLimit = 8192;
    specs.adcSamplesDivisor = 8;

    switch scannerType
        % 3T       
        case 'Siemens 3T Cima.X'
            specs.type = scannerType;            
            specs.maxGrad = 200; 
            specs.maxSlew = 200; 
            specs.B0 = 2.89;               
            specs.forbiddenBandsEchoSpacingLimits = [0.77e-3, 1.07e-3;... %[1113,344]Hz
                                                     1.62e-3, 1.93e-3 ];  %[567,100] Hz                    
        case 'Siemens 3T Connectom'
            specs.type = scannerType;            
            specs.maxGrad = 300; 
            specs.maxSlew = 200; 
            specs.B0 = 2.89;            
            specs.forbiddenBandsEchoSpacingLimits = [0.66e-3, 1.0e-3;... %[1250,500]Hz
                                                     1.54e-3, 1.83e-3 ]; %[596,100] Hz
        % 5T 
        case 'United 5T uMR Jupiter'
            specs.type = scannerType;            
            specs.maxGrad = 120; 
            specs.maxSlew = 200; 
            specs.B0 = 4.95;            
            specs.forbiddenBandsEchoSpacingLimits = [0.62e-3, 0.70e-3;...
                                                     0.80e-3, 0.88e-3;...
                                                     0.98e-3, 1.05e-3;...
                                                     1.09e-3, 1.17e-3 ];    
        % 7T 
        case 'Siemens 7T Terra SC72CD'
            specs.type = scannerType;            
            specs.maxGrad = 40; 
            specs.maxSlew = 200; 
            specs.B0 = 6.98;
            specs.forbiddenBandsEchoSpacingLimits = [0.63e-3, 0.74e-3;...
                                                     1.20e-3, 1.47e-3 ];

        case 'Siemens 7T Terra.X'
            specs.type = scannerType;            
            specs.maxGrad = 135; 
            specs.maxSlew = 250; 
            specs.B0 = 6.98;
            specs.forbiddenBandsEchoSpacingLimits = [0.357e-3, 0.385e-3;...
                                                     0.400e-3, 0.556e-3;...
                                                     0.833e-3, 1.000e-3;...
                                                     1.266e-3, 1.476e-3 ];
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