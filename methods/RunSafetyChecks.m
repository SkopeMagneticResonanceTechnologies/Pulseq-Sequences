function RunSafetyChecks(seq, specs)
% RunSafetyChecks(obj,specs)
% Executes available safety checks:
%   - PNS
%   - CarNS (if supported by ASC)
%   - Gradient spectrum / forbidden bands

ascfile = fullfile('dependencies','asc', specs.HWfilename);
useASCForbiddenBands = false;
if isfile(ascfile)
    % PNS and CarNS
    fprintf('PNS and CarNS computation: using hardware file %s \n', ascfile);
    [ok, pns_norm, pns_comp, t_axis] = seq.calcPNS(ascfile);
    % PNS
    maxPNS = max(pns_norm(1,:));
    if ok(1)
        fprintf('PNS check passed successfully (%.1f%% < 100%%)\n',100*maxPNS);
    else
        fprintf('PNS check failed (%.1f %% >= 100%%)\n', 100*maxPNS);
    end
    % CarNS
    if size(pns_norm,1) >= 2
        maxCNS = max(pns_norm(2,:));
        if ok(2)
            fprintf('CarNS check passed successfully (%.1f%% < 100%%)\n', 100*maxCNS);
        else
            fprintf('CarNS check failed (%.1f %% >= 100%%)\n', 100*maxCNS);
        end
    else 
        warning('CarNS model not available in hardware .asc file. CarNS safety check was skipped.');
    end
    % Forbidden bands check
    ascData = mr.Siemens.readasc(ascfile);
    if isfield(ascData.asGPAParameters(1).sGCParameters, 'aflAcousticResonanceFrequency') 
        freqs = ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency;
        useASCForbiddenBands = any(freqs > 0);
    end
else
    warning('Hardware .asc file not found in dependencies/asc. PNS/CarNS safety checks were skipped.')
end
% Gradient spectrum check
if useASCForbiddenBands
    fprintf('Gradient Spectrum computation using acoustic resonances from hardware file %s \n', ascfile);
    [R,Rax,F] = seq.gradSpectrum(ascfile);
elseif ~isempty(specs.forbiddenBandsEchoSpacingLimits)
    warning('Acoustic resonance frequencies not found in hardware .asc file. Using scanner forbidden bands from specs.');
    FB = struct([]);
    for ii = 1:size(specs.forbiddenBandsEchoSpacingLimits,1)
        esMin = specs.forbiddenBandsEchoSpacingLimits(ii,1);
        esMax = specs.forbiddenBandsEchoSpacingLimits(ii,2);
        fHigh = 1/esMin;
        fLow  = 1/esMax;
        FB(ii).freq = (fHigh + fLow)/2;
        FB(ii).bw   = fHigh - fLow;
    end
    [R,Rax,F] = seq.gradSpectrum(FB);
else
    warning('Acoustic resonance frequencies not found in hardware .asc file or in scanner specs. Plotting gradient spectrum without them.');
    [R,Rax,F] = seq.gradSpectrum([]);
end