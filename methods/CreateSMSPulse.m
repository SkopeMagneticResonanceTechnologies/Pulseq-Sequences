function [rf, gz, gzAmplitude, t_rfCenter] = CreateSMSPulse(alpha, slThick, tbw, dur, nSlices, sliceSep, sys, varargin)
% Create simultaneous multi-slice (SMS) RF pulse as Pulseq events.
%
% Designs a single-slice SLR pulse, then superposes frequency-shifted,
% phase-cycled copies to produce the multi-band waveform (Wong ISMRM 2012
% p2209). Outputs are Pulseq mr.makeArbitraryRf / mr.makeArbitraryGrad
% events ready for a Siemens .seq file.
%
% Usage:
%   [rf, gz, gzAmplitude, t_rfCenter] = createsmspulse(alpha, slThick, tbw, dur, ...
%                                               nSlices, sliceSep, sys)
%   [rf, gz, gzAmplitude, t_rfCenter] = createsmspulse(..., 'Name', Value, ...)
%
% Required inputs
%   alpha      [1]   flip angle (deg)
%   slThick    [1]   slice thickness (m)
%   tbw        [1]   time-bandwidth product (integer, e.g. 6 or 8)
%   dur        [1]   RF pulse duration (s)
%   nSlices    [1]   multi-band factor (1-8)
%   sliceSep   [1]   centre-to-centre slice separation (m)
%   sys        struct   Pulseq system struct from mr.opts()
%                Fields used: gamma, maxGrad, maxSlew, rfRasterTime,
%                             gradRasterTime, rfDeadTime
%
% Options (name-value pairs)
%   'type'        'st' (default) | 'ex' | 'se' | 'sat' | 'inv'
%                  SLR pulse type (passed to dzrf).
%   'ftype'       'ls' (default) | 'min' | 'max' | 'pm' | 'ms'
%                  SLR filter type (passed to dzrf).
%   'noRfOffset'  false (default) | true
%                  Place all sub-pulses at isocentre (useful for mb=1 debug).
%   'doSim'       false (default) | true
%                  Run Bloch simulation and plot slice profile.
%
% Outputs
%   rf           Pulseq rf event (mr.makeArbitraryRf)
%   gz           Pulseq gradient event (mr.makeArbitraryGrad, z-axis)
%   gzAmplitude  [1]  slice-select gradient plateau amplitude (Hz/m); use as
%                     obj.gzSMSAmplitude * slicePosition to compute the RF
%                     frequency offset, analogous to obj.gz.amplitude.
%   t_rfCenter   [1]  time from rf event start to RF pulse centre (s)
%
% Dependencies
%   + dzrf (John Pauly SLR toolbox, +jpauly package). The folder that
%     contains +jpauly must be on the MATLAB path, e.g.:
%       addpath('<repo>/sequence/toppe/+toppe/+utils/+rf');
%     so that the call toppe.utils.rf.jpauly.dzrf works, or restructure
%     the +jpauly folder as a standalone package on the path.
%   + Pulseq MATLAB toolbox (mr.makeArbitraryRf, mr.makeArbitraryGrad)
%   + MATLAB Signal Processing Toolbox (resample)
%
% Example
%   sys = mr.opts('maxGrad', 28, 'gradUnit', 'mT/m', ...
%                 'maxSlew', 150, 'slewUnit', 'T/m/s', ...
%                 'rfDeadTime', 100e-6, 'rfRingdownTime', 60e-6, ...
%                 'adcDeadTime', 40e-6);
%   [rf, gz] = getsmspulse(70, 5e-3, 6, 8e-3, 4, 20e-3, sys, ...
%                          'type', 'st', 'doSim', true);

% MIT License
% 
% Copyright (c) 2023 Jon-Fredrik Nielsen, <jfnielse@umich.edu>
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% ---------------------------------------------------------------------------
% Test / self-check call
% ---------------------------------------------------------------------------
if ischar(alpha) && strcmp(alpha, 'test')
    sub_test();
    rf = []; gz = []; gzAmplitude = []; t_rfCenter = [];
    return;
end


% ---------------------------------------------------------------------------
% Parse name-value options
% ---------------------------------------------------------------------------
arg.type        = 'st';
arg.ftype       = 'ls';
arg.noRfOffset  = false;
arg.doSim       = false;
arg = sub_parseargs(arg, varargin);

% ---------------------------------------------------------------------------
% Constants and unit conversions
% ---------------------------------------------------------------------------
GAMMA_HZ_G  = 4.2576e3;     % Hz/G  (= 42.576e3 Hz/mT = 42.576e6 Hz/T)
DT_MS       = 4e-3;         % ms  — internal GE 4 µs design raster
DT_S        = DT_MS * 1e-3; % s

slThick_cm  = slThick  * 100;   % m  → cm
dur_ms      = dur      * 1e3;   % s  → ms
sliceSep_cm = sliceSep * 100;   % m  → cm

% Derive gradient limits for SLR design from the Pulseq sys struct.
% sys.maxGrad [Hz/m], sys.maxSlew [Hz/m/s], sys.gamma [Hz/T]
mxg_Gcm   = sys.maxGrad  / sys.gamma * 100;       % G/cm
mxs_Gcmms = sys.maxSlew  / sys.gamma * 100 / 1e3; % G/cm/ms

% ---------------------------------------------------------------------------
% Design single-slice (base) SLR pulse on 4 µs raster
% ---------------------------------------------------------------------------
[rf1, gz1, gPlateau] = sub_makeslr(alpha, slThick_cm, tbw, dur_ms, ...
    mxg_Gcm, mxs_Gcmms, DT_MS, arg.type, arg.ftype);

% ---------------------------------------------------------------------------
% Locate centre of RF pulse (peak of |rf|)
% ---------------------------------------------------------------------------
I         = find(abs(rf1) > max(abs(rf1(:))) - eps);
iRfCenter = mean(I);

N = length(rf1);
t = ((1:N)' * DT_S) - (iRfCenter * DT_S);   % time axis, centred on RF peak
t_rfCenter = iRfCenter * DT_S;

% ---------------------------------------------------------------------------
% Build SMS waveform: sum of phase-cycled, frequency-shifted sub-pulses
% ---------------------------------------------------------------------------
PHS   = sub_smsphase(nSlices);    % Wong 2012 phase table (rad), Table 1 p2209
rfSMS = zeros(N, 1);

for sl = 1:nSlices
    if arg.noRfOffset
        slOff_cm = 0;
    else
        slOff_cm = (-nSlices/2 + sl - 1) * sliceSep_cm;  % cm
    end
    % f [Hz] = gamma [Hz/G] * G [G/cm] * z [cm]  (units cancel correctly)
    f     = GAMMA_HZ_G * gPlateau * slOff_cm;
    rfSMS = rfSMS + rf1 .* exp(1i*2*pi*f*t) * exp(1i*PHS(sl));
end

% Slice-select gradient plateau amplitude in Pulseq units (Hz/m), returned
% to caller so it can be used as: rfSMS.freqOffset = gzAmplitude * slicePosition
gzAmplitude = gPlateau * GAMMA_HZ_G * 100;   % Hz/m  (×100: G/cm → Hz/m)

% ---------------------------------------------------------------------------
% Optional Bloch simulation / slice-profile display
% ---------------------------------------------------------------------------
if arg.doSim
    fov = 20;
    z   = -fov/2 : 0.05 : fov/2;   % cm
    figure;
    % Note that we simulate here for a negative gz gradient because that's
    % how we will play out the gradient for transversal imaging
    sub_slicesim([0 0 1], rfSMS, -gz1, DT_MS, z, 1000, 100, true);
    title(sprintf('SMS slice profile for transversal (mb=%d, sep=%.1f cm, type=%s)', ...
        nSlices, sliceSep_cm, arg.type));
end

% ---------------------------------------------------------------------------
% Convert to Pulseq events (Siemens units)
% ---------------------------------------------------------------------------

% RF: Gauss → Hz, resample to sys.rfRasterTime
rfp = sub_rf2pulseq(rfSMS, DT_S, sys.rfRasterTime);

% Gradient: G/cm → Hz/m, resample to sys.gradRasterTime
gzp = sub_g2pulseq(gz1, DT_S, sys.gradRasterTime);

% --- Trim leading zeros (keep on gradRaster boundary) ---
I      = find(abs(rfp) > 0, 1, 'first');
if isempty(I), I = 1; end
nDelay = (I-1) - mod(I-1, round(sys.gradRasterTime / sys.rfRasterTime));
rfp    = rfp(nDelay+1 : end);
delay  = nDelay * sys.rfRasterTime;

% --- Trim trailing zeros ---
J   = find(abs(rfp) > 0, 1, 'last');
if ~isempty(J)
    rfp = rfp(1:J);
end

% --- Zero-pad to nearest gradRaster boundary ---
wavdur  = numel(rfp) * sys.rfRasterTime;
ttarget = ceil(wavdur / sys.gradRasterTime) * sys.gradRasterTime;
rfp     = [rfp(:); zeros(round((ttarget - wavdur) / sys.rfRasterTime), 1)];

% --- Enforce RF dead time ---
if delay < sys.rfDeadTime
    gdelay = sys.rfDeadTime - delay;
    delay  = sys.rfDeadTime;
else
    gdelay = 0;
end

% --- Create Pulseq events ---
% mr.makeArbitraryRf expects the integrated flip in rad*Hz (= rad/s / 2pi * ... )
% It internally does: signal = signal/abs(sum(signal*dt)) * flip/(2*pi)
% We pass flip_rad * |sum(rfp * dt)| * 2*pi so that the amplitude is preserved.
flip_rad = alpha / 180 * pi;
rf = mr.makeArbitraryRf(rfp, ...
    flip_rad * abs(sum(rfp * sys.rfRasterTime)) * (2*pi), ...
    'delay',  delay, ...
    'system', sys);
rf.signal = rf.signal / max(abs(rf.signal)) * max(abs(rfp));  % preserve amplitude

gz = mr.makeArbitraryGrad('z', gzp, sys, 'delay', gdelay);
gz.first = 0;
gz.last  = 0;

end % createsmspulse


% ===========================================================================
%  LOCAL FUNCTIONS  (no TOPPE dependency)
% ===========================================================================

% ---------------------------------------------------------------------------
% sub_makeslr   —   design SLR pulse + balanced slice-select gradient
%
%   All units: G, cm, ms
%   Dependency: dzrf from the +jpauly package (John Pauly SLR toolbox)
% ---------------------------------------------------------------------------
function [rf1, gz1, gPlateau] = sub_makeslr(flip, slthick, tbw, dur, mxg, mxs, dt, type, ftype)

GAMMA_KHZ_G = 4.2576;    % kHz/G  (= 4257.6 Hz/G)
GAMMA_HZ_G  = 4.2576e3;  % Hz/G

% --- RF waveform (SLR design) ---
res = round(dur / dt);
dur = res * dt;         % snap duration to raster

nrf    = 200;           % design points before resampling
rfBase = dzrf(nrf, tbw, type, ftype);   % John Pauly toolbox — only external dep.
rfBase = real(rfBase);
rfBase = resample(rfBase, res, nrf);    % Signal Processing Toolbox

% Scale to Gauss: sum(rf) = flip angle in rad
rfBase = flip/180*pi * rfBase(:) / sum(rfBase);
rfBase = rfBase / (GAMMA_KHZ_G * dt * 2*pi);  % Gauss

% Even-length
rfBase = [rfBase; zeros(mod(numel(rfBase), 2), 1)];
npix   = numel(rfBase);

% Index of RF pulse centre (peak)
I    = find(abs(rfBase) > max(abs(rfBase)) - eps);
iref = round(mean(I));

% --- Slice-select gradient trapezoid ---
bw       = tbw / dur;                    % kHz
gPlateau = bw / (GAMMA_KHZ_G * slthick);  % G/cm

if gPlateau > mxg
    error(['sub_makeslr: gPlateau (%.2f G/cm) > mxg (%.2f G/cm).\n', ...
           'Increase pulse duration or reduce slice thickness.'], gPlateau, mxg);
end

s        = mxs * dt * 0.995;  % max G/cm change per raster sample
gss_ramp = s : s : gPlateau;
if isempty(gss_ramp), gss_ramp = 0; end
% Fix boundary when gPlateau is not an exact multiple of s
if gPlateau - gss_ramp(end) > s
    gss_ramp = [gss_ramp, (gPlateau + gss_ramp(end))/2];
end
gss_plat = gPlateau * ones(1, npix);
gss_trap = [gss_ramp, gss_plat, fliplr(gss_ramp)];

iref = iref + numel(gss_ramp);   % adjust for ramp prepended to plateau

% --- Rephaser ---
% No prephaser (balancing lobe) is included here; the SMS pulse is
% embedded in a sequence where the prephaser is handled externally.
% This matches the behaviour of the old getsmspulse, which called
% makeslr with nSpoilCycles=1e-6 (>0) to suppress the prephaser.
switch type
    case {'ex','st','sat'}
        arearep = (gss_trap(iref)/2 + sum(gss_trap((iref+1):end))) * dt * 1e-3;  % G/cm·s
        gzrep   = -sub_trapwave2(arearep, mxg, mxs, dt);
    otherwise   % 'se', 'inv': no rephaser
        arearep = 0;
        gzrep   = [];
end

% --- Assemble (no prephaser) ---
gex  = [gss_trap, gzrep(:)'];
rf1  = [zeros(numel(gss_ramp),1); ...
        rfBase; ...
        zeros(numel(gss_ramp) + numel(gzrep), 1)];

% --- 'ex' only: refine rephaser for flat phase profile across slice ---
if strcmp(type, 'ex')
    Z  = linspace(-0.5*slthick/2, 0.5*slthick/2, 50);
    m  = sub_slicesim([0 0 1], rf1(:), gex(:), dt, Z, 1000, 100, false);
    ph = angle(m);
    P  = polyfit(Z, ph, 1);
    extraarea = P(1) / (2*pi*GAMMA_HZ_G);   % G/cm·s
    gzrep = -sub_trapwave2(arearep - extraarea, mxg, mxs, dt);
    gex   = [gss_trap, gzrep(:)'];
end

% --- 'se' only: bookend with crushers ---
% For SMS the 'se' type is rarely used; crushers are omitted here.
% To add them: gex = [gcrush; gex(:); gcrush]; rf1 = [zeros; rf1; zeros].

% --- Finalize: bookend with zeros, equal length, pad to 4-sample boundary ---
rf1 = [0; rf1(:); 0];
gex = [0; gex(:); 0];
n   = max(numel(rf1), numel(gex));
rf1 = [rf1(:); zeros(n - numel(rf1), 1)];
gex = [gex(:); zeros(n - numel(gex), 1)];

rf1 = sub_makeGElength(rf1);
gz1 = sub_makeGElength(gex);

end % sub_makeslr


% ---------------------------------------------------------------------------
% sub_trapwave2   —   gradient trapezoid with a given area
%
%   Inputs  area  G/cm·s  (signed OK)
%           mxg   G/cm
%           mxs   G/cm/ms
%           dt    ms (raster time)
%   Output  waveform  [1 n] G/cm  (starts and ends at 0)
%
%   Adapted from toppe.utils.trapwave2 (no external dependency)
% ---------------------------------------------------------------------------
function waveform = sub_trapwave2(area, mxg, mxs, dt)

area_ms = area * 1e3;   % G/cm·ms (working in ms throughout)

% Preserve sign; work with positive area
if area_ms < 0
    invert = true;
    area_ms = -area_ms;
else
    invert = false;
end

mxg = 0.995 * mxg;
mxs = 0.995 * mxs;

dg      = mxs * dt;        % max ΔG per raster step
tr      = mxg / mxs;       % ramp duration to full amplitude (ms)
Acrit   = mxs * tr^2;      % area of largest triangle (G/cm·ms)

if area_ms <= Acrit
    % Triangle: ramp up then ramp down
    rtime = sqrt(area_ms / mxs);
    n     = ceil(rtime / dt);
    ramp  = 0 : dg : (n * dg);
    waveform = [ramp, fliplr(ramp)];
else
    % Trapezoid: ramp up, plateau, ramp down
    nr       = ceil(tr / dt);
    ramp     = (0:(nr-1)) * mxs * dt;
    areaRamps = 2 * sum(ramp) * dt;
    np       = ceil((area_ms - areaRamps) / mxg / dt);
    plat     = mxg * ones(1, np);
    waveform = [ramp, plat, fliplr(ramp)];
end

% Scale to exact desired area
wavArea = sum(waveform) * dt;
if wavArea < area_ms
    error('sub_trapwave2: cannot achieve requested area (%.4f). Bug in code.', area);
end
waveform = waveform * (area_ms / wavArea);

if invert
    waveform = -waveform;
end

end % sub_trapwave2


% ---------------------------------------------------------------------------
% sub_makeGElength   —   pad waveform to a multiple of 4 samples
% ---------------------------------------------------------------------------
function g = sub_makeGElength(g)
r = mod(size(g, 1), 4);
if r ~= 0
    g = [g; zeros(4 - r, size(g, 2))];
end
end % sub_makeGElength


% ---------------------------------------------------------------------------
% sub_smsphase   —   Wong ISMRM 2012 phase table (Table 1, p2209)
%
%   Input   mb    multiband factor (1..8)
%   Output  PHS   [1 mb] phase values in radians
% ---------------------------------------------------------------------------
function PHS = sub_smsphase(mb)
if mb < 1 || mb > 8
    error('sub_smsphase: mb must be 1..8 (got %d)', mb);
end
P = zeros(8, 8);
P(1, 1)   = 0;
P(2, 1:2) = [0,     pi    ];
P(3, 1:3) = [0,      0.730,   4.602];
P(4, 1:4) = [0,      3.875,   5.940,   6.197];
P(5, 1:5) = [0,      3.778,   5.335,   0.872,   0.471];
P(6, 1:6) = [0,      2.005,   1.674,   5.012,   5.736,   4.123];
P(7, 1:7) = [0,      3.002,   5.998,   5.909,   2.624,   2.528,   2.440];
P(8, 1:8) = [0,      1.036,   3.414,   3.778,   3.215,   1.756,   4.555,   2.467];
PHS = P(mb, 1:mb);
end % sub_smsphase


% ---------------------------------------------------------------------------
% sub_slicesim   —   Bloch simulation of slice-selective pulse
%
%   Inputs
%     m0        [1 3]     initial magnetisation [Mx My Mz]
%     rf        [n 1]     RF waveform (Gauss, complex)
%     gz        [n 1]     gradient waveform (G/cm)
%     dt        [1]       raster time (ms)
%     Z         [nz 1]    spatial positions (cm)
%     T1, T2    [1]       relaxation times (ms)
%     doDisplay bool      plot? (default true)
%   Output
%     m         [nz 1]    complex transverse magnetisation (normalised)
% ---------------------------------------------------------------------------
function m = sub_slicesim(m0, rf, gz, dt, Z, T1, T2, doDisplay)

if nargin < 8, doDisplay = true; end

rf = rf(:);
gz = gz(:);
Z  = Z(:);

Bz    = gz * Z' * 1e-4;   % [nstep nz] Tesla; 1e-4: G→T, cm/cm cancel
nstep = numel(rf);
m     = zeros(numel(Z), 1);

for iz = 1:numel(Z)
    Beff = [real(rf)*1e-4, -imag(rf)*1e-4, Bz(:,iz)];   % [nstep 3] Tesla
    M    = sub_blochsim(m0, Beff, T1, T2, dt, nstep);
    m(iz) = M(end,1) + 1i*M(end,2);
end

if doDisplay
    T = dt * (1:nstep);
    subplot(1,3,1); hold off;
    plot(T, abs(rf)*1e4); hold on;
    plot(T, gz/max(abs(gz)+eps)*max(abs(rf)*1e4)*0.8);
    legend('|rf| (a.u.)', 'gz (scaled)'); xlabel('time (ms)');
    subplot(1,3,2); plot(Z, abs(m)); xlabel('z (cm)'); ylabel('|m_{xy}|');
    subplot(1,3,3); plot(Z, angle(m)); xlabel('z (cm)'); ylabel('\angle m_{xy} (rad)');
end

end % sub_slicesim


% ---------------------------------------------------------------------------
% sub_blochsim   —   Bloch equation simulator (rotation method)
%
%   Inputs
%     Mi     [1 3]        initial magnetisation [Mx My Mz]
%     beff   [nstep 3]    effective field (Tesla), columns: [Bx By Bz]
%     T1,T2  [1]          relaxation times (ms)
%     dt     [1]          raster time (ms)
%     nstep  [1]          number of time steps
%   Output
%     M      [nstep 3]    magnetisation trajectory
%
%   Adapted from D. Noll's blochsim (BME 516), used in the TOPPE toolbox.
% ---------------------------------------------------------------------------
function M = sub_blochsim(Mi, beff, T1, T2, dt, nstep)

gambar = 42.57e3;          % gamma/(2*pi) in kHz/T
gam    = gambar * 2*pi;    % rad/ms/T

T1inv  = dt / T1;          % loss per step (longitudinal recovery fraction)
T2fac  = 1 - dt / T2;      % decay factor per step

beff = beff * (dt * gam);  % convert to rad (dimensionless rotation angle per step)

M      = zeros(nstep, 3);
M(1,:) = Mi;

for lp = 2:nstep
    B    = beff(lp-1, :);
    Bmag = sqrt(sum(B.^2));

    if Bmag == 0
        M(lp,:) = M(lp-1,:);
    else
        Btrans = sqrt(B(1)^2 + B(2)^2);

        ct = B(3) / Bmag;
        st = sqrt(max(1 - ct^2, 0));

        cphi = 1; sphi = 0;
        if Btrans > 0
            cphi = B(1) / Btrans;
            sphi = sqrt(max(1 - cphi^2, 0)) * sign(B(2));
        end

        cpsi = cos(Bmag);
        spsi = sin(Bmag);

        Mx0 = M(lp-1,1);  My0 = M(lp-1,2);  Mz0 = M(lp-1,3);

        Mx1 = cphi*(ct*(cpsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) ...
              + spsi*(cphi*My0-sphi*Mx0)) + st*(ct*Mz0+st*(sphi*My0+cphi*Mx0))) ...
              - sphi*(-spsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) ...
              + cpsi*(cphi*My0-sphi*Mx0));

        My1 = sphi*(ct*(cpsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) ...
              + spsi*(cphi*My0-sphi*Mx0)) + st*(ct*Mz0+st*(sphi*My0+cphi*Mx0))) ...
              + cphi*(-spsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) ...
              + cpsi*(cphi*My0-sphi*Mx0));

        Mz1 = ct*(ct*Mz0 + st*(sphi*My0+cphi*Mx0)) ...
              - st*(cpsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) ...
              + spsi*(cphi*My0-sphi*Mx0));

        M(lp,1) = Mx1 * T2fac;
        M(lp,2) = My1 * T2fac;
        M(lp,3) = Mz1 + (1 - Mz1) * T1inv;
    end
end

end % sub_blochsim


% ---------------------------------------------------------------------------
% sub_rf2pulseq   —   RF: Gauss → Hz, resample to output raster
% ---------------------------------------------------------------------------
function rfOut = sub_rf2pulseq(rf, rasterIn, rasterOut)
GAMMA_HZ_G = 4.2576e3;
rfHz  = rf(:) * GAMMA_HZ_G;
dur   = numel(rfHz) * rasterIn;
ttIn  = (1:numel(rfHz))' * rasterIn - rasterIn/2;
ttOut = (rasterOut/2 : rasterOut : dur)';
rfOut = interp1(ttIn, rfHz, ttOut, 'linear', 'extrap');
rfOut = rfOut(:);
end % sub_rf2pulseq


% ---------------------------------------------------------------------------
% sub_g2pulseq   —   gradient: G/cm → Hz/m, resample to output raster
% ---------------------------------------------------------------------------
function gOut = sub_g2pulseq(g, rasterIn, rasterOut)
if rasterOut < rasterIn
    error('sub_g2pulseq: rasterOut must be >= rasterIn');
end
GAMMA_HZ_G = 4.2576e3;
gHzm  = g(:) * GAMMA_HZ_G * 100;   % Hz/m  (×100: cm→m)
ttIn  = (1:numel(gHzm))' * rasterIn - rasterIn/2;
ttOut = (rasterOut/2 : rasterOut : ttIn(end))';
gOut  = interp1(ttIn, gHzm, ttOut);
if any(isnan(gOut))
    error('sub_g2pulseq: NaN after gradient interpolation');
end
gOut = gOut(:);
end % sub_g2pulseq


% ---------------------------------------------------------------------------
% sub_parseargs   —   minimal name-value argument parser (no toolbox needed)
% ---------------------------------------------------------------------------
function arg = sub_parseargs(arg, varargs)
if mod(numel(varargs), 2) ~= 0
    error('createsmspulse: options must be name-value pairs');
end
fields = fieldnames(arg);
for k = 1:2:numel(varargs)
    name = varargs{k};
    val  = varargs{k+1};
    if ~ischar(name)
        error('createsmspulse: option name must be a string');
    end
    if ~any(strcmp(name, fields))
        error('createsmspulse: unknown option ''%s''', name);
    end
    arg.(name) = val;
end
end % sub_parseargs


% ---------------------------------------------------------------------------
% sub_test   —   quick smoke-test (mb=4, type='st', display slice profile)
% ---------------------------------------------------------------------------
function sub_test()
fprintf('createsmspulse self-test: mb=4, st pulse ...\n');

alpha    = 70;
slThick  = 5e-3;    % m
sliceSep = 20e-3;   % m
tbw      = 6;
dur      = 8e-3;    % s
mb       = 4;

sys = mr.opts('maxGrad', 28, 'gradUnit', 'mT/m', ...
              'maxSlew', 150, 'slewUnit', 'T/m/s', ...
              'rfDeadTime',    100e-6, ...
              'rfRingdownTime', 60e-6, ...
              'adcDeadTime',    40e-6);

[rf, gz, gzAmplitude, t_rfCenter] = CreateSMSPulse(alpha, slThick, tbw, dur, mb, ...
    sliceSep, sys, 'type', 'st', 'doSim', true);

fprintf('  gz amplitude = %.1f Hz/m\n', gzAmplitude);
fprintf('  t_rfCenter   = %.3f ms\n', t_rfCenter*1e3);
fprintf('  RF duration  = %.3f ms\n', numel(rf.signal)*sys.rfRasterTime*1e3);
fprintf('createsmspulse self-test complete.\n');

end % sub_test
