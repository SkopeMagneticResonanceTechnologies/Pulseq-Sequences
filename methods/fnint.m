%% Copyright notice
%
% Retrieved from https://www.mathworks.com/matlabcentral/fileexchange/71225-splinefit
%
% Copyright (c) 2009, Jonas Lundgren
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of none nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function output = fnint(pp,a,b)
%PPINT Integrate piecewise polynomial.
%   QQ = PPINT(PP,A) returns the indefinite integral from A to X of a
%   piecewise polynomial PP. PP must be on the form evaluated by PPVAL.
%   QQ is a piecewise polynomial on the same form. Default value for A is
%   the leftmost break of PP.
%
%   I = PPINT(PP,A,B) returns the definite integral from A to B.
%
%   Example:
%       x = linspace(-pi,pi,7);
%       y = sin(x);
%       pp = spline(x,y);
%       I = ppint(pp,0,pi)
%
%       qq = ppint(pp,pi/2);
%       xx = linspace(-pi,pi,201);
%       plot(xx,-cos(xx),xx,ppval(qq,xx),'r')
%
%   See also PPVAL, SPLINE, SPLINEFIT, PPDIFF
%   Author: Jonas Lundgren <splinefit@gmail.com> 2009

if nargin < 1, help ppint, return, end
if nargin < 2, a = pp.breaks(1); end
% Get coefficients and breaks
coefs = pp.coefs;
[m,n] = size(coefs);
xb = pp.breaks;
pdim = prod(pp.dim);
% Interval lengths
hb = diff(xb);
hb = repmat(hb,pdim,1);
hb = hb(:);
% Integration
coefs(:,1) = coefs(:,1)/n;
y = coefs(:,1).*hb;
for k = 2:n
    coefs(:,k) = coefs(:,k)/(n-k+1);
    y = (y + coefs(:,k)).*hb;
end
y = reshape(y,pdim,[]);
I = cumsum(y,2);
I = I(:);
coefs(:,n+1) = [zeros(pdim,1); I(1:m-pdim)];
% Set preliminary indefinite integral
qq = pp;
qq.coefs = coefs;
qq.order = n+1;
% Set output
if nargin < 3
    % Indefinite integral from a to x
    if a ~= xb(1)
        I0 = ppval(qq,a);
        I0 = I0(:);
        I0 = repmat(I0,m/pdim,1);
        qq.coefs(:,n+1) = qq.coefs(:,n+1) - I0;
    end
    output = qq;
else
    % Definite integral from a to b
    output = ppval(qq,b) - ppval(qq,a);
end
