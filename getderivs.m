function dP = getderivs(P,fr)
% dP = getderivs(P,ellipseCon,fr)
% Calculates derivatives at (lat,lon) and alpha for use with Runge-Kutta integration.
% Does not depend on x, or the time variable.
%
% Input:
% P = [ alpha lat lon ]
% ellipseCon = Structure,
%   ellipseCon.req, 
%   ellipseCon.finv, 
%   ellipseCon.eps2
% fr = frequency [Hz] 
%
% Returns:
% dP = [a1 a2 a3]
% a1:  d(alpha)/dt
% a2:  d(lat)/dt
% a3:  d(lon)/dt
%
% All angles in radians.
%

alpha = P(1);
lat = P(2);
lon = P(3);

% First get the ellipsoid factors for this latitude.
[m, v]=getmu(lat);

% Get speed and derivatives at lat,lon.
[C, D1, D2] = getc(lat,lon,fr);

% minimize the sins, cosines and divisions to be made...
sa=sin(alpha);
ca=cos(alpha);
ivcl=1/(v*cos(lat));
im=1/m;

%%%%%%% EQUATIONS SAME AS MUNK, OREILLY, REID 1988
a1= sa*D1*im - ca*D2*ivcl + C*sa*sin(lat)*ivcl;

a2= C*ca*im;

a3= C*sa*ivcl;

dP = [a1 a2 a3];



end