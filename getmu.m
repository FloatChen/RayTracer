function [m, v]=getmu(lat)
%[m, v]=mu(lat,ellipseCon);
% Given a latitude (rad), returns mu and v for
% integration of geodesic equations.  
% Uses WGS84 constants.
% B. Dushaw  November 2007
% Inputs:
% lat - latitude in radians
% ellipseCon - Structure containing
% ellipseCon.req    =   6378137.0
% ellipseCon.finv   =   298.257223563
% ellipseCon.eps2   =   2/ellipseCon.finv
global ellipseCon

fac     = 1-ellipseCon.eps2*sin(lat)^2;
sfac    = sqrt(fac);

m   = ellipseCon.req*(1-ellipseCon.eps2)/(fac*sfac);

v   = ellipseCon.req/sfac;
end
