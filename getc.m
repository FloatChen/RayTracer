function [c, D1, D2] = getc(lat,lon,f)
% function [C D1 D2] = getc(lat,lon)
% lat: [radians]
% lon: [radians]
% D1: dC/d(lat) [m/s per radian]
% D2: dC/d(lon) [m/s per radian]
% F is global griddedinterpolant F(lon,lat) in radians
% Esimate the derivatve dF/ds = F(s + ds)-F(s-ds) / 2ds

global F

% Get phase speed
h = F(lon,lat);
c = getcp(f,h); 

% Step Size (rad)
ds = 10e-5;   %About 50m

% Estimate phase speed gradients
h_ahead = F(lon,lat+ds);
h_behind = F(lon,lat-ds);
c_ahead = getcp(f,h_ahead);
c_behind = getcp(f,h_behind);
D1 = (c_ahead-c_behind)/(2*ds);

h_ahead = F(lon+ds,lat);
h_behind = F(lon-ds,lat);
c_ahead = getcp(f,h_ahead);
c_behind = getcp(f,h_behind);
D2 = (c_ahead-c_behind)/(2*ds);


end
