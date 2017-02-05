function [c, D1, D2] = getc(lat,lon)
% function [C D1 D2] = getc(lat,lon)
% lat: [radians]
% lon: [radians]
% fr: frequency [Hz]  
% C: sound speed [m/s]
% D1: dC/d(lat) [m/s per radian]
% D2: dC/d(lon) [m/s per radian]
% F is global griddedinterpolant F(lon,lat) in radians
% Esimate the derivatve dF/ds = F(s + ds)-F(s-ds) / 2ds

global F

% Get phase speed
c = F(lon,lat);

% Step Size (rad)
ds = 6e-5;   %About 10m

% dF/dlat
D1 = (F(lon,lat+ds)-F(lon,lat-ds))/(2*ds);

% dF/dlon
D2 = (F(lon+ds,lat)-F(lon-ds,lat))/(2*ds);



end
