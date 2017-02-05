function [ tfac ] = gettfac( fr,h )
% [ tfac ] = gettfac( fr,h )
%   Get shoaling coefficient from deep water for location

g=9.81;
k_local = getk(fr,h);
cg_local = 0.5*(g*tanh(k_local*h)+g*(k_local*h)*(sech(k_local*h)^2))/sqrt(g*k_local*tanh(k_local*h));            

% 5000 meters = deep
k_deep = getk(fr,5000);
cg_deep = 0.5*(g*tanh(k_deep*h)+g*(k_deep*h)*(sech(k_deep*h)^2))/sqrt(g*k_deep*tanh(k_deep*h));            

tfac = (k_local/k_deep)*(cg_deep/cg_local);


end

