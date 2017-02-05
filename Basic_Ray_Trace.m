% Backwards Ray Tracer 
% Based on refraction equations for ellipsoid,
% See Eq (2) in Munk, et al (1988) JPO
% Crosby, S. C.

clearvars

% Global variables
global F % GriddedInterpolant of phase speed

% Constants Used in various functions
global ellipseCon
global myconst
ellipseCon.req  = 6378137.0;
ellipseCon.finv = 298.257223563;
ellipseCon.eps2 = 2/ellipseCon.finv;
myconst.d2r     = pi/180;
myconst.r2d     = 180/pi;


% Starting location, lat North, lon East
lat_start = 34.229;
lon_start = -119.796;

%%%%%% Create escape bounds %%%%%%%%
lat_ulim = 34.80*myconst.d2r;
lat_dlim = 32.00*myconst.d2r;
lon_lim = -121*myconst.d2r;

% Use slanting line for western edge
lon_ax = [-121 -119];
lat_ax = 34.5-(1.3)*(lon_ax+121);
slantborder = @(lat,lon) (lat-(34.5-(1.3)*(lon+121)));

% Load plottable map
bathy = load('/Users/sccrosby/Documents/CDIP/Transforms/elev_2/bathygrid_500m.mat');
clf
pcolor(bathy.lon,bathy.lat,bathy.h1)
shading flat
caxis([-500 0])
colorbar
hold on
contour(bathy.lon,bathy.lat,bathy.h1,[0 0],'k')

% Choose wave frequency, which loads appropriate phase speed map
fr = 0.04;
fname = sprintf('/Users/sccrosby/Documents/CDIP/Transforms/PhaseMaps/PhaseInMLLW/phaseSpeed_F%4.4fHz.mat',fr);
load(fname);

% Ultimately we need a grid of phase speed with Lon/Lat
[LON, LAT] = meshgrid(lonV,latV);
lon = LON(1,:);
lat = LAT(:,1);



% Load phase speed map
% cp, lonV, latV

%%%%%%%%% Created griddedInterpolant %%%%%%%%%
% This is global variable so that ray funcs can use it
% Spline doesn't seem to work with NANs
% Use linear or Cubic
% Linear is faster
F = griddedInterpolant({lon*myconst.d2r,lat*myconst.d2r},cp');
F.ExtrapolationMethod = 'linear';
F.Method = 'linear';

%%%%%%%%% Setup for Ray Tracing %%%%%%%%%%%%%

%Ray Tracing parameters
htry = 1*60;            % First step size to try
eps = 10e-5;           % 3e-5 seems to very accurate, 100e-5 is reasonable
Yscale = [1 1 1];       % Scale errors if needed, alpha error is dominate!
minStep = 5;            % Minimum step (seconds), don't want to take all day
maxStep = 60;         % Maximum step (seconds), don't want to jump islands/shoals

%Intialize everything for loop
stepcount = 1;      % Every ray step is counted, cummulative (to save ray paths)
raycount = 0;       % Counts the number of rays traced (length of Gam fxn)
alpha_start = 120:1:360;           %Initial sweep
alpha_end = NaN(size(alpha_start)); %Initial sweep

%Set aside space to save Gam func, rays, and time
raytime = NaN(1,length(alpha_start));      % Amount of seconds integrated over
Rays = NaN(2,200000);       % Rays(lat,lon), NAN separates rays



for rr = 1:length(alpha_start); % Begin loop of initial rays range
    raycount = raycount+1;
    stepcount = stepcount+1;
    
    % Set first Position Vector, P [ alpha lat lon ], radians
    P0 = [alpha_start(rr) lat_start lon_start]*myconst.d2r; %Position and angle (radians)
    t0 = 0;
    
    % Initialize & begin tracing
    P = P0;
    t = t0;
    hnext = htry;
    Rays(1,stepcount) = P(2)*myconst.r2d; %Record initial position
    Rays(2,stepcount) = P(3)*myconst.r2d;
    
    % Continue tracing until hit land (NAN) or out of bounds
    while ~isnan(P(1)) && P(2) > lat_dlim && P(2) < lat_ulim && P(3) > lon_lim && slantborder(P(2)*myconst.r2d,P(3)*myconst.r2d) > 0
        stepcount = stepcount+1;
        
        dPdt = getderivs(P);
        [ P, t, hnext, maxcount ] = rkqc(P,dPdt,t,hnext,eps,Yscale,minStep);
        
        hnext = min(hnext,maxStep);
        hnext = max(hnext,minStep);
        
        plot(P(3)*myconst.r2d,P(2)*myconst.r2d,'r.')
        Rays(1,stepcount) = P(2)*myconst.r2d;
        Rays(2,stepcount) = P(3)*myconst.r2d;
    end
    
    
    % If any part of P is nan, set all to NAN
    if isnan(sum(P))
        P(1:3)=NaN;
    end
    
    % Save Last ray point
    Rays(1,stepcount) = P(2)*myconst.r2d;
    Rays(2,stepcount) = P(3)*myconst.r2d;
    
    
    % Mark Ending of Rays (for later examination)
    stepcount = stepcount+1;
    Rays(1,stepcount) = -999;
    Rays(2,stepcount) = -999;
    
    % Mark ending angle or blocked angle (NaN)
    alpha_end(rr) = P(1)*myconst.r2d;
    %disp([P0(1)*myconst.r2d P(1)*myconst.r2d])

    
    % Save time integrated
    raytime(1,raycount) = t;
    
end



%%

