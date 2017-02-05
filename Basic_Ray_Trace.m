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

% Load bathymetry
Bathy = load('bathygrid_500m.mat');

% Starting location, lat North, lon East
lat_start   = 33.60;
lon_start   =-118.30;

% wave energy freq
fr = 0.075;

%%%%%% Create escape bounds %%%%%%%%
lat_ulim    = 34.80*myconst.d2r;
lat_dlim    = 32.5*myconst.d2r;
lon_lim     = -121*myconst.d2r;

% Use slanting line for western edge
lon_ax      = [-121 -119];
lat_ax      = 34.5-(1.3)*(lon_ax+121);
slantborder = @(lat,lon) (lat-(34.5-(1.3)*(lon+121)));

%%%%%%%%% Created griddedInterpolant %%%%%%%%%
% This is global variable so that ray funcs can use it
% Linear is fastest
Bathy.h1_complete = Bathy.h1;
Bathy.h1(Bathy.h1>0) = NaN;
F = griddedInterpolant({Bathy.lon*myconst.d2r,Bathy.lat*myconst.d2r},Bathy.h1');
F.ExtrapolationMethod = 'linear';
F.Method = 'linear';


%%


%%%%%%%%% Setup for Ray Tracing %%%%%%%%%%%%%

%Ray Tracing parameters
htry = 1*60;            % Initial step size
eps = 3e-4;             % Step error threshold (3e-5 is accurate)
yscale = [1 1 1];       % Scale errors if needed, 
minStep = 15;            % Minimum step (seconds), don't want to take all day
maxStep = 2*60;           % Maximum step (seconds), don't want to jump islands/shoals

%Intialize counters
stepcount = 0;          % Ray location
raycount = 0;           % Number of rays traced

%Initial ray angles to sweep through
alpha_start = 160:1:300;            %Starting ray angle
alpha_end = NaN(size(alpha_start)); %Ending ray angle

%Initialize space for rays and time
ray_time = NaN(1,length(alpha_start));  
ray_loc = NaN(2,200000);                % Rays(lat,lon), NAN separates rays


% Loop over initial directions
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
    ray_loc(1,stepcount) = P(2)*myconst.r2d; %Record initial position
    ray_loc(2,stepcount) = P(3)*myconst.r2d;
    
    % Continue tracing until hit land (NAN) or out of bounds
    while all(~isnan(P)) && P(2) > lat_dlim && P(2) < lat_ulim && P(3) > lon_lim 
        stepcount = stepcount+1;
        
        dPdt = getderivs(P,fr);
        [ P, t, hnext, maxcount ] = rkqc(P,dPdt,fr,t,hnext,eps,yscale,minStep);
        
        hnext = min(hnext,maxStep);
        hnext = max(hnext,minStep);
        
        plot(P(3)*myconst.r2d,P(2)*myconst.r2d,'r.')
        ray_loc(1,stepcount) = P(2)*myconst.r2d;
        ray_loc(2,stepcount) = P(3)*myconst.r2d;
        
        %plot(ray_loc(2,stepcount),ray_loc(1,stepcount),'.b')
        %pause(1e-4)
    end
    
    
    % If any part of P is NAN, set all to NAN
    if isnan(sum(P))
        P(1:3)=NaN;
    end
    
    % Save Last ray point
    ray_loc(1,stepcount) = P(2)*myconst.r2d;
    ray_loc(2,stepcount) = P(3)*myconst.r2d;
        
    % Mark Ending of Rays 
    stepcount = stepcount+1;
    ray_loc(1,stepcount) = -999;
    ray_loc(2,stepcount) = -999;
    
    % Mark ending angle or blocked angle (NaN)
    alpha_end(rr) = P(1)*myconst.r2d;
        
    % Save time integrated
    ray_time(1,raycount) = t;
    
end

%%
% Plot map

%%
clf
contourf(Bathy.lon,Bathy.lat,Bathy.h1_complete,[0 0])
hold on
contour(Bathy.lon,Bathy.lat,Bathy.h1,[-20 -50 -100 -250 -500 -1000],'Color',[.7 .7 .7])
shading flat
colormap(copper)

%%
plot(ray_loc(2,:),ray_loc(1,:),'.b')
xlim([-121 -117])
ylim([32.5 34.6])


return
%%

% Plot relation between starting and ending angles
subplot(2,1,1)
plot(alpha_end,alpha_start)



% % Integrate Gamma function
% K2 = zeros(360,1);
% theta_edge = 0:5:360;
% theta_cen = 2.5:5:357.5;
% dGam = diff(alpha_start); %slope
% CGam = conv(Gam(2,:),ones(1,2)/2,'valid'); %midpt
%     for ii=1:length(theta_cen)
%         inds = CGam >= theta_edge(ii) & CGam < theta_edge(ii+1);
%         K2(ii) = sum(dGam(inds))*gettfac(fr(ff),bdepth);   %Include depth-shoaling factor tfact 
%         Td(ii,ff,1) = sum(dGam(inds).*cosd(theta_cen(ii)))*gettfac(fr(ff),bdepth);  
%         Td(ii,ff,2) = sum(dGam(inds).*sind(theta_cen(ii)))*gettfac(fr(ff),bdepth); 
%         Td(ii,ff,3) = sum(dGam(inds).*cosd(2*theta_cen(ii)))*gettfac(fr(ff),bdepth);  
%         Td(ii,ff,4) = sum(dGam(inds).*sind(2*theta_cen(ii)))*gettfac(fr(ff),bdepth);  
%     end
%     
%     Ta0_all(:,ff) = K2;


