function [alpha,delta,lon,lat,Y] = Groundtracks(kep_coord,long_green,t,mu,we,t0)
%Groundtracks computes RAAN, declination, latitude and longitude to have the projection of 
%the orbit of s/c onto Earth’s surface

%PROTOTYPE: 
%     [alpha,delta,lon,lat] = Groundtracks(kep_coord,long_green,t,mu,we,t0)
% 
% INPUT:
%     kep_coord [6x1]     keplerian coordinates = [a,e,i,omega,w,f], [km,~,rad,rad,rad,rad]
%     long_green [1]      longitude of Greenwich meridian [rad]
%     t                   Time [s]
%     mu [1]              Gravitational constant of the Earth [km^3/s^2]
%     we [1]              Angular velocity of the Earth [rad/s]
%     t0 [1]              Initial time [s]
%     
% OUTPUT:
%     alpha               Right Ascension [deg]     
%     delta               Declination [deg]
%     lon                 Longitude [deg]
%     lat                 Latitude [deg]
%
% CONTRIBUTORS
%       Bertolini Edoardo
%       Busi Silvia
%       Muylle Julia
%       Pellegrini Matias
%
% VERSIONS
%
% 30/11/2020: First Version


[r0,v0] = kep2car(kep_coord,mu);
y0      = [r0,v0];
options = odeset('RelTol',1e-13,'AbsTol',1e-14);
[~,Y]   = ode113(@(t,y) odefun(t,y,mu),t,y0,options);
Y       = [Y(:,1) Y(:,2) Y(:,3)];

%Right ascension of the ascending node and declination
for i =1:length(Y(:,1))
    
    delta(i) = asin(Y(i,3)/norm(Y(i,:)));
    alpha(i) = acos(Y(i,1)/(norm(Y(i,:))*cos(delta(i))));
    
    if (Y(i,2)/norm(Y(i,:))) <= 0
        alpha(i) = 2*pi-alpha(i);
    end
    
    %Longitude of Greenwich meridian
    thetag(i) = long_green+we*(t(i)-t0);
    
    %Longitude
    lon(i)    = alpha(i)-thetag(i);
    lon(i)    = wrapToPi(lon(i));
    
    %Latitude
    lat(i)    = delta(i);
    
end

lat = rad2deg(lat);
lon = rad2deg(lon);

