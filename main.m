clc
clear all
close all

        
%Keplerian parameters  
%kepl_coord(1) = a		semi-major axis	[km]
%kepl_coord(2) = e		eccentricity		[-]
%kepl_coord(3) = i		inclination		[rad]
%kepl_coord(4) = omega	RAAN			[rad]
%kepl_coord(5) = w		argument of perigee	[rad]
%kepl_coord(6) = f 		true anomaly		[rad]


kepl_coord = [34814, 0.5054,deg2rad(51.6177),deg2rad(156.0701),deg2rad(253.1633),0];

mu         = astroConstants(13);    %Planetary constant [km^3/s^2]
Re         = astroConstants(23);    %Earth’s radius [km]
we         = deg2rad(15.04/3600);   %Angular speed of Earth [rad/s]
long_green = 0;                     %Greenwich longitude [rad]
J2         = 0.00108263;            %J2 perturbation
%n_orbit    = 10;

%Orbital period
T = 2*pi*sqrt(kepl_coord(1)^3/mu);

%Different cases
n = input('enter\n 1 for one orbit period\n 2 for a day\n 3 for 10 days\n')
switch n
case 1
    tmax = T;
    
case 2
    tmax = 86400;
    
case 3
    tmax = 8640000;

end 

%Set time span
n_orbit = tmax/T;
t0 = 0;
t  = linspace(t0,tmax,100000);

%Compute groundtrack
[alpha,delta,long,lat,Y] = Groundtracks(kepl_coord,long_green,t,mu,we,t0);


%Repeating groundtrack
m = 3; %rot of planet
k = 4; % we have to choose them. %rev of sc

[a_new] = rep_groundtrack(k,m,we,mu);
kepl_coord_rep = kepl_coord;
kepl_coord_rep(1) = a_new;
T_new = 2*pi*sqrt(a_new^3/mu);
t_rep = linspace(0,n_orbit*T_new,100000); 
[alpha_new,delta_new,long_new,lat_new] = Groundtracks(kepl_coord_rep,long_green,t_rep,mu,we,t0);

%Groundtrack with J2 perturbation
[alpha_J2,delta_J2,long_J2,lat_J2,YJ2] = Groundtracks_J2(kepl_coord,long_green,t,mu,we,t0,J2);

%Repeating groundtrack with J2 perturbation
a_J2_new = rep_groundtrack_J2(k,m,we,mu,J2,Re,kepl_coord(1),kepl_coord(2),kepl_coord(3));
kepl_coord_rep_J2 = kepl_coord;
kepl_coord_rep_J2(1) = a_J2_new;
T_J2_new = 2*pi*sqrt(a_J2_new^3/mu);
t_rep_J2 = linspace(0,n_orbit*T_J2_new,100000); %dubbio: non so se usare questo o sempre lo stesso
[alpha_J2_new,delta_J2_new,long_J2_new,lat_J2_new] = Groundtracks_J2(kepl_coord_rep_J2,long_green,t_rep_J2,mu,we,t0,J2);


% Some conditions to avoid having horizontal lines connecting opposite
% longitude

for k = 2:length(long)
    if long(k)*long(k-1)<0
        long(k) = NaN;
    end
    if long_new(k)*long_new(k-1)<0;
        long_new(k) = NaN;
    end
    if long_J2(k)*long_J2(k-1)<0
        long_J2(k) = NaN;
    end
    if long_J2_new(k)*long_J2_new(k-1)<0
        long_J2_new(k) = NaN;
    end
end

%Plot

%Plot of the orbit and of the Earth
figure(1)
plot3(Y(:,1),Y(:,2),Y(:,3),'b')
xlabel('rx[km]')
ylabel('ry[km]')
zlabel('rz[km]')
title('Orbit')
hold on
grid on
run Plot_Earth

%Plot of the orbit and of the Earth with J2 effect
figure(2)
plot3(YJ2(:,1),YJ2(:,2),YJ2(:,3),'b')
xlabel('rx[km]')
ylabel('ry[km]')
zlabel('rz[km]')
title('Orbit with J2 effect')
hold on
grid on
run Plot_Earth


figure(3)
c = imread('EarthTexture.jpg');
c1 = imrotate(c,-180);
image([180,-180],[-90,+90],c1)
ax = gca;
ax.YDir = 'normal';
hold on
plot(long,lat,'g')
hold on 
plot(long(1),lat(1),'go','Linewidth',2)
hold on
plot(long(end),lat(end),'gs','Linewidth',2)
legend('Orbit','Start','End')
title('Orbit without perturbation')
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')



figure(4)
c = imread('EarthTexture.jpg');
c1 = imrotate(c,-180);
image([180,-180],[-90,+90],c1)
ax = gca;
ax.YDir = 'normal';
hold on
plot(long,lat,'g')
hold on 
plot(long(1),lat(1),'go','Linewidth',2)
hold on
plot(long(end),lat(end),'gs','Linewidth',2)
hold on
plot(long_J2,lat_J2,'y')
hold on
plot(long_J2(1),lat_J2(1),'yo','Linewidth',2)
hold on
plot(long_J2(end),lat_J2(end),'ys','Linewidth',2)
legend('Orbit','Start','End')
title('Orbit with J2 perturbation')
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')



figure(5)
c = imread('EarthTexture.jpg');
c1 = imrotate(c,-180);
image([180,-180],[-90,+90],c1)
ax = gca;
ax.YDir = 'normal';
hold on
plot(long,lat,'g')
hold on 
plot(long(1),lat(1),'go','Linewidth',2)
hold on
plot(long(end),lat(end),'gs','Linewidth',2)
hold on
plot(long_new,lat_new,'r')
hold on
plot(long_new(1),lat_new(1),'ro','Linewidth',2)
hold on
plot(long_new(end),lat_new(end),'rs','Linewidth',2)
legend('Orbit','Repeating orbit','Start','End','Start rep','End rep')
title('Repeating groundtrack')
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')

figure(6)
c = imread('EarthTexture.jpg');
c1 = imrotate(c,-180);
image([180,-180],[-90,+90],c1)
ax = gca;
ax.YDir = 'normal';
hold on
plot(long,lat,'g')
hold on 
plot(long(1),lat(1),'go','Linewidth',2)
hold on
plot(long(end),lat(end),'gs','Linewidth',2)
hold on
plot(long_new,lat_new,'r')
hold on 
plot(long_new(1),lat_new(1),'ro','Linewidth',2)
hold on
plot(long_new(end),lat_new(end),'rs','Linewidth',2)
hold on
plot(long_J2_new,lat_J2_new,'y')
hold on
plot(long_J2_new(1),lat_J2_new(1),'yo','Linewidth',2)
hold on
plot(long_J2_new(end),lat_J2_new(end),'ys','Linewidth',2)
legend('Original orbit','Start','End','Repeating Unperturbed','Start rep','End rep','Repeating perturbed','Start rep pert','End rep pert')
title('Repeating groundtrack')
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')


