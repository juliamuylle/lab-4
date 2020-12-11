clc
clear all
close all

%how to choose the range for rp0?

%Data
mu_S = astroConstants(4);       %[km^3/s^2] Planetary constant of Sun
mu_E = astroConstants(13);      %[km^3/s^2] Planetary constant of Earth
Re = astroConstants(23);        %[km] Radius of Earth

%Heliocentric velocity of the s/c before flyby
V_M = [31.5;4.69;0];            %[km/s]
%Heliocentric velocity of the s/c after flyby
V_P = [38.58;0;0];              %[km/s]

%Position of Earth
AU = astroConstants(2);             %[km] Astronomical unit
r_e=[0;-AU;0];                      %[km]
%Velocity of Earth
Vpl_norm=sqrt(mu_S/norm(r_e));            %[km/s]
Vpl=cross([0;0;1],r_e/norm(r_e)).*Vpl_norm; %[km/s] Velocity vector of Earth

%location of the incoming asymptote:
%u: vector around which the velocity vector gets rotated
%Earth is always assmed to be in the ecliptic


%Velocities relative to the planet
%entry velocity in SOI
vinfM = V_M-Vpl;        %[km/s]
%exit velocity in SOI
vinfP = V_P-Vpl;        %[km/s]

[deltav_perig,rp,delta,arcs] = flybyPow(vinfM,vinfP,mu_E,Re);

u = cross(vinfM,vinfP);
u = u/norm(u);          % normal to orbital plane


%%
%Total cost
deltaV_tot = V_P-V_M;

%Apse line negative direction
dir_rp = -deltaV_tot/norm(deltaV_tot);  %since magnitudes of vinfM and vinfP are different it is not exactly on the apse line

rp = dir_rp.*rp;
vp_dir = cross(u,dir_rp);

%Velocities at perigee
vpM = arcs(1,2);
vpM = vpM*vp_dir;
vpP = arcs(2,2);
vpP = vpP*vp_dir;

%Set options
options = odeset ( 'RelTol', 1e-13,'AbsTol', 1e-14 );
%State vector
y0M = [rp;vpM];
y0P = [rp;vpP];
%Set time span
tspanM = linspace(0,-4000,10000);
tspanP = linspace(0,4000,10000);
%Propagation
[t,YM] = ode113(@(t,y) odefun(mu_E,y,t),tspanM,y0M,options);
YM = YM/Re;
[t,YP] = ode113(@(t,y) odefun(mu_E,y,t),tspanP,y0P,options);
YP =YP/Re;

figure(1)
run Plot_Earth
hold on
plot3(YM(:,1),YM(:,2),YM(:,3))
plot3(YP(:,1),YP(:,2),YP(:,3))
xlabel('x[Re]')
ylabel('y[Re]')
grid on
axis equal
view(0,90)





