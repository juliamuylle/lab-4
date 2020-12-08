
clear
clc
close all

    
%% Data
vinfM = [15.1; 0; 0];              %Vinf before flyby, planet-centric frame, in km/s
Delta = 9200;                      %Impact parameter modulus, in km

AU      = astroConstants(2);       %Austronomical unit, in km
r_pl = [1; 0; 0]*AU;               %Position vector of the planet, in km

mu_earth = astroConstants(13);     %Planetary constant of Earth, in km^3/s^2
mu_sun   = astroConstants(4);      %Planetary constant of Sol, in km^3/s^2

Re = astroConstants(23);           %Earth's radius


%% Computing vinf (Vinf rotation)

%Determination of side of the flyby

flyby_location = 'front';
switch flyby_location
case 'front'
    u = [0;0;-1];
    case 'behind'
        u = [0;0;1];
    case 'under'
        u = cross ([0;0;-1],vinfM/norm(vinfM));
    case 'over'
        u = cross([0;0;1],vinfM/norm(vinfM));
end

[rp_norm,delta,a,e,vinfP,h] = flybyUnpow_ok(vinfM,mu_earth,Delta,u);

delta_v      = vinfP - vinfM;
delta_v_norm = norm(delta_v);



%% Heliocentric orbit propagation

%Initial conditions
sc_v0    = vinfM;                                  %Transformation of v0 to sun-centric frame
Vpl_norm = sqrt(mu_sun/norm(r_pl));          %velocity of the planet, in Km/s
Vpl      = Vpl_norm.* cross([0;0;1],r_pl)/norm(r_pl);

V_M = sc_v0 + Vpl;
V_P = vinfP + Vpl;

minus_y0 = [r_pl; V_M];
plus_y0  = [r_pl; V_P];

sc_r0    = [1e6; 0; 0];
FlyBy_y0 = [sc_r0; sc_v0];

%Options for the ODE solver
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

%Timespans
minus_tspan = linspace(0, -20000000, 10000); 
plus_tspan = linspace(0, 20000000, 10000);
FlyBy_tspan = linspace(0, 100000, 1000);

%Propagation
[minus_t, minus_y] = ode113(@(t,y) TwoBP_OdeFn(y,mu_sun), minus_tspan, minus_y0, options);
minus_y = minus_y.*(1/AU);
[plus_t, plus_y] = ode113(@(t,y) TwoBP_OdeFn(y,mu_sun), plus_tspan, plus_y0, options);
plus_y = plus_y.*(1/AU);

[FlyBy_t, FlyBy_y] = ode113(@(t,y) TwoBP_OdeFn(y,mu_earth), FlyBy_tspan, FlyBy_y0, options);

%% Planetocentric orbit propagation

rdir = -delta_v/delta_v_norm;
rp = rp_norm*rdir;

vp_dir = cross(u,rdir);     %Direction of the velocity at pericenter
vp_magn = mu_earth/h*(1+e);
vp = vp_magn*vp_dir;

y0 = [rp;vp];
options = odeset ( 'RelTol', 1e-13,'AbsTol', 1e-14 );

tspan = linspace(0,2000,10000);
[t,YM] = ode113(@(t,y) odefun(mu_earth,y,t),tspan,y0,options);
YM = YM/Re;


tspan = linspace(0,-2000,10000);
[t,YP] = ode113(@(t,y) odefun(mu_earth,y,t),tspan,y0,options);
YP =YP/Re;





%% Plotting

figure(1)

[x,y,z] = sphere(50);
surf(x,y,z)
x=x*1;
y=y*1;
z=z*1;
surf(x,y,z,'EdgeColor','none','FaceAlpha',0.1);
map=imread('EarthTexture.jpg');
map=imrotate(map,180);
warp(x,y,z,map);
hold on
plot3(YM(:,1),YM(:,2),YM(:,3))
hold on
plot3(YP(:,1),YP(:,2),YP(:,3))
view([0,-90])
grid on
axis equal


figure(2)
Orb_minus = plot3(minus_y(:,1), minus_y(:,2), minus_y(:,3),'b','LineWidth',2);
hold on
Orb_plus = plot3(plus_y(:,1), plus_y(:,2), plus_y(:,3),'r','LineWidth',2);
hold on
run Plot_Sun
legend('Before Fly-by', 'After Fly-by');
%axis equal
grid on
title('Orbit')
xlabel('X [AU]')
ylabel('Y [AU]')


