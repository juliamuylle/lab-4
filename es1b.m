
clear
clc
close all


%% Data
vinfM = [15.1; 0; 0];              %Vinf before flyby, planet-centric frame, in km/s

AU      = astroConstants(2);       %Austronomical unit, in km
r_pl = [1; 0; 0]*AU;               %Position vector of the planet, in km

mu_earth = astroConstants(13);     %Planetary constant of Earth, in km^3/s^2
mu_sun   = astroConstants(4);      %Planetary constant of Sol, in km^3/s^2

Re = astroConstants(23);           %Earth's radius

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

%%Impact parameter variation
Delta_vec = [9200:1000:13200];

for i = 1:length(Delta_vec)
    Delta = Delta_vec(i);
    [rp_norm,delta,a,e,vinfP(:,i),h] = flybyUnpow_ok(vinfM,mu_earth,Delta,u);
    
    delta_v = vinfP(:,i)-vinfM;
    delta_v_norm = norm(delta_v);
    
    rdir = -delta_v/delta_v_norm;
    rp = rp_norm*rdir;

    vp_dir = cross(u,rdir);
    h = sqrt(mu_earth*rp_norm*(1+e));
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



    figure(1)
    hold on
    plot3(YM(:,1),YM(:,2),YM(:,3))
    hold on
    plot3(YP(:,1),YP(:,2),YP(:,3))
    view([0,90])
    grid on
    axis equal

       
end

figure(1)
hold on
[x,y,z] = sphere(50);
surf(x,y,z)
x=x*1;
y=y*1;
z=z*1;
surf(x,y,z,'EdgeColor','none','FaceAlpha',0.1);
map=imread('EarthTexture.jpg');
map=imrotate(map,180);
warp(x,y,z,map);
view([0,-90])

%% Heliocentric orbit propagation

%Initial conditions
sc_v0 = vinfM;                                  %Transformation of v0 to sun-centric frame
Vpl_norm = sqrt(mu_sun/norm(r_pl));          %velocity of the planet, in Km/s
V_pl      = Vpl_norm.* cross([0;0;1],r_pl)/norm(r_pl);

minus_v0 = sc_v0 + V_pl;
minus_y0 = [r_pl; minus_v0];

plus1_v0 = vinfP(:,1)+V_pl;
plus2_v0 = vinfP(:,2)+V_pl;
plus3_v0 = vinfP(:,3)+V_pl;
plus4_v0 = vinfP(:,4)+V_pl;
plus5_v0 = vinfP(:,5)+V_pl;

plus1_y0 = [r_pl; plus1_v0];
plus2_y0 = [r_pl; plus2_v0];
plus3_y0 = [r_pl; plus3_v0];
plus4_y0 = [r_pl; plus4_v0];
plus5_y0 = [r_pl; plus5_v0];

%Options for the ODE solver
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

%Timespans
minus_tspan = linspace(0, -20000000, 10000); %Recordatorio: el tercer termino es el nro de segmentos deseados
plus_tspan = linspace(0, 20000000, 10000);

%Propagation
[minus_t, minus_y] = ode113(@(t,y) TwoBP_OdeFn(y,mu_sun), minus_tspan, minus_y0, options);
minus_y = minus_y.*(1/AU);

[plus1_t, plus1_y] = ode113(@(t,y) TwoBP_OdeFn(y,mu_sun), plus_tspan, plus1_y0, options);
[plus2_t, plus2_y] = ode113(@(t,y) TwoBP_OdeFn(y,mu_sun), plus_tspan, plus2_y0, options);
[plus3_t, plus3_y] = ode113(@(t,y) TwoBP_OdeFn(y,mu_sun), plus_tspan, plus3_y0, options);
[plus4_t, plus4_y] = ode113(@(t,y) TwoBP_OdeFn(y,mu_sun), plus_tspan, plus4_y0, options);
[plus5_t, plus5_y] = ode113(@(t,y) TwoBP_OdeFn(y,mu_sun), plus_tspan, plus5_y0, options);
plus1_y = plus1_y.*(1/AU);
plus2_y = plus2_y.*(1/AU);
plus3_y = plus3_y.*(1/AU);
plus4_y = plus4_y.*(1/AU);
plus5_y = plus5_y.*(1/AU);

%% Plotting


%plot(hy_impact, hy_delta,'LineWidth',2);
%yyaxis right
%hold on

figure(2)
Orb_minus = plot3(minus_y(:,1), minus_y(:,2), minus_y(:,3),'b','LineWidth',2);
hold on
Orb_plus1 = plot3(plus1_y(:,1), plus1_y(:,2), plus1_y(:,3),'r','LineWidth',2);
hold on
Orb_plus2 = plot3(plus2_y(:,1), plus2_y(:,2), plus2_y(:,3),'y','LineWidth',2);
hold on
Orb_plus3 = plot3(plus3_y(:,1), plus3_y(:,2), plus3_y(:,3),'g','LineWidth',2);
hold on
Orb_plus4 = plot3(plus4_y(:,1), plus4_y(:,2), plus4_y(:,3),'c','LineWidth',2);
hold on
Orb_plus5 = plot3(plus5_y(:,1), plus5_y(:,2), plus5_y(:,3),'m','LineWidth',2);
hold on
[x,y,z] = sphere(10);
sun_r = 0.1;
sun = surf(x*sun_r, y*sun_r, z*sun_r);
set(sun,'facecolor',[1 1 0]);
view(0,90)
legend('Before Fly-by', 'After Fly-by w/per@2000km', 'After Fly-by w/per@4000km', 'After Fly-by w/per@6000km', 'After Fly-by w/per@8000km', 'After Fly-by w/per@10000km');
%axis equal
grid on
title('Orbit')
xlabel('X [AU]')
ylabel('Y [AU]')

