clc
clear all
close all

%how to choose the range for rp0?


V_M = [31.5;4.69;0];
V_P = [38.58;0;0];
AU = astroConstants(2);
r = [0;-1;0]*AU;
mu_S = astroConstants(4);
mu_E = astroConstants(13);
Re = astroConstants(23);

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

Vpl_norm = sqrt(mu_S/norm(r)); %take it from Eph in assignment
Vpl = cross([0;0;1],r/norm(r)).*Vpl_norm;

vinfM = V_M-Vpl;
vinfP = V_P-Vpl;

[deltav_perig,rp,delta,arcs] = flybyPow(vinfM,vinfP,mu_E,Re);



%%
deltaV_tot = V_P-V_M;
dir_rp = -deltaV_tot/norm(deltaV_tot);
rp = dir_rp.*rp;
vp_dir = cross(u,dir_rp);

vpM = arcs(1,2);
vpM = vpM*vp_dir;
vpP = arcs(2,2);
vpP = vpP*vp_dir;

options = odeset ( 'RelTol', 1e-13,'AbsTol', 1e-14 );

y0M = [rp;vpM];
y0P = [rp;vpP];

tspan = linspace(0,-2000,10000);
[t,YM] = ode113(@(t,y) odefun(mu_E,y,t),tspan,y0M,options);
YM = YM/Re;


tspan = linspace(0,2000,10000);
[t,YP] = ode113(@(t,y) odefun(mu_E,y,t),tspan,y0P,options);
YP =YP/Re;

figure(1)
% [x,y,z] = sphere(50);
% surf(x,y,z)
% x=x*1;
% y=y*1;
% z=z*1;
% surf(x,y,z,'EdgeColor','none','FaceAlpha',0.1);
% map=imread('earth.jpeg');
% map=imrotate(map,180);
% warp(x,y,z,map);
% hold on
plot3(YM(:,1),YM(:,2),YM(:,3))
hold on
plot3(YP(:,1),YP(:,2),YP(:,3))
view([0,90])
grid on





