function [rp,delta,a,e,vinfP,h] = flybyUnpow_ok(vinfM,mu_E,Delta,u)

v_inf = norm(vinfM);

% Solving the hyperbola
a = -mu_E/v_inf^2;                  %semi-mayor axis, in km
delta = 2*atan(-a/Delta);           %turning angle, in rad
deltaDeg = rad2deg(delta);         %turning angle, in deg
e = 1/(sin(delta/2));               %eccentricity
rp = a*(1-e);                       %radius of pericentre, in km
h = sqrt(mu_E*rp*(1+e));            %angular momentum

%rotation of Vinf to obtain Vinf at SOI exit
R = axang2rotm([u',delta]);
vinfP = R*vinfM;
