function dy = TwoBP_OdeFn(y, mu)

%INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[2x1] State of the orbit (position and velocity) [ R, R/T ]
% mu[1] Gravitational constant
% r[1] Orbit radius modulus
%
% OUTPUT:
% dy[2x1] Derivative of the state [ L/T^2, L/T^3 ]
r_v = [y(1) y(2) y(3)];
r_mod  = norm(r_v);

dy = [y(4)
      y(5)
      y(6)
      (-mu/r_mod^3)* y(1)
      (-mu/r_mod^3)* y(2)
      (-mu/r_mod^3)* y(3)];
end