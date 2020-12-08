function dy = odefun (t,y,mu)

%odefun gives the derivative of the state vector using the dynamic equation of 2BP
%
%PROTOTYPE: 
%     dy = odefun (t,y,mu)
% 
% INPUT:
%     t           Time [s]
%     y [6x1]     State Vector containing position and velocity vectors concatanate
%     mu [1]      Gravitational constant of the Earth [km^3/s^2]
%     
% OUTPUT:
%     dy [6x1]     Derivative fo the state vector
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


r  = y(1:3);
dy = [y(4:6); -mu/(norm(r)^3)*y(1:3)];

end