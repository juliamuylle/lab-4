function [a] = rep_groundtrack(k,m,we,mu)

% rep_groundtrack computes the value of the semi-major axis to obtain repeating groundtrack

%PROTOTYPE: 
%     [a,T] = rep_groundtrack(k,m,we,mu)
%
% INPUT:
%     k [1]       Revolution of the satelite
%     m [1]       Rotations of the planet
%     we [1]      Angular velocity of the Earth [rad/s]
%     mu [1]      Gravitational constant of the Earth [km^3/s^2]
%     
% OUTPUT:
%    a [1]        Semimajor axis [km]
%    T [1]        Period [s]
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

%Earth’s rotational period 
Te = 2*pi/we;

%s/c’s orbital period
T  = Te*m/k;

%Semi-major axis for repeating groundtrack
a  = ((T/(2*pi))^2*mu)^(1/3);

