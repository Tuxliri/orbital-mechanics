function [rp,a,e,delta,vinfP] = flybyunpow(vinfM,mu,Delta,u)
%FLYBYUNPOW function that computes an unpowered flyby manoeuvre
% INPUT:
%   vinfM[3x1]  incoming velocity at infinity in the planetocentric 
%               reference frame    [km/s]
%   mu[1]       gravitational parameter of the planet           [km^3/s^2]
%   Delta[3x1]  impact parameter vector                         [km]
%   u[3x1]      unit axis of rotation vector around which the 
%               velocity is rotated     [-]

% OUTPUT:
%   rp[3x1]     perigee radius wrt planet                       [km]
%   a[1]        semi-major axis of the hyperbola                [km]
%   e[1]        eccentricity of the hyperbola                   [-]
%   delta[1]    turning angle of hyperbola                      [rad]
%   vinfP[3x1]  outgoing velocity at infinity in the planetocentric
%               reference frame                                 [km/s]
% 
% AUTHORS:
%   Davide Iafrate
% 
% REVISIONS:
%   04-12-2020  first revision

vinf = norm(vinfM);

% Calculate semi-major axis
a = -mu/vinf^2;

% Calculate the turning angle
delta = 2*atan2(-a,Delta);

% Calculate eccentricity
e = 1/sin(delta/2);

% Calculate magnitude of perigee radius
rp = a*(1-e);

% Rotation about the u axis using Rodrigues' formula

vinfP = vinfM*cos(delta) + sin(delta).*cross(u,v) + dot(u,vinfM).*u.*(1-cos(delta));


end

