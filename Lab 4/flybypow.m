function [dv_p,delta,a,e,Delta,rp,vp0] = flybypow(vinfM,vinfP,mu,rp_min)
%FLYBYUNPOW function that computes an unpowered flyby manoeuvre
% 
% PROTOTYPE [dv_p,delta,Delta] = flybypow(vinfM,vinfP,mu,rp_min)
% 
% INPUT:
%   vinfM[3x1]  incoming velocity at infinity in the planetocentric 
%               reference frame    [km/s]
%   vinfP[3x1]  outgoing velocity at infinity in the planetocentric 
%               reference frame    [km/s]
%   mu[1]       gravitational parameter of the planet           [km^3/s^2]
%   rp_min[1]   minimum perigee radius in magnitude             [km]
% 
% OUTPUT:
%   dv_p[1]     cost of the perigee burn maneuvre               [km/s]
%   delta[1]    total turning angle                             [rad]
%   a[2]        semi-major axis before and after perigee burn   [km]
%   e[2]        eccentricities before and after perigee burn    [-]
%   Delta[1]    arrival impact parameter                        [km]
%   vp0[1]      velocity at perigee passage                     [km/s]
% 
% AUTHORS:
%   Davide Iafrate
% 
% REVISIONS:
%   04-12-2020  first revision

vinfm = norm(vinfM);
vinfp = norm(vinfP);

% Calculate the turning angle of the maneuvre
delta = acos(dot(vinfM,vinfP)/(vinfm*vinfp));

% Solve the implicit equation for x
eq =  @(x) asin(1/(1 + x*vinfm^2/mu)) + asin(1/(1 + x*vinfp^2/mu)) - delta;

rp = fzero(eq,6378);

% Check validity of the solution, that is if it is larger than the minimum
% perigee radius
if rp < rp_min
    rp = NaN;
    
end

% Calculate the arrival semi-major axis
am = - mu/vinfm^2;
ap = -mu/vinfp^2;

a = [am ap];

% Calculate eccentricities
em = 1 + rp*(vinfm^2)/mu;
ep = 1 + rp*(vinfp^2)/mu;

e = [em ep];

% Calculate the impact parameter

Delta = - am/tan(delta/2);

% Compute the velocities of the two hyperbolic arcs at perigee
vpm = sqrt(mu*(1/am + 2/rp));
vpp = sqrt(mu*(1/ap + 2/rp));

vp0=vpm;

% Compute the delta v required at perigee burn
dv_p = vpp - vpm;