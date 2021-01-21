function [dv_p,delta,a,e,Delta,rp,vp0] = flybypow(vinfM,vinfP,mu)
%FLYBYUNPOW function that computes an unpowered flyby manoeuvre
%
% PROTOTYPE:
%   [dv_p, delta, Delta] = flybypow(vinfM, vinfP, mu, rp_min)
%
% INPUT:
%   vinfM [ 3 x 1 ]  incoming velocity at infinity in the planetocentric
%                    reference frame                           [ km/s ]
%   vinfP [ 3 x 1 ]  outgoing velocity at infinity in the planetocentric
%                    reference frame                           [ km/s ]
%   mu    [ 1 x 1 ]  gravitational parameter of the planet     [ km^3/s^2 ]
%
% OUTPUT:
%   dv_p  [ 1 x 1 ]  cost of the perigee burn maneuvre         [ km/s ]
%   delta [ 1 x 1 ]  total turning angle                       [ rad ]
%   a     [ 2 x 1 ]  semi-major axis before and after perigee burn [ km ]
%   e     [ 2 x 1 ]  eccentricities before and after perigee burn  [ - ]
%   Delta [ 1 x 1 ]  arrival impact parameter                  [ km ]
%   vp0   [ 1 x 1 ]  velocity at perigee passage               [ km/s ]
%
% CONTRIBUTORS:
%   Alkady Marwan
%   Bossi Núñez Pedro
%   Demartini Davide
%   Iafrate Davide
%
% VERSIONS
%   2020-12-04: First version

vinfm = norm(vinfM);
vinfp = norm(vinfP);

% Calculate the turning angle of the maneuvre
delta = acos(dot(vinfM,vinfP)/(vinfm*vinfp));

% FIX: USE FLYBY-UNPOW?
if delta < 1e-6
    dv_p = 0;
    delta = 0;
    a = 0;
    e = 0;
    Delta = 0;
    rp = 0;
    vp0 = 0;
else
    % Solve the implicit equation for x
    eq =  @(x) asin(1./(1 + x*vinfm.^2./mu)) + asin(1./(1 + x*vinfp^2/mu)) - delta;
    rp = fzero(eq,[0.0000001 6e15]);
    
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
    vpm = sqrt(mu*(1/abs(am) + 2/rp));
    vpp = sqrt(mu*(1/abs(ap) + 2/rp));
    vp0=vpm;
    
    % Compute the delta v required at perigee burn
    dv_p = abs(vpp - vpm);
end
end