function [dv_p,delta_tot,a,e,delta,Delta,rp,vp] = flybypow(vinfM,vinfP,mu)
%FLYBYUNPOW function that computes an unpowered flyby manoeuvre
%
% PROTOTYPE:
%   [dv_p, delta_tot, a, e, delta, Delta, rp, vp] = flybypow(vinfM, vinfP, mu, rp_min)
%
% INPUT:
%   vinfM     [ 3 x 1 ]  incoming velocity at infinity in the planetocentric
%                    reference frame                           [ km/s ]
%   vinfP     [ 3 x 1 ]  outgoing velocity at infinity in the planetocentric
%                    reference frame                           [ km/s ]
%   mu        [ 1 x 1 ]  gravitational parameter of the planet [ km^3/s^2 ]
%
% OUTPUT:
%   dv_p      [ 1 x 1 ]  cost of the perigee burn maneuvre         [ km/s ]
%   delta_tot [ 1 x 1 ]  total turning angle                       [ rad ]
%   a         [ 2 x 1 ]  semi-major axis before and after perigee burn [ km ]
%   e         [ 2 x 1 ]  eccentricities before and after perigee burn  [ - ]
%   delta     [ 2 x 1 ]  hyperbola's turning angls                 [ rad ]
%   Delta     [ 1 x 1 ]  arrival impact parameter                  [ km ]
%   rp        [ 1 x 1 ]  perigee radius                            [ km ]
%   vp        [ 2 x 1 ]  velocity at perigee passage               [ km/s ]
%
% CONTRIBUTORS:
%   Alkady Marwan
%   Bossi Nunez Pedro
%   Demartini Davide
%   Iafrate Davide
%
% VERSIONS
%   2020-12-04: First version


vinfm = norm(vinfM);
vinfp = norm(vinfP);

% Calculate the turning angle of the maneuvre
delta_tot = acos(dot(vinfM,vinfP)/(vinfm*vinfp));

% FIX: USE FLYBY-UNPOW?
if delta_tot < 1e-6
    dv_p = 0;
    delta_tot = 0;
    a = 0;
    e = 0;
    Delta = 0;
    rp = 0;
    vp0 = 0;
else
    % Solve the implicit equation for x
    eq =  @(x) asin(1./(1 + x*vinfm.^2./mu)) + asin(1./(1 + x*vinfp^2/mu)) - delta_tot;
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
    Delta = abs(a).*sqrt(e.^2-1);
    
    % Hyperbola turning angles
    delta = 2.*asin(1./e);
    
    % Compute the velocities of the two hyperbolic arcs at perigee
    vpm = sqrt(mu*(1/abs(am) + 2/rp));
    vpp = sqrt(mu*(1/abs(ap) + 2/rp));
    vp = [vpm vpp];
    
    % Compute the delta v required at perigee burn
    dv_p = abs(vpp - vpm);
end
end