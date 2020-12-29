  
function [alpha, delta, lon, lat] = GroundTrack(state_vec,gw_lon0,t,omega_e,mu,t0);
%  function groundTrack that computes the ground track
%  of an orbit
%
% PROTOTYPE:
%    [alpha, delta, lon, lat] = groundtrack(state_vec, gw_lon0, t, omega, mu, t0);
%   
% INPUT:
%   state_vec[2 | 6]        Cartesian or keplerian elements at initial time [ - ]
% 
%       [a, e_mag, i, RAAN, omega, f] if Keplerian:
%
%           a       [ 1 ]        Semi-major axis                       [ km ]
%           e_mag   [ 1 ]        Eccentricity                          [ - ]
%           i       [ 1 ]        Inclination                           [ rad ]
%           RAAN    [ 1 ]        Right ascension of the ascending node [ rad ]
%           omega   [ 1 ]        Argument of perigee                   [ rad ]
%           f       [ 1 ]        True anomaly                          [ rad ]
%       
%       [r,v] if Cartesian:
% 
%           r_vec   [ 3 x 1 ]    Position vector                       [ km/s ]
%           v_vec   [ 3 x 1 ]    Velocity vector                       [ km/s ]
%  
%   gw_longitude0   [ 1 ]        Longitude of Greenwich meridian at initial 
%                                  time                                [ rad ]
%   t               [ n ]        Vector of times at which the ground track 
%                                   will be computed                   [ s ]
%   omega_e         [ 1 ]        Earth's angular velocity              [ rad/s ]
%   mu              [ 1 ]        Gravitational paramer                 [ km^3/s^2 ]
%   t0              [ 1 ]        Initial time                          [ s ]
% 
% OUTPUT:
%   alpha           [ n x 1 ]    Right ascension in Earth Centered 
%                                   Equatorial Inertial frame          [ rad ]
%   delta           [ n x 1 ]    Declination vector                    [ rad ]
%   lon             [ n x 1 ]    Longitude with respect to rotating 
%                                   Earth                              [ rad ]
%   lat             [ n x 1 ]    Latitude with respect to rotating 
%                                   Earth  [ rad ]
% 
% CONTRIBUTORS:
%   Davide Iafrate
%   Davide Demartini
%
% VERSIONS
%   16-1012020: First version
%   15-12-2020: Second version

% Input check
if length(state_vec) == 6
    a = state_vec(1);
    e = state_vec(2);
    i = state_vec(3);
    RAAN = state_vec(4);
    omega = state_vec(5);
    f = state_vec(6);
    [r0, v0] = kep2car(a,e,i,RAAN,omega,f,mu);
else
    r0 = state_vec(1);
    v0 = state_vec(2);
end

% Orbit propagation
[r] = propagator(r0,v0,mu,t);

% Conversion to RA and declination
[delta, alpha] = r2RAdelta(r);

% Conversion to latitude and longitude
[lon, lat] = RAdelta2latlon(alpha, delta, gw_lon0, omega_e, t, t0);
end