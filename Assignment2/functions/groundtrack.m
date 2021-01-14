function [alpha, delta, lon, lat] = groundtrack_J2(state_vec, gw_longitude0, t, omega_e, mu, t0,j2,R_e)
%  function groundTrack that computes the ground track
%  of an orbit
%
% PROTOTYPE:
%    [alpha, delta, lon, lat] = groundtrack(state_vec, gw_longitude, t, omega_e, mu, t0);
%   
% INPUT:
%   state_vec[2 | 6]        Cartesian or keplerian elements at initial time [ - ]
% 
%       [a, e_mag, i, RAAN, omega, f] if Keplerian:
%
%           a[1]        semi-major axis     [ km ]
%           e_mag[1]        eccentricity    [ - ]
%           i[1]        inclination         [ rad ]
%           RAAN[1]     right ascension of the ascending node [ rad ]
%           omega[1]    argument of perigee [ rad ]
%           f[1]        true anomaly        [ rad ]
%       
%       [r,v] if Cartesian:
% 
%           r_vec[3]   position vector     [ km/s ]
%           v_vec[3]   velocity vector     [ km/s ]
%  
%   gw_longitude0[1]        Longitude of Greenwich meridian at initial time [ rad ]
%   t[n]        Vector of times at which the ground track will be computed [ s ]
%   omega_e[1]  Earth's angular velocity      [ rad/s ]
%   mu[1]       Gravitational paramer         [ km^3/s^2 ]
%   t0[1]       Initial time                  [ s ]
%   J2[1]       Second zonal harmonic         [-]
%   R_e[1]      Equatorial radius of Earth    [km]
%
% OUTPUT:
%   alpha[1]   right ascension in Earth Centered Equatorial Inertial frame     [ rad ]
%   delta[]    velocity vector     [ rad ]
%   lon[1]     longitude with respect to rotating Earth [ rad ]
%   lat[1]     latitude with respect to rotating Earth  [ rad ]
% 
% CONTRIBUTORS:
%   
%   Davide Iafrate
%
% VERSIONS
%   2020-10-16: First version
%

%% check if the initial state vector is given in Cartesian or Keplerian elements

if length(state_vec) == 6   % Keplerian elements
    a = state_vec(1);
    e = state_vec(2);
    i = state_vec(3);
    RAAN = state_vec(4);
    omega = state_vec(5);
    f0 = state_vec(6);
    
    [r0,v0] = kep2car(a,e,i,RAAN,omega,f0,mu); 
else
    
    r0 = state_vec(1);
    v0 = state_vec(2);
end

%% Orbit propagation

[r,~] = propagator(r0,v0,mu,t,j2,R_e);

% Conversion to RA and declination
[delta, alpha] = r2RAdelta(r);

% Conversion to latitude and longitude
[lon, lat] = RAdelta2latlon(alpha, delta, gw_longitude0, omega_e, t, t0);


end

