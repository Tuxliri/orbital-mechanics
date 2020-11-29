function a = repeating_ground_track_j2(m, k, mu, omega_e)
% REPEATING_GROUND_TRACK function that computes the required a (semi-major axis)
%   for a repeating ground track with k satellite revolutions and m Earth revolutions.
%   
% PROTOTYPE:
%   a = repeating_ground_track(m, k, mu, omega)
% 
% INPUT:
%   m[1]        rotations of the planet (positive integer)      [-]
%   k[1]        revolutions of the satellite (positive integer) [-]
%   mu[1]       gravitational parameter                         [km^3/s^2]
%   omega_e[1]    rotation rate of the planet                     [rad/s]
%   J2[1]       Second zonal harmonic         [-]
%   R_e[1]      Equatorial radius of Earth    [km]
%   e[1]        Eccentricity                  [-]
%   i[1]        Inclination                   [rad]
%
% OUTPUT:
%   a[1]        required semi-major axis                        [km]  

% calculate an initial gues for the semi-major axis
% from the undisturbed problem

a = mu^(1/3) * (m/(k*omega_e))^(2/3);

end

