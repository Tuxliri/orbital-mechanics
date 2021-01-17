function a = repeating_ground_track(m, k, mu, omega_e, J2,R_e, e, i)
% REPEATING_GROUND_TRACK function that computes the required a (semi-major axis)
%   for a repeating ground track with k satellite revolutions and m Earth
%   revolutions, with the option of accounting for the SECULAR PERTURBATIONS due to J2
%
% PROTOTYPE:
%   a = repeating_ground_track(m, k, mu, omega_e, J2,R_e, e, i)
%
% INPUT:
%   m[1]        rotations of the planet (positive integer)      [-]
%   k[1]        revolutions of the satellite (positive integer) [-]
%   mu[1]       gravitational parameter                         [km^3/s^2]
%   omega_e[1]    rotation rate of the planet                   [rad/s]
%   J2[1]       Second zonal harmonic                           [-]
%   R_e[1]      Equatorial radius of Earth                      [km]
%   e[1]        Eccentricity                                    [-]
%   i[1]        Inclination                                     [rad]
%
% OUTPUT:
%   a[1]        required semi-major axis                        [km]
%
% CONTRIBUTORS:
%   Davide Iafrate      
%   Alkady Marwan
%   Pedro Bossi Núñez
%   Davide Demartini
%
% VERSIONS
%   2020-10-30: First version
%   2021-01-14: Second version

% calculate an initial guess for the semi-major axis
% from the undisturbed problem

if nargin == 4
    a = mu^(1/3) * (m/(k*omega_e))^(2/3);
    
end

if nargin == 8
    a0 = mu^(1/3) * (m/(k*omega_e))^(2/3);
    
    %% expression of mean motion n
    
    n = @(a) sqrt(mu/a^3);
    
    %% expressions of the secular evolution of RAAN, omega, M0 = M(t0)
    coeff = @(a) -((3/2)*sqrt(mu)*J2*R_e^2 / ((1-e^2)^2 * a^(7/2)));
    
    % Regression of the node
    RAAN_dot = @(a) coeff(a) * cos(i);
    
    % Advance of the perigee
    omega_dot = @(a) coeff(a)*((5/2)*(sin(i)^2) - 2);
    
    % Rate of change of mean anomaly at t0
    
    M0_dot = @(a) (1.5*sqrt(mu)*J2*R_e^2 / ((1-e^2)^1.5 * a^(7/2))) * (1-1.5*(sin(i)^2));
    
    %% Secular effects of J2 modify the nodal periods of satellite and Earth
    
    options = optimset('TolFun', 1e-6, 'Display', 'off');
    % Satellite nodal period
    
    fun = @(a) ((omega_e - RAAN_dot(a))/(n(a) + omega_dot(a) + M0_dot(a)) - m/k);
    
    a = fzero(fun,a0, options);
end

end

