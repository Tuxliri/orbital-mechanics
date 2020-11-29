function Y = an_harmonic_oscill(t, y0, omega0, gamma)
%an_harmonic_oscill analytical solution for the UNDERDAMPED harmonic
%oscillator
% PROTOTYPE
%   [x_t, v_t] = an_harmonic_oscill(t, y0(1), y0(2), omega0, gamma)
%
% INPUT:
%   t[1]            Time [T]
%   y0(1)[2x1]      Initial position and velocity [L, L/T]
%   omega0[1]       Natural frequency of the damped oscillator [1/T]
%   gamma[1]        Damping coefficient [1/T]
%
% OUTPUT:
%   Y[2x1]          Position and velocity [L, L/T]
%
% CONTRIBUTORS:
%   Davide Iafrate
%
% VERSIONS
%   2020-09-24: First version
%

% Calculate omega [1/T]

omega = sqrt(omega0.^2 - gamma.^2);

% Calculate the position and velocity

Y = [ (exp(-gamma.*t)).*(y0(1).*cos(omega.*t) + (y0(2) + gamma.*y0(1)).*sin(omega.*t)/omega)
      (exp(-gamma.*t)).*(y0(2).*cos(omega.*t) - ((omega0.^2)*y0(1) + gamma.*y0(2)).*sin(omega.*t)/omega)];
end

