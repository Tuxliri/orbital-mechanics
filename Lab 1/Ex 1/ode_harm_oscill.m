function dy = ode_harm_oscill(~, y, omega0, gamma)
% ode_harm_oscill ODE 1st order system for the damped harmonic oscillator
%
% PROTOTYPE:
%   dy = ode_harm_oscill(t, y, omega0, gamma)
%
% INPUT:
%
%   t[1]            Time (can be omitted, the system is autonomous) [T]
%   y[2x1]          State of the oscillator (position and velocity) [L, L/T]
%   omega0[1]       Natural frequency of the damped oscillator [1/T]
%   gamma[1]        Damping coefficient [1/T]
%
% OUTPUT:
%   dy[2x1]         Derivative of the state [L, L/T]
%
% CONTRIBUTORS:
%   Davide Iafrate
%
% VERSIONS
%   2020-09-24: First version
%

% Set the derivatives of the state

dy = [ y(2)
       -2*gamma*y(2)-omega0^2*y(1) ];

end

