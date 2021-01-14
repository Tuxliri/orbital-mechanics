function dx = twobody( ~, x, mu )
% ODE system for the keplerian orbit
%
% PROTOTYPE:
%   dx = keplerian_orbit_ode( t, x, mu )
%
% INPUT:
%   t[1] time                      [ s ]
%   x[6x1] state vector            [ km, km/s ]
%   mu planetary constant          [ km^3/s^2 ]
%
% OUTPUT:
%   dx[6x1] Derivative of the state vector [ km/s, km/s^2 ]
%
% CONTRIBUTORS:
%   Davide Demartini
%
% VERSIONS
%   12-12-2020: First version

% Find r
r = sqrt(x(1)^2+x(3)^2+x(5)^2);

% Set the derivatives of the state
dx = [ x(2); - mu*x(1)/r^3;
       x(4); - mu*x(3)/r^3;
       x(6); - mu*x(5)/r^3; ];
end