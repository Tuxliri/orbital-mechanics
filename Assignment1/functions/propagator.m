function [r,v] = propagator(r0,v0,mu,tspan)
% PROPAGATOR function that predicts body’s orbital characteristics 
%                   at a future date given the current orbital 
%                   characteristics.
%   
% PROTOTYPE:
%    [v,r] = orbit_propagation(r0,v0,mu,tspan)
% 
% INPUT:
%   r0[3]            Initial radius vector             [ km ]
%   v0[3]            Initial velocity vector           [ km/s ]
%   mu[1]            Planetary gravitational parameter [ km^3/s^2 ]
%   tspan[1 x n]     Vector of times at which the ground track will
%                    be computed                       [ s ]
% OUTPUT:
%   r [ n x 3]      Position vector computed at each time of the tspan 
%                   vector                             [km]
%   v [ n x 3]      Velocity vector computed at each time of the tspan 
%                   vector                             [km/s]
%
% CONTRIBUTORS:
%   Alkady Marwan
%   Bossi Núñez Pedro
%   Davide Demartini
%   Davide Iafrate
%
% VERSIONS
%   2020-10-16: First version
%

% Set options for ODE solver
opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','on');

% Initial state vector
y0 = [ r0(1) v0(1) r0(2) v0(2) r0(3) v0(3) ];   

% Compute the integration of the ode

[t, y] = ode113(@(t,y) twobodyode(t,y,mu), tspan, y0, opts);

% Calculate numerically the velocities, radiuses and specific energy

r = [y(:,1) y(:,3) y(:,5)];                   % Radiuses         [ km ]
v = [y(:,2) y(:,4) y(:,6)];                   % Velocities       [ km/s ]

end
