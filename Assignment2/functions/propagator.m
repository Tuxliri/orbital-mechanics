function [r,v,eps] = propagator(r0,v0,mu,tspan)
% PROPAGATOR function that predicts body’s orbital characteristics 
%   at a future date given the current orbital characteristics.
%   
% PROTOTYPE:
%    [v,r] = orbit_propagation(r0, v0, mu, tspan)
% 
% INPUT:
%   r0[3]            Initial radius vector               [ km ]
%   v0[3]            Initial velocity vector             [ km/s ]
%   mu[1]            Planetary gravitational parameter   [ km^3/s^2 ]
%   tspan[1 x n]     Vector of times at which the ground track will
%                    be computed                         [ s ]
% OUTPUT:
%   r   [ n x 3]    Position vector computed at each time of the tspan 
%                   vector                               [ km ]
%   v   [ n x 3]    Velocity vector computed at each time of the tspan 
%                   vector                               [km/s]
%   eps [ n x 1]    Specific energy vector computed at each time of the
%                   tspan vector                         [ km^2/s^2 ]
%
% CONTRIBUTORS:
%   Davide Iafrate
%   Davide Demartini
%
% VERSIONS
%   2020-10-16: First version
%   2020-12-15: Second version

% Set ODE options
opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','on');

% Iitial state
x0 = [r0(1) v0(1) r0(2) v0(2) r0(3) v0(3)];

% ODE integration
[~, x] = ode113(@(t,x) twobody(t,x,mu), tspan, x0, opts);

% Calculate numerically the velocities, radiuses and specific energy
% vectors
r = [x(:,1) x(:,3) x(:,5)];                           % Radiuses        [ km ]
v = [x(:,2) x(:,4) x(:,6)];                           % Velocities      [ km/s ]
eps = zeros(1,length(v));
for i=1:length(v)
    eps(i) = 0.5*(norm(v(i,:)))^2 - mu/norm(r(i,:));  % Specific Energy [km^2/s^2]
end

end