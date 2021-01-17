function [r,v] = propagator(r0,v0,mu,tspan,J2,R_e,date0,Cr,Psr,Am)
% PROPAGATOR function that predicts body’s orbital characteristics 
%                   at a future date given the current orbital 
%                   characteristics with optional parameters for J2
%                   effect and SRP perturbation
%                   J2: necessary parameters [J2,R_e]
%                   J2+SRP: necessary parameters [J2,R_e,date0,Cr,Psr,Am]
% 
%   
% PROTOTYPE:
%    [r,v] = propagator(r0,v0,mu,tspan,J2,R_e,date0,Cr,Psr,Am)
% 
% INPUT:
%   r0[3]            Initial radius vector             [ km ]
%   v0[3]            Initial velocity vector           [ km/s ]
%   mu[1]            Earth's gravitational parameter   [ km^3/s^2 ]
%   tspan[1 x n]     Vector of times at which the ground track will
%                    be computed                       [ s ]
%   J2[1]      Second zonal harmonic                   [-]
%   R_E[1]     Equatorial radius of Earth              [km]
%   date0[6]   initial date vector                              [date]
%   Cr[1]      reflectivity coefficient                  [-]
%   Psr[1]     Solar radiation pressure at 1AU           [N/m^2]
%   Am[1]      Area to mass ratio of the spacecraft      [m^2/km]
%
% OUTPUT:
%   r [ n x 3]      Position vector computed at each time of the tspan 
%                   vector                             [km]
%   v [ n x 3]      Velocity vector computed at each time of the tspan 
%                   vector                             [km/s]
% CONTRIBUTORS:
%   Davide Iafrate      
%   Alkady Marwan
%   Pedro Bossi Núñez
%   Davide Demartini
%
% VERSIONS
%   2020-10-16: First version
%   2021-01-14: Implemented optional parameters
%

%% Set options for ODE solver
opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','off');

%% Initial state vector
y0 = [ r0(1) v0(1) r0(2) v0(2) r0(3) v0(3) ];   

%% Compute the integration of the appropriate ode

if nargin == 4
    [~, y] = ode113(@(t,y) twobodyode(t,y,mu), tspan, y0, opts);

% J2 case
elseif nargin == 6
    [~, y] = ode113(@(t,y) twobodyode(t,y,mu,J2,R_e), tspan, y0, opts);

% J2 + SRP case
elseif nargin == 10
    [~, y] = ode113(@(t,y) twobodyode(t,y,mu,J2,R_e,date0,Cr,Psr,Am), tspan, y0, opts);
end
%% Calculate numerically the velocities, radiuses and specific energy

r = [y(:,1) y(:,3) y(:,5)];                   % Radiuses         [ km ]
v = [y(:,2) y(:,4) y(:,6)];                   % Velocities       [ km/s ]

end

