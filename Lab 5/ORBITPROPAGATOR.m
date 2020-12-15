function [t,kep] = ORBITPROPAGATOR(t0,s0,tspan)
%ORBITPROPAGATOR Propagator for an earth orbit accounting for gaussian
%                perturbations effects
% 
% PROTOTYPE:
%   [t,kep] = ORBITPROPAGATOR(t0,s0,tspan)
% 
% INPUT:
%   t0[1]      Initial time                                     [s]
%   s0[6x1]    Initial state of the keplerian elements
%              [a,e,i,RAAN,omega,f]            [km,-,rad,rad,rad,rad]
%   tspan[N]   Vector of times for function integration         [s]
%
% OUTPUT:
%   t[N]       Vector of propagated times                       [s]
%   kep[Nx6]   Matrix of keplerian elements at each evaluated time instant
%              [a,e,i,RAAN,omega,f]            [km,-,rad,rad,rad,rad]
%
% CONTRIBUTORS:
%   Davide Iafrate
%
% VERSIONS
%   2020-12-14: First version
%

%% useful constants
const = astroConstants([23 13 9]);
R_E = const(1);
muE = const(2);
J2 = const(3);

% Set options for the solver
options = odeset('Reltol',1e-16,'Stats','off');

% Create a function handle for the perturbing accelerations
a_per_rsw = @(t,s) a_tot_rsw(t,s,muE,J2,R_E);

% Create a function handle for the ode113 solver
odefun = @(t,kep) EOM_RSW(t,kep,muE,a_per_rsw);

% Numerically solve the EOM
[t,kep] = ode113(odefun,tspan,s0,options);


