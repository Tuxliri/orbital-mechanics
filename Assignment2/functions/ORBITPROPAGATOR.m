function [t,kep] = ORBITPROPAGATOR(t0,s0,tspan,date0)
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
%   2021-01-02: Implemented SRP acceleration

%% useful constants
const = astroConstants([23 13 9]);
R_E = const(1);
muE = const(2);
J2 = const(3);

% Spacecraft characteristics
Cr = 1.2;                       % [-]
Psr = 4.5e-6;                   % [N/m^2]   at 1AU
Am = 4;                         % [m^2/km]

% This should be mjd2000_0, at the initial time of integration, then we
% should puy into the function to calculate the SRP acceleration the
% mjd2000 + t/(24*3600) to get the actual exact mjd2000 date at each time
% instant of integration -> Done inside the a_tot_rsw function

mjd2000_0 = date2mjd2000(date0);

% Set options for the solver
options = odeset('Reltol',1e-16,'Stats','off');

% Create a function handle for the perturbing accelerations
a_per_rsw = @(t,s) a_tot_rsw(t,s,muE,J2,R_E,Cr,Psr,Am,mjd2000_0);

% Create a function handle for the ode113 solver
odefun = @(t,kep) EOM_RSW(t,kep,muE,a_per_rsw);

% Numerically solve the EOM
[t,kep] = ode113(odefun,tspan,s0,options);


