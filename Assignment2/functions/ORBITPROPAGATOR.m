function [t,kep] = ORBITPROPAGATOR(t0,s0,tspan,date0,J2,Cr,Psr,Am)
%ORBITPROPAGATOR Propagator for an earth orbit accounting for gaussian
%                perturbations effects
% 
% PROTOTYPE:
%   [t,kep] = ORBITPROPAGATOR(t0,s0,tspan,date0,J2,Cr,Psr,Am)
% 
% INPUT:
%   t0[1]      Initial time                                     [s]
%   s0[6x1]    Initial state of the keplerian elements
%              [a,e,i,RAAN,omega,f]            [km,-,rad,rad,rad,rad]
%   tspan[N]   Vector of times for function integration         [s]
%   date0[6]   Initial date vector                              
%   J2[1]      Second zonal harmonic                            [-]
%   Cr[1]      reflectivity coefficient                  [-]
%   Psr[1]     Solar radiation pressure at 1AU           [N/m^2]
%   Am[1]      Area to mass ratio of the spacecraft      [m^2/km]
%
% OUTPUT:
%   t[N]       Vector of propagated times                       [s]
%   kep[Nx6]   Matrix of keplerian elements at each evaluated time instant
%              [a,e,i,RAAN,omega,f]            [km,-,rad,rad,rad,rad]
%
% CONTRIBUTORS:
%   Davide Iafrate      
%   Alkady Marwan
%   Pedro Bossi Núñez
%   Davide Demartini
%
% VERSIONS
%   2020-12-14: First version
%   2021-01-02: Implemented SRP acceleration

%% useful constants
const = astroConstants([23 13 9]);
R_E = const(1);
muE = const(2);

mjd2000_0 = date2mjd2000(date0);

% Set options for the solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14,'Stats','off');

% Create a function handle for the perturbing accelerations
a_per_rsw = @(t,s) a_tot_rsw(t,s,muE,J2,R_E,Cr,Psr,Am,mjd2000_0);

% Create a function handle for the ode113 solver
odefun = @(t,kep) EOM_RSW(t,kep,muE,a_per_rsw);

% Numerically solve the EOM
[t,kep] = ode113(odefun,tspan,s0,options);


