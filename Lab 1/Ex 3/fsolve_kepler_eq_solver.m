function [f] = fsolve_kepler_eq_solver(t, e, a, mu, t0, f0)
% Solver for Kepler's equation
% 
% PROTOTYPE:
%   [f] = kepler_eq(t, e, a, mu, t0, f0)
%
% INPUT:
%   t [1] time             [ T ]
%   e [1] eccentricity      [ - ]
%   a [1] semimajor axis    [ L ]
%   mu [1] gravitational parameter of the primary [ L^3/T^2 ]
%   t0 [1] reference initial time of f0 [ T ]
%   f0 [1] reference initial true anomaly [rad]
%
% OUTPUT:
%   f [1] true anomaly      [ deg ]
%
% LIMITATIONS:
%   elliptical orbits only! have to fix it for hyperbolic and parabolic
%   orbits ( See UNIVERSAL VARIABLES (Chapter 3.7 of Curtis' book)

if nargin < 5
    % Initial values of t0 and f0 assumed to be at perigee passage: 
    % t0 = tp and f0 = 0
    f0 = 0;
end

n = sqrt(mu/a^3);               % mean motion       [ - ]
M = n*(t-t0);                   % Mean anomaly      [ T ]

%% compute the initial Eccentric Anomaly at time t0
E0 = 2*atan2(tan(f0/2),sqrt((1+e)/(1-e)));

%% Reduce M to M_prime ∈ [0,2*pi]rad plus a whole number of revolutions
% k ∈ Z, so that M = M_prime + k*2*pi

M_prime = wrapTo2Pi(M);
k = floor(M/(2*pi));
%% Numerically solve Kepler's equation for M_prime

E_guess = M_prime + (e*sin(M_prime))/(1- sin(M_prime + e) + sin(M_prime));
eqn = @(E) E - e*sin(E) - (E0 - e*sin(E0)) - M_prime;

options = optimoptions('fsolve','Display','none');
E_sol = fsolve(eqn,E_guess, options);

%% Compute the corresponding true anomaly f ∈ [0,2*pi] rad
f = 2*atan2d(tan(E_sol/2),sqrt((1-e)/(1+e)));

if f < 0
    f = f + 360;
end

%% Add the whole number of revolutions

f = f + k*360 ;


end

