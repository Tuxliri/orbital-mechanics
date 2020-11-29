function [f] = kepler_eq_solver(t, e, a, mu, t0, f0)
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
%   t0 [1] reference initial time [ T ]
%   f0 [1] reference initial true anomaly [-]
%
% OUTPUT:
%   f [1] true anomaly      [ deg ]

n = sqrt(mu/a^3);               % mean motion       [ - ]
M = n*(t-t0);                   % Mean anomaly      [ T ]

%% Reduce M to M_prime ∈ [0,2*pi]rad plus a whole number of revolutions
% 𝑘 ∈ ℤ, so that M = M_prime + k*2*pi

M_prime = wrapTo2Pi(M);
k = floor(M/(2*pi));

%% Numerically solve Kepler's equation for M_prime

E_guess = M_prime + (e*sin(M_prime))/(1- sin(M_prime + e) + sin(M_prime));
syms E
eqn = E - e*sin(E) == M_prime;
E_sol = vpasolve(eqn,E,E_guess);

%% Compute the corresponding true anomaly f ∈ [0,2*pi] rad
f = 2*atan2d(tan(E_sol/2),sqrt((1-e)/(1+e)));

if f < 0
    f = f + 360;
end

%% Add the whole number of revolutions

f = f + k*360;


end

