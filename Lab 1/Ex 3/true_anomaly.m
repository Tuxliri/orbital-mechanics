function f = true_anomaly(e, a, mu, N, k, t0,f0)
% computes the true anomaly for an orbit with
% fixed e,a and mu, for an N-points array of times covering k periods of
% the orbit
% 
% PROTOTYPE:
%   f = true_anomaly(e, a, mu, N, k)
%
% INPUT:
%   e [1] eccentricity                            [ - ]
%   a [1] semimajor axis                          [ L ]
%   mu [1] gravitational parameter of the primary [ L^3/T^2 ]
%   N [1] Points of calculation                   [ - ]
%   k [1] number of periods                       [ - ]
%   t0 [1] reference initial time                 [ T ]
%   f0 [1] reference initial true anomaly         [ - ]
%
% OUTPUT:
%   f [1xN] true anomaly      [ deg ]

% Calculate period

T = 2*pi*sqrt(a^3/mu);      % Period                            [ s ]

t = linspace(t0, k*T ,N);

f = zeros(1,length(t));

for i=1:length(t)
    f(i) = fsolve_kepler_eq_solver(t(i),e,a,mu,t0,f0);
end


end

