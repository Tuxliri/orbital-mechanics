function dy = twobodyode(t,y,mu)
% Restructed two body problem ODE function
% 
% PROTOTYPE:
%   dy = twobodyode(t,y)
% 
% INPUT:
%   t[1]       Time                                             [T]
%   y[6x1]     State of the oscillator 
%              (position and velocity vectors)                  [L, L/T]
%   mu[1]      Planetary constant
%
% OUTPUT:
%   dy[6x1]         Derivative of the state [L, L/T]
%
% CONTRIBUTORS:
%   Alkady Marwan
%   Bossi Núñez Pedro
%   Davide Demartini
%   Davide Iafrate
%
% VERSIONS
%   2020-09-24: First version

% Calculate radius
r = sqrt(y(1)^2 + y(3)^2 + y(5) ^2);

% Set the derivatives of the state
dy = [ y(2); - mu*y(1)/r^3;
       y(4); - mu*y(3)/r^3;
       y(6); - mu*y(5)/r^3; ];
end
