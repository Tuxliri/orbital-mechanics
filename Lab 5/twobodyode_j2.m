function dy = twobodyode_j2(t,y,mu,J2,R_e)
% Restructed two body problem ODE function
% 
% PROTOTYPE:
%   dy = twobodyode(t,y)
% 
% INPUT:
%   t[1]       Time                                             [T]
%   y[6x1]     State of the oscillator 
%              (position and velocity vectors)                  [L, L/T]
%   mu[1]      Planetary constant                               [L^3/T^2]
%   J2[1]      Second zonal harmonic                            [-]
%   R_e[1]     
%
% OUTPUT:
%   dy[6x1]         Derivative of the state [L, L/T]
%
% CONTRIBUTORS:
%   Davide Iafrate
%
% VERSIONS
%   2020-09-24: First version
%

% Calculate radius
r = norm(y(1:2:6));

% Set the derivatives of the state
dy = [ y(2); - mu*y(1)/r^3 + ((1.5*J2*mu*R_e^2)/r^4) * y(1) / r * (5*y(5)^2/r^2 - 1);
       y(4); - mu*y(3)/r^3 + ((1.5*J2*mu*R_e^2)/r^4) * y(3) / r * (5*y(5)^2/r^2 - 1);
       y(6); - mu*y(5)/r^3 + ((1.5*J2*mu*R_e^2)/r^4) * y(5) / r * (5*y(5)^2/r^2 - 3); ];
end

