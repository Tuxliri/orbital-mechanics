function dy = twobodyode_j2_SRP(t,y,mu,J2,R_e,date0)
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
%   R_e[1]     Equatorial radius of Earth                       [km]
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

DAY2SECS = 24*3600;

% Spacecraft characteristics
Cr = 1.2;                       % [-]
Psr = 4.5e-6;                   % [N/m^2]   at 1AU
Am = 4;                         % [m^2/km]

% Calculate radius
r = norm(y(1:2:6));

r_vec = y(1:2:6);
v_vec = y(2:2:6);

[s(1),s(2),s(3),s(4),s(5),s(6)] = car2kep(r_vec,v_vec,mu);

% This should be mjd2000_0, at the initial time of integration, then we
% should puy into the function to calculate the SRP acceleration the
% mjd2000 + t/(24*3600) to get the actual exact mjd2000 date at each time
% instant of integration -> Done inside the a_tot_rsw function

mjd2000_0 = date2mjd2000(date0);

[~,a_ECI] = a_SRP(s,Cr,Psr,R_e,Am,mjd2000_0 + t / DAY2SECS);

% Set the derivatives of the state
dy = [ y(2); - mu*y(1)/r^3 + ((1.5*J2*mu*R_e^2)/r^4) * y(1) / r * (5*y(5)^2/r^2 - 1) + a_ECI(1);
       y(4); - mu*y(3)/r^3 + ((1.5*J2*mu*R_e^2)/r^4) * y(3) / r * (5*y(5)^2/r^2 - 1) + a_ECI(2);
       y(6); - mu*y(5)/r^3 + ((1.5*J2*mu*R_e^2)/r^4) * y(5) / r * (5*y(5)^2/r^2 - 3) + a_ECI(3)];
end
