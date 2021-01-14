function dy = twobodyode(t,y,mu,J2,R_e,date0,Cr,Psr,Am)
% Restructed two body problem ODE function with optional parameters for J2
% effect and SRP perturbation
% J2: necessary parameters [J2,R_e]
% J2+SRP: necessary parameters [J2,R_e,date0,Cr,Psr,Am]
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
%   Cr[1]      reflectivity coefficient                         [-]
%   date0[6]   initial date vector                              [date] 
%   Psr[1]     Solar radiation pressure at 1AU                  [N/m^2]
%   Am[1]      Area to mass ratio of the spacecraft             [m^2/km]
%
% OUTPUT:
%   dy[6x1]         Derivative of the state [L, L/T]
%
% CONTRIBUTORS:
%   Davide Iafrate
%
% VERSIONS
%   2020-09-24: First version
%   2021-01-14: Implemented optional parameters
%

DAY2SECS = 24*3600;

% Calculate radius
r = norm(y(1:2:6));

r_vec = y(1:2:6);
v_vec = y(2:2:6);

[s(1),s(2),s(3),s(4),s(5),s(6)] = car2kep(r_vec,v_vec,mu);

% Check if J2 arguments are passed
if nargin == 3
    J2 = 0;
    R_e = 0;
    
end
    
% Perform check on SRP required
if nargin == 9
    mjd2000_0 = date2mjd2000(date0);
    
    % Calculate the acceleration in the ECI frame
    [~,a_ECI] = a_SRP(s,Cr,Psr,R_e,Am,mjd2000_0 + t / DAY2SECS);

else
    a_ECI = [0 0 0];
end

% Set the derivatives of the state
dy = [ y(2); - mu*y(1)/r^3 + ((1.5*J2*mu*R_e^2)/r^4) * y(1) / r * (5*y(5)^2/r^2 - 1) + a_ECI(1);
       y(4); - mu*y(3)/r^3 + ((1.5*J2*mu*R_e^2)/r^4) * y(3) / r * (5*y(5)^2/r^2 - 1) + a_ECI(2);
       y(6); - mu*y(5)/r^3 + ((1.5*J2*mu*R_e^2)/r^4) * y(5) / r * (5*y(5)^2/r^2 - 3) + a_ECI(3)];
end
