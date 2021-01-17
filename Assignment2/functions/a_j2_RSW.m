function a = a_j2_RSW(s,mu,J2,R_E)
%a_j2_RSW this function calculates the acceleration in the RSW
%   (Radial-transveral-out of plane) reference frame, returning a vector of
%   the three components
% 
% PROTOTYPE:
%   a = a_j2_RSW(mu,s,J2,R_E)
% 
% INPUT:
%   s[6]        state vector containing the keplerian elements
%               [a,e,i,RAAN,omega,f]                 [km,-,rad,rad,rad,rad]
%   mu[1]  Earth gravitational parameter       [km^3/s^2]
%   J2[1]   Second zonal harmonic               [-]

%   R_E[1]  mean radius of the Earth            [km]
% 
% OUTPUT:
%   a[3]    vector of perturbing accelerations [ar,as,aw]       [km/s^2]
% 
% CONTRIBUTORS:
%   Davide Iafrate      
%   Alkady Marwan
%   Pedro Bossi Núñez
%   Davide Demartini
%
% VERSIONS
%   2020-10-24: First version
%   14-12-2020: Second version

% Extract the keplerian elements from the state vector
a = s(1);
e = s(2);
i = s(3);
RAAN = s(4);
omega = s(5);
f = s(6);

% Calculate the needed orbital values
p = a*(1 - e^2);
r = p / (1+e*cos(f));

%  Calculate the vector of acceleration in the RSW frame
a = -3/2 * J2*mu*R_E^2/r^4 .* [1-3*sin(i)^2 *sin(f+omega)^2
                            sin(i)^2*sin(2*(f+omega))
                            sin(2*i)*sin(f+omega)];

end

