function [r_vec,v_vec] = kep2car(a,e,i,RAAN,omega,f,mu)

% Solution to: [a,e,f,RAAN,omega,i,mu]~~>[r,v]
%
% ATTENTION: 
% 1)   CHECK if it works with parabolic and hyperbolic orbits;
% PROTOTYPE:
%    [r_vec,v_vec] = kep2car(a, e, i, RAAN, omega, f, mu) 
%   
% INPUT:
%   a[1]        semi-major axis     [ km ]
%   e[1]        eccentricity        [ - ]
%   i[1]        inclination         [ rad ]
%   RAAN[1]     right ascension of the ascending node [ rad ]
%   omega[1]    argument of perigee [ rad ]
%   f[1]        true anomaly        [ rad ]
%   mu[1]      gravitational paramer [ km^3/s^2 ]
%
% OUTPUT:
%   r_vec[3]   position vector     [ km ]
%   v_vec[3]   velocity vector     [ km/s ]
%   
% CONTRIBUTORS:
%   Davide Demartini
%   Davide Iafrate
%
% VERSIONS
%   2020-10-12: First version
%   2020-11-24: Bug fixing, code cleanup
%

% Rotation Matrixes

% rotation around K of an angle RAAN

R_RAAN=[cos(RAAN)   sin(RAAN) 0
        -sin(RAAN)  cos(RAAN) 0
         0          0         1];

% rotation around the line of nodes of an angle i
R_i=[1 0 0;
    0 cos(i) sin(i)
    0 -sin(i) cos(i)]; 

% rotation around h of amount of an angle omega
R_omega=[cos(omega) sin(omega) 0
    -sin(omega) cos(omega) 0
    0 0 1]; % Perigee anomaly

% Trasformation matrix RF Perifocal ~~> Geocentric
% R = R_omega*R_i*R_RAAN  GE~~>PF
% T = R'                  PF~~>GE
T = (R_omega*R_i*R_RAAN)';

% Parameter p
p=a*(1-e^2);

% Compute position components in PF RF
r_pf = p/(1+e*cos(f));

r_pf_vect= r_pf.*[cos(f); sin(f); 0];

% Compute velocity components in PF RF

v_pf_vect= sqrt(mu/p) .* [-sin(f); (e+cos(f)); 0];

% Compute position and velocity components in Geocentric RF
% to convert vectors from the PF ref frame to the geocentric one
% we multiply the vectors in the PF frame by the matrix T
r_vec = T*r_pf_vect;
v_vec = T*v_pf_vect;

r_vec = r_vec';
v_vec = v_vec';

end
