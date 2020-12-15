function acc = a_tot_rsw(t,s,muE,J2,R_E)
%a_tot_rsw Computation of the vector of disturbing accelerations in the RSW
%           (Radial-transveral-out of plane) reference frame
% PROTOTYPE:
%   acc = a_tot_rsw(t,s,R_E,muE,J2)
% 
% INPUT:
%   t[1]    time                                [s]
%   s[6]    state vector containing the keplerian elements
%               [a,e,i,RAAN,omega,f]            [km,-,rad,rad,rad,rad]
%   R_E[1]  mean radius of the Earth            [km]
%   muE[1]  Earth gravitational parameter       [km^3/s^2]
%   J2[1]   Second zonal harmonic               [-]
% 
% OUTPUT:
%   acc[3]  vector of perturbing accelerations [ar,as,aw]
% 
% CONTRIBUTORS:
%   Davide Iafrate      14-12-2020


acc = a_j2_RSW(s,muE,J2,R_E); % + .....
end

