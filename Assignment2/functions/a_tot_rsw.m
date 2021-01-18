function acc = a_tot_rsw(t,s,muE,J2,R_E,Cr,Psr,Am,mjd2000)
%a_tot_rsw Computation of the vector of disturbing accelerations in the RSW
%           (Radial-transveral-out of plane) reference frame
% PROTOTYPE:
%   acc = a_tot_rsw(t,s,muE,J2,R_E,Cr,Psr,Am,mjd2000)
% 
% INPUT:
%   t[1]    time                                [s]
%   s[6]    state vector containing the keplerian elements
%               [a,e,i,RAAN,omega,f]            [km,-,rad,rad,rad,rad]
%   R_E[1]  mean radius of the Earth            [km]
%   muE[1]  Earth gravitational parameter       [km^3/s^2]
%   J2[1]   Second zonal harmonic               [-]
%   Cr[1]       Refelctivity coefficient        [-]
%   Psr[1]      Solar radiation pressure               [N/m^2]
%   Am[1]       area to mass ratio of the S/C       [m^2/kg]
%   mjd2000[1]  mjd2000 date                        [-]
%
% 
% OUTPUT:
%   acc[3]  vector of perturbing accelerations [ar,as,aw]
% 
% CONTRIBUTORS:
%   Davide Iafrate      
%   Alkady Marwan
%   Pedro Bossi Núñez
%   Davide Demartini
%
% VERSIONS
%   14-12-2020: First version

DAY2SECS = 24*3600;
date_mjd = mjd2000 + t /DAY2SECS;
[acc_SRP,~] = a_SRP(s,Cr,Psr,R_E,Am,date_mjd);
acc = a_j2_RSW(s,muE,J2,R_E) +  acc_SRP;    % + .....
end

