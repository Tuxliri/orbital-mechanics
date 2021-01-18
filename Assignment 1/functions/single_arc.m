function [DV, DV1, DV2, v_t_start, v_t_end]= single_arc(t1,t2,mu,ri,rf,vi,vf)
% Get the DV needed for the transfer between two orbits
% and the vectors of arrival position and velocity
%
% PROTOTYPE:
%    [DV, DV1, DV2, r_t_end, v_t_end] = single_arc(t1,t2,mu,ri,rf,vi,vf)
%   
% INPUT :
%   t1[1]      departure in mjd2000 [days]
%   t2[1]      arrival time in mjd2000   [days]
%   mu[1]      gravitational constant of  the central body [km^3/s^2]
%   ri[1]      initial position [km]
%   vi[1]      initial velocity [km/s]
%   rf[1]      final position [km]
%   vf[1]      final velocity [km/s]
%
% OUTPUT:
%   DV[1]       total cost of maneuvre [km/s]
%   DV1[1]      cost of departure maneuvre [km/s]
%   DV2[1]      cost of arrival maneuvre [km/s]
%   v_t_start[3]     departure velocity vector of the transfer orbit [km]
%   v_t_end[3]     arrival velocity vector of the transfer orbit [km]
%
% CALLED FUNCTIONS:
%   lambertMR
% 
% CONTRIBUTORS:
%   Alkady Marwan
%   Bossi Núñez Pedro
%   Davide Demartini
%   Davide Iafrate
%
% VERSIONS
%   2020-11-20: First version

DAY2SECS = 24*3600;
TOF = (t2 - t1)*DAY2SECS;         % [s]

% Call Lambert solver
[~,~,~,~,v1L,v2L,~,~] =...
    lambertMR(ri,rf,TOF,mu,0,0,0,1);

% Compute the total cost of the maneuvre
DV1 = norm(vi - v1L);
DV2 = norm(v2L - vf);
DV = DV1 + DV2;

% Output the departure and arrival velocities along the transfer orbit
v_t_start = v1L;
v_t_end = v2L;
end