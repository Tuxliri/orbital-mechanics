function [DV, T1, T2, DV1, DV2, r_tf0, v_tf0]...
    = tfixed_design(ri,vi,rf,vf,t1,t2,mu)

% Get the DV needed for the transfer between two orbits
% and the respective times of departure and arrival
% PROTOTYPE:
%    [DV, T1, T2] = tfixed_design(ri,vi,rf,vf,t1,t2)
%   
% INPUT :
%   ri[1]      initial position [km]
%   vi[1]      initial velocity [km/s]
%   rf[1]      final position [km]
%   vf[1]      final velocity [km/s]
%   mu[1]      gravitational constant of  the central body [km^3/s^2]
%   t1         departure time [s]
%   t2         arrival time   [s]
%
% OUTPUT:
%   DV[1]       total cost of maneuvre [km/s]
%   DV1[1]      cost of departure maneuvre [km/s]
%   DV2[1]      cost of arrival maneuvre [km/s]
%   T1[1]       initial time [s]
%   T2[1]       final time   [s]
%   r_tf_orbit[npoints_orbit] position vector of the transfer orbit [km]
%   v_tf_orbit[npoints_arc]     position vector of the transfer orbit [km]
% 
% CONTRIBUTORS:
%   
%   Davide Iafrate
%
% VERSIONS
%   2020-11-20: First version


TOF = t2 - t1;         % [s]
%% Call Lambert solver

[~,~,~,~,v1L,v2L,~,~] =...
    lambertMR(ri,rf,TOF,mu,0,0,0,1);

%% Compute the total cost of the maneuvre
DV1 = norm(vi-v1L);
DV2 = norm(v2L-vf);

DV = DV1 + DV2;

T1 = t1;
T2 = t2;

r_tf0 = ri;
v_tf0 = v1L;