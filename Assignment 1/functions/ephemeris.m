function [kep_el,R,V] = ephemeris(timevecMJD2000,PLANETNUM)
%EPHEMERIS Calculate the keplerian elements of a planet
%          for a given time vector (MJD2000 format)
%   INPUT:
%       timevecMJD2000[N]   vector of times for which to calculate the
%                           ephemeris       [days]
%       PLANETNUM[1]        ID of the planet (1-9) [-]
%       
%   OUTPUT:
%       kep_el[N,6]         vector of keplerian elements
%                           [a,e,i,RAAN.omega,theta]

%       R[N,3]              position vectors of planetary body
%                           [km]

%       V[N,3]              velocities of planetary body [km/s]
%       mu[1]               gravitational parameter of the central body
%                           [km^3/s^2]

% PREALLOCATION for speed
N = length(timevecMJD2000);
kep_el = zeros(N,6);
R = zeros(N,3);
V = zeros(N,3);

for i = 1:N
    [kep_el(i,:),mu] = uplanet(timevecMJD2000(i),PLANETNUM);
    [R(i,:),V(i,:)] = kep2car(kep_el(i,1),kep_el(i,2),kep_el(i,3),kep_el(i,4),kep_el(i,5),kep_el(i,6),mu);
end

