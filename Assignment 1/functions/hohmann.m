function [DV1,DV2,TOF] = hohmann(ID_A,ID_B,date)
% Calculator for hohmann transfers between planets
% 
% ATTENTION:
%   this program works with circular orbits. If the orbit of interest is
%   elliptical, the program approximates such orbit with a circular one.
%   The radiuses of the circular orbits are the distance between the Sun and
%   the planets calculated on given date
%
% PROTOTYPE:
%   [DV1, DV2, TOF] = hohmann(ID_A, ID_B, date)
% 
% INPUT:
%   ID_A [ 1 ]     id of departure planet
%   ID_B [ 1 ]     id of arrival planet
%   Date [ 6 x 1 ] date of interest
% 
% OUTPUT: 
%   DV1  [ 1 ]     delta-v for departure
%   DV2  [ 1 ]     delta-v for arrival
%   TOF  [ 1 ]     time of flight for the transfer
%
% CALLED FUNCTIONS:
%   date2mjd2000
%   uplanet
%   kep2ca
%
% CONTRIBUTORS:
%   Alkady Marwan
%   Bossi Núñez Pedro
%   Demartini Davide
%   Davide Iafrate
%
% VERSIONS
%   2020-12-13: First version

% Determining planets positions in given date
datemjd2000 = date2mjd2000(date);
[kepA ,~] = uplanet(datemjd2000, ID_A);
[kepB, muSun] = uplanet(datemjd2000, ID_B);

% Fictitious circular orbits
[rA, vA] = kep2car(kepA(1), kepA(2), kepA(3), kepA(4), kepA(5), kepA(6), muSun);
[rB, vB] = kep2car(kepB(1), kepB(2), kepB(3), kepB(4), kepB(5), kepB(6), muSun);

% Hohmann orbit characterization
a_hoh = 0.5*(norm(rA) + norm(rB));
if norm(rA) < norm(rB)
    rp_hoh = norm(rA);
    ra_hoh = norm(rB);
    v_p = norm(vA);
    v_a = norm(vB);
else
    rp_hoh = norm(rB);
    ra_hoh = norm(rA);
    v_p = norm(vB);
    v_a = norm(vA);
end

va_hoh = sqrt(muSun*(2/ra_hoh - 1/a_hoh));
vp_hoh = sqrt(muSun*(2/rp_hoh - 1/a_hoh));

DV1 = abs(va_hoh - v_a);
DV2 = abs(vp_hoh - v_p);

TOF = pi*sqrt(a_hoh^3 / muSun)/(365*24*3600);eb