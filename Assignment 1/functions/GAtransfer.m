function [DV,DV_dep,DV_arr,DV_ga,FLYBY,TRANSFER1,TRANSFER2]= GAtransfer(planetA,planetB,planetC,t_dep,t_ga,t_arr)
% Get the different DV needed for the transfer between two celestial bodies
% and the respective times of departure and arrival
%
% PROTOTYPE:
%    [DV, DV_dep, DV_arr, DV_ga, FLYBY,TRANSFER1,TRANSFER2]= GAtransfer(planetA, planetB, planetC, t_dep, t_ga, t_arr)
%
%  INPUT :
%	planetA [ 1 ] Integer number identifying the departure celestial body (< 11)
%                   1:   Mercury
%                   2:   Venus
%                   3:   Earth
%                   4:   Mars
%                   5:   Jupiter
%                   6:   Saturn
%                   7:   Uranus
%                   8:   Neptune
%                   9:   Pluto
%                   10:  Sun
%	planetB [ 1 ] Integer number identifying the flyby celestial body (< 11)
%	planetC [ 1 ] Integer number identifying the arrival celestial body (< 11)
%   t_dep         Departure time in mjd2000
%   t_ga          Gravity assist time in mjd2000
%   t_arr         Arrival time in mjd2000
% 
%
% OUTPUT:
%   DV      [ 1 ] Total delta V for transfer      [ km/s ]
%   DV_dep  [ 1 ] Cost of departure maneuvre      [ km/s ]
%   DV_arr  [ 1 ] Cost of arrival maneuvre        [ km/s ]
%   DV_ga   [ 1 ] Cost of gravity assist maneuvre [ km/s ]
%
% CALLED FUNCTIONS:
%   astroConstants
%   ephemeris
%   single_arc
%   flybypow
%
% CONTRIBUTORS:
%   Alkady Marwan
%   Bossi Nunez Pedro
%   Demartini Davide
%   Davide Iafrate
%
% VERSIONS
%   2020-11-23: First version


muSun = astroConstants(4);

% Compute the positions and velocities of the planet with the ephemeris
[~, r_A, v_A] = ephemeris(t_dep, planetA.ID);
[~, r_B, v_B] = ephemeris(t_ga, planetB.ID);   % Orbital radius of the GA planet
[~, r_C, v_C] = ephemeris(t_arr, planetC.ID);

% Calculate the first transfer arc and the cost in deltav to depart from 
% the initial orbit of the planet to the transfer orbit
% the velocity at arrival of planet B is the one to be used as the V_minus
% in the computation of the powered flyby manoeuvre (v_GA_minus)
[~, DV_dep , ~, TRANSFER1.v1, v_GA_minus] = single_arc(t_dep, t_ga,muSun, r_A, r_B, v_A, v_B);
TRANSFER1.r1 = r_A;
TRANSFER1.r2 = r_B;
TRANSFER1.v2 = v_GA_minus;
TRANSFER1.t_dep = t_dep;
TRANSFER1.t_arr = t_ga;

% Calculate the second transfer arc and the cost to get from the second 
% transfer orbit to the orbit of the arrival planet
[~, ~, DV_arr, v_GA_plus, TRANSFER2.v2] = single_arc(t_ga, t_arr, muSun, r_B, r_C, v_B, v_C);
TRANSFER2.r1 = r_B;
TRANSFER2.r2 = r_C;
TRANSFER2.v1 = v_GA_plus;
TRANSFER2.t_dep = t_ga;
TRANSFER2.t_arr = t_arr;

% Calcualte the cost of the gravity assist
[DV_ga, delta, a, e, Delta, rp, vp0] = flybypow(v_GA_minus - v_B, v_GA_plus - v_B,planetB.mu);
DV = DV_dep + DV_arr + abs(DV_ga);
FLYBY.DV = DV_ga;
FLYBY.delta = delta;
FLYBY.a = a;
FLYBY.e = e;
FLYBY.Delta = Delta;
FLYBY.rp = rp;
FLYBY.vp0 = vp0;

end