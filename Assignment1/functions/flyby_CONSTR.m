function [c,ceq] = flyby_CONSTR(planetA,planetB,planetC,t_dep,t_ga,t_arr,rp_min)
% FLYBY_CONSTR Constraint function for the genetic algorithm function
%             FIRST CONSTRAINT: rp_min < rp the radius of the flyby
%             must be higher than the minimum radius to avoid crashing into
%             the planet's atmosphere
%
% PROTOTYPE:
%   [c, ceq] = flyby_CONSTR(planetA, planetB, planetC, t_dep, t_ga, t_arr, rp_min)
%
%   INPUT:
%       PlanetA [ N ]     
%       PlanetB [ 1 ]
%       PlanetC [   ]
%       t_dep
%       t_arr
%       rp_min
%       
%   OUTPUT:
%    see documentation for fmincon - nonlinear
%       c
%       c_eq
%
% CONTRIBUTORS:
%   Davide Demartini
%   Davide Iafrate
%   Marwan Alkady
%   Pedro Bossi Núñez
%
% VERSIONS
%   2020-12-26: First version
%
% CALLED FUNCTIONS:
%   GAtransfer

[~, ~, ~, ~,FLYBY]= GAtransfer(planetA,planetB,planetC,t_dep,t_ga,t_arr);

% needed in the fmincon
c(1) = rp_min - FLYBY.rp;
ceq = [];

end