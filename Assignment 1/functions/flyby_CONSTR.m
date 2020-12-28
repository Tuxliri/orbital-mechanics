function [c,ceq] = flyby_CONSTR(planetA,planetB,planetC,t_dep,t_ga,t_arr, rp_min)
%FLYBY_CONSTR Constraint function for the genetic algorithm function
%             FIRST CONSTRAINT: rp_min < rp     the radius of the flyby
%             must be higher than the minimum radius to avoid crashing into
%             the planet atmosphere
[~, ~, ~,FLYBY]...
    = GAtransfer(planetA,planetB,planetC,t_dep,t_ga,t_arr);

% c(1) = DV1 - max_vinf;
% c(2) = FLYBY.DV + DV2 - max_vinf_sc;
c(1) = rp_min - FLYBY.rp;
ceq = [];

end

