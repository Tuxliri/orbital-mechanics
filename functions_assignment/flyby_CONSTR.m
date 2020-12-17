function [c,ceq] = flyby_CONSTR(planetA,planetB,planetC,t_dep,t_ga,t_arr, max_vinf, max_vinf_sc,rp_min)
%FLYBY_CONSTR Constraint function for the genetic algorithm function
%             FIRST CONSTRAINT: DV_dep < max_vinf      the deltav for the 
%                               departure must be less than the maximum
%                               the launch vehicle is capable of
% 
%             SECOND CONSTRAINT: DV_GA + DV_ARR < max_vinf_sc   the combined
%                               deltav for the perigee burn at the flyby
%                               and the one for the final orbit insertion
%                               is in the capabilities of the probe
%             THIRD CONSTRAINT: rp_min < rp     the radius of the flyby
%             must be higher than the minimum radius to avoid crashing into
%             the planet atmosphere
[~, DV1, DV2, DV_ga, rp]...
    = GAtransfer(planetA,planetB,planetC,t_dep,t_ga,t_arr);

c(1) = DV1 - max_vinf;
c(2) = DV_ga + DV2 - max_vinf_sc;
c(3) = rp_min - rp;
ceq = [];

end

