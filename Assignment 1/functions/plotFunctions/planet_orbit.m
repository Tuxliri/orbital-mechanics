function [kep_el,R,V,R_dep,R_GA,R_arr] = planet_orbit(planetID,OPT_DEP,OPT_GA,OPT_ARR,muSun)
% PLANET ORBIT calculate orbit of a planet

% Conversion constant
DAY2SECS = 24*60*60;

[kepA,~,~] = ephemeris(OPT_DEP,planetID);
T_A = 2*pi*sqrt(kepA(1)^3 / muSun);              % orbital period of planetA [s]

tspan = linspace(OPT_DEP,OPT_DEP + T_A/DAY2SECS,1000);

[kep_el,R,V] = ephemeris(tspan,planetID);
[~,R_dep,~] = ephemeris(OPT_DEP,planetID);
[~,R_GA,~] = ephemeris(OPT_GA,planetID);
[~,R_arr,~] = ephemeris(OPT_ARR,planetID);

end

