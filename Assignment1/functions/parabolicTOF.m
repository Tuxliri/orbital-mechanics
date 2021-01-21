function t_p = parabolicTOF(r1,r2,mu)
% Quick function to calculate the parabolic ToF between two planets, in the
% case of Hohmann transfer
%
% PROTOTYPE:
%    t_p = parabolicTOF(r1,r2,mu)
% 
%  INPUT :
%       mu[1]   planetary constant of the central   [km^3/s^2]
%       r1[1]   orbital radius of first planet      [km]
%       r2[1]   orbital radius of the second planet [km]
% 
%   OUTPUT:
%       t_p[1]  parabolic time of flight            [years]
% 
%   CONTRIBUTORS:
%       Davide Iafrate
% 
%   REVISIONS:
%       12/01/2021  First version

% Conversion constant
SECS2YEARS = 365*24*3600;

c = abs(r1 - r2);
s = (c + r1 + r2)/2;

t_p = sqrt(2)/3*sqrt(s^3/mu)*(1 - ((s - c)/s)^(3/2));

t_p = t_p/SECS2YEARS;