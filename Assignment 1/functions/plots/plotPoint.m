function plotPoint(a,e,i,RAAN,om,theta,mu,val)
% Plots position given keplerian elements
%
% PROTOTYPE:
%    plotPoint(a, e, i, RAAN, om, theta, mu, val) 
%   
% INPUT:
%   a[1]        semi-major axis     [ km ]
%   e[1]        eccentricity        [ - ]
%   i[1]        inclination         [ rad ]
%   RAAN[1]     right ascension of the ascending node [ rad ]
%   omega[1]    argument of perigee [ rad ]
%   f[1]        true anomaly        [ rad ]
%   mu[1]      gravitational paramer [ km^3/s^2 ]
%
% OUTPUT:
%   
% CONTRIBUTORS:
%   Alkady Marwan
%   Bossi Nunez Pedro
%   Davide Demartini
%   Davide Iafrate
%
% VERSIONS
%   2020-05-15: First version
%   2021-01-12: Second version
%
% See kep2car instructions for help

% Input check
if nargin<8
    val='ko'
end
if nargin<7
    mu=398600;
    val='ko';
end
if i==0 && e~=0
    RAAN=0;
elseif e==0 && i~=0
    om=0;
elseif i==0 && e==0
    om=0;
    RAAN=0;
end

% Rotation Matrixes
R_Ohm = [cos(RAAN) sin(RAAN) 0
    -sin(RAAN) cos(RAAN) 0
    0 0 1];
R_i = [1 0 0 
    0 cos(i) sin(i)
    0 -sin(i) cos(i)];
R_omega=[cos(om) sin(om) 0
    -sin(om) cos(om) 0
    0 0 1];
T=R_Ohm'*R_i'*R_omega';

% p parameter
p=a*(1-e^2);
h=sqrt(p*mu);

% Position in perifocal frame
xpf=h^2/mu*(1/(1+e*cos(theta)))*cos(theta);
ypf=h^2/mu*(1/(1+e*cos(theta)))*sin(theta);

% Reference frame switch
rge=T*[xpf;ypf;0];

% Computation of position vectors 
X=rge(1);
Y=rge(2);
Z=rge(3);

% Point representation
plot3(X,Y,Z,val,'LineWidth',1.5);
end
