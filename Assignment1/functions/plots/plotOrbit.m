function plotOrbit(a,e,i,RAAN,om,mu,color)
% Plots planet's orbit
%
% PROTOTYPE:
%    plotOrbit(a, e, i, RAAN, om, mu, color) or plotOrbit(obj, date)
%   
% INPUT:
%
%   a[1]        semi-major axis     [ km ]
%   e[1]        eccentricity        [ - ]
%   i[1]        inclination         [ rad ]
%   RAAN[1]     right ascension of the ascending node [ rad ]
%   omega[1]    argument of perigee [ rad ]
%   f[1]        true anomaly        [ rad ]
%   mu[1]       gravitational paramer [ km^3/s^2 ]
%   color[1]    display attribute (see plot3) [ - ]
%
%   obj[1]      object of interest [ - ]
%   date[6]     date of interest in Gregorian calendar [ yyyy mm dd hh min ss ]
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
%   2021-01-17: Third version
%
% CALLED FUNCTIONS:
%   ephemeris
%
% See kep2car & Plot_Celestial_Object instructions for help

if nargin == 2 % Second case. Celestial bodies only!
    color = '-';
    
    % Rename vectors
    obj = a;
    date = e;
    
    % Periods list [ days ]
    T = [88 224 365 686.8 4332.55 10749.25 30685.55 60188.5];
    
    % Date of interest
    t0 = date2mjd2000(date);
    t = linspace((t0 - T(obj)/2*1.01), (t0 + T(obj)/2*1.01), 2000*obj);
    
    % Ephemeris determination
    [~, R, ~] = ephemeris(t, obj);
    X = R(:,1); Y = R(:,2); Z = R(:,3);
    
else % First case
    
    % Input check
    if nargin < 6
            mu = 398600;
    end
    if i==0 && e~=0
        RAAN=0;
    elseif e==0 && i~=0
        om=0;
    elseif i==0 && e==0
        om=0;
        RAAN=0;
    end

    % Rotation matrixes

    R_Ohm = [cos(RAAN) sin(RAAN) 0
        -sin(RAAN) cos(RAAN) 0
        0 0 1]; 
    R_i = [1 0 0 
        0 cos(i) sin(i)
        0 -sin(i) cos(i)];
    R_omega = [cos(om) sin(om) 0
        -sin(om) cos(om) 0
        0 0 1];
    T = R_Ohm'*R_i'*R_omega';

    % p parameter
    p = a*(1 - e^2);
    h = sqrt(p*mu);

    % Position vectors
    X = []; Y = []; Z = [];

    for theta=0:0.01:2*pi
        % Position computation in perifocal reference frame
        xpf = h^2/mu*(1/(1 + e*cos(theta)))*cos(theta);
        ypf = h^2/mu*(1/(1 + e*cos(theta)))*sin(theta);
        % Switching reference frame
        rge = T*[xpf; ypf; 0];
        % Position vectors
        x = rge(1);
        y = rge(2);
        z = rge(3);
        X = [X; x];
        Y = [Y; y];
        Z = [Z; z];
    end
end

% Data representation
plot3(X, Y, Z, color, 'LineWidth', 2);
end
