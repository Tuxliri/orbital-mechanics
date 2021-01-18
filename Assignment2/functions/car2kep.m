function [a, e, i, RAAN, omega, f] = car2kep(r_vec,v_vec,mu)
% Convert from state vector in cartesian coordinates to keplerian coordinates
% SPECIAL CASES: 
%
%       [i = 0]
%           the NODE LINE is taken as the reference direction Gamma (vernal
%           equinox
%           RAAN is 0
%           argument of periapsis omega is undefined (NaN)
%       
%       [e = 0 & i != 0]
%           f is calculated as the argument of latitude, angular
%           distance between the line of nodes and the position vector
%           omega is 0
%       [e = 0 & i = 0]
%           f is calculated as the true longitude, angular
%           distance between the gamma direction and the position vector
% PROTOTYPE:
%   [a, e, i, RAAN, omega, f] = car2kep(r_vec,v_vec,mu)
%   
% INPUT:
%
%   r_vec[3]   position vector     [ km/s ]
%   v_vec[3]   velocity vector     [ km/s ]
%   mu[1]      gravitational paramer [ km^3/s^2 ]
%
% OUTPUT:
%
%   a[1]        semi-major axis     [ km ]
%   e[1]        eccentricity        [ - ]
%   i[1]        inclination         [ rad ]
%   RAAN[1]     right ascension of the ascending node [ rad ]
%   omega[1]    argument of perigee [ rad ]
%   f[1]        true anomaly        [ rad ]
%
% CONTRIBUTORS:
%   Davide Iafrate      
%   Alkady Marwan
%   Pedro Bossi Núñez
%   Davide Demartini
%
% VERSIONS
%   2020-10-12: First version

eps = 1.e-10;

% calculate magnitude of r_vec, v_vec

r = norm(r_vec);
v = norm(v_vec);

%% compute the radial component of velocity

v_radial = dot(v_vec,r_vec)/r;

%% compute the angular momentum vector

h_vec = cross(r_vec,v_vec);

%% compute the magnitude of the angular momentum vector

h = norm(h_vec);

%% compute the inclination of the orbit
i = acos(h_vec(3)/h);

%% compute the node line

N_vec = cross([0 0 1], h_vec);

n = norm(N_vec);

%% compute the RAAN     % solving the case i = 0

if i < eps
    RAAN = 0;
    
elseif (N_vec(2)>= 0)
    RAAN = acos(N_vec(1)/n);
    
else
    RAAN = 2*pi - acos(N_vec(1)/n);

end

%% compute the eccentricity vector and its magnitude

e_vec = (r_vec.*(v^2 - mu/r) - r.*v_radial.*v_vec)./mu;
e = norm(e_vec);


%% compute the argument of perigee omega

if e < eps
    
    omega = 0;

elseif (e_vec(3) >= 0)
    omega = acos(dot(N_vec,e_vec)/(n*e));

else
    omega = 2*pi - acos(dot(N_vec,e_vec)/(n*e));
    
end



%% compute the true anomaly f

% REFERENCE:
% https://en.wikipedia.org/wiki/True_anomaly
% used for special cases

if e < eps     
    if i < eps               % Output the true longitude, omega is set to 0:
        f = acos(r_vec(1) / r);
        
        if v_vec(1) > 0
            
            f = 2*pi - acos(r_vec(1) / r);
        
        end
        
    else
        
        % Output the argument of latitude
        % https://en.wikipedia.org/wiki/True_anomaly#Circular_orbit
        % omega is assumed as 0
        if r_vec(3) < 0
            
            f = 2*pi - acos(dot(N_vec , r_vec) / (n*r));
            
        else
            f = acos(dot(N_vec , r_vec) / (n*r));
        end
    end
elseif (v_radial >= 0)
    f = acos(dot(e_vec,r_vec)/(e*r));
else
    f = 2*pi - acos(dot(e_vec,r_vec)/(e*r));
end

%% semi-major axis (scattered among various lectures)

if e < eps % Circular orbits
    a = r;
elseif e == 1 % Parabolic orbit
    a = 'infinite';
else % Elliptical/Hyperbolic orbit
    a = (h^2/mu)/(1 - e^2);
end

