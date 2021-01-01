function a = a_SRP_RSW(s,Cr,Psr,R_E,Am,mjd2000)
%a_SRP_RSW this function calculates the acceleration in the RSW
%   (Radial-transveral-out of plane) reference frame, returning a vector of
%   the three components
%
% PROTOTYPE:
%   a = a_j2_RSW(mu,s,J2,R_E)
%
% INPUT:
%   s[6]        state vector containing the keplerian elements
%               [a,e,i,RAAN,omega,f]                 [km,-,rad,rad,rad,rad]
%   cr[1]       Refelctivity coefficient       [-]
%   Psr[1]      Solar radiation pressure               [N/m^2]

%   R_E[1]      mean radius of the Earth            [km]
%   Am[1]       area to mass ratio of the S/C       [m^2/kg]
%   mjd2000[1]  mjd2000 date                        [-]
%
% OUTPUT:
%   a[3]        vector of perturbing accelerations [ar,as,aw]       [km/s^2]
%
% CONTRIBUTORS:
%   Davide Iafrate      14-12-2020

% Extract the keplerian elements from the state vector
a = s(1);
e = s(2);
i = s(3);
RAAN = s(4);
omega = s(5);
f = s(6);

% Calculate the radius of the S/C wrt to sun in the Inertial frame:
% VECTORIAL SUM: R_sc-sun = - (R_sun-e + R_e-sc

% Get keplerian elements of the Earth
[kepSun,muSun] = uplanet(mjd2000, 3); 
muE = astroConstants(13);
AU = astroConstants(2);

% get position vector R_e-sun IN ECLIPTIC FRAME
[rE,~] = kep2car(kepSun(1),kepSun(2),kepSun(3),kepSun(4),kepSun(5),kepSun(6),muSun);

% get position vector R_sc-e IN ECI FRAME
[rSC,vSC] = kep2car(a,e,i,RAAN,omega,f,muE);

% To get the position vector of the Earth wrt the sun in the ECI frame we
% must use a rotation matrix to rotate around the X axis by an amount equal
% to the negative inclination of the spin axis of the Earth

epsilon = astroConstants(8);
ROT = [cos(epsilon) sin(epsilon) 0;
                  -sin(epsilon) cos(epsilon) 0;
                   0    0   1;];
rE = ROT*rE;

% Get position vector R-sc-sun
R_sc_sun = - (rE + rSC);

%% IMPROVEMENT/NEEDED: IMPLEMENT ECLIPSE CONDITION!!!

%  Calculate the vector of acceleration in the RSW frame
a_SRP_ECI = - Psr * AU^2 / norm(R_sc_sun)^3 * Cr * Am * R_sc_sun;

% Convert to the RSW reference frame

r_vers = rSC/norm(rSC);
w_vers = cross(rSC,vSC)./norm(cross(rSC,vSC));
s_vers = cross(w_vers,r_vers);

a = [r_vers' s_vers' w_vers']' * a_SRP_ECI;

end


