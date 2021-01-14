function [a_RSW, a_ECI] = a_SRP(s,Cr,Psr,R_E,Am,mjd2000)
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
%   a_RSW[3x1]        vector of perturbing accelerations [ar,as,aw]       [km/s^2]
%   a_[3x1]        vector of perturbing accelerations [ar,as,aw]       [km/s^2]
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

% Get keplerian elements of the Earth in Sun-centred ecliptic 
%     system

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

% Vector from Earth to Sun in sun-centered ecliptic frame
rS = -rE;

ROT = [1 0 0;
       0 cos(epsilon) -sin(epsilon);
       0 sin(epsilon) cos(epsilon);];
   
% Vector in ECI coordinates

rS = ROT*(rS');

% Get position vector R-sc-sun
% R_sc_sun = -(rE + rSC');

%% IMPROVEMENT/NEEDED: IMPLEMENT ECLIPSE CONDITION!!!

% ECLIPSE CHECK algorithm from curtis page 526
% Vector from Earth to sun is opposite of vecotr from Sun to Earth
% rS = -rE;

theta = acos(dot(rS,rSC)/(norm(rS)*norm(rSC)));
theta_1 = acos(R_E/norm(rSC));
theta_2 = acos(R_E/norm(rS));
if theta_1+theta_2 < theta
    vi = 0;
else
    vi=1;
end

%  Calculate the vector of acceleration in the RSW frame AND CONVERT FROM
%  m/s to km/s for the propagators
% a_ECI = - vi*Psr * AU^2 / norm(R_sc_sun)^3 * Cr * Am * R_sc_sun /1000;
a_ECI = - vi*Psr * AU^2 / norm(rS)^3 * Cr * Am * rS /1000;

% Convert to the RSW reference frame

r_vers = rSC/norm(rSC);
w_vers = cross(rSC,vSC)./norm(cross(rSC,vSC));
s_vers = cross(w_vers,r_vers);

a_RSW = [r_vers' s_vers' w_vers']' * a_ECI;

end


