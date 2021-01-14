% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function [V2, ierr] = gibbs(R1, R2, R3)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
  This function uses the Gibbs method of orbit determination to
  to compute the velocity corresponding to the second of three
  supplied position vectors.
 
  mu            - gravitational parameter (km^3/s^2 
  R1, R2, R3    - three coplanar geocentric position vectors (km)
  r1, r2, r3    - the magnitudes of R1, R2 and R3 (km)
  c12, c23, c31 - three independent cross products among
                  R1, R2 and R3
  N, D, S       - vectors formed from R1, R2 and R3 during
                  the Gibbs' procedure
  tol           - tolerance for determining if R1, R2 and R3
                  are coplanar
  ierr          - = 0 if R1, R2, R3 are found to be coplanar
                  = 1 otherwise
  V2            - the velocity corresponding to R2 (km/s)
 
  User M-functions required: none
%}
% ---------------------------------------
 
global mu
tol  = 1e-4;
ierr = 0;
 
%...Magnitudes of R1, R2 and R3:
r1 = norm(R1);
r2 = norm(R2);
r3 = norm(R3);
 
%...Cross products among R1, R2 and R3:
c12 = cross(R1,R2);
c23 = cross(R2,R3);
c31 = cross(R3,R1);
 
%...Check that R1, R2 and R3 are coplanar; if not set error flag:
if abs(dot(R1,c23)/r1/norm(c23)) > tol
    ierr = 1;
end
 
%...Equation 5.13:
N = r1*c23 + r2*c31 + r3*c12;
 
%...Equation 5.14:
D = c12 + c23 + c31;
 
%...Equation 5.21:
S = R1*(r2 - r3) + R2*(r3 - r1) + R3*(r1 - r2);
 
%...Equation 5.22:
V2 = sqrt(mu/norm(N)/norm(D))*(cross(D,R2)/r2 + S);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  end %gibbs

[End MATLAB script]

Script file: Example_5_01.m
[Begin MATLAB script]
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example_5_01
% ~~~~~~~~~~~~
%{
  This program uses Algorithm 5.1 (Gibbs method) and Algorithm 4.2
  to obtain the orbital elements from the data provided in Example 5.1.
 
  deg        - factor for converting between degrees and radians
  pi         - 3.1415926...
  mu         - gravitational parameter (km^3/s^2)
  r1, r2, r3 - three coplanar geocentric position vectors (km)
  ierr       - 0 if r1, r2, r3 are found to be coplanar
               1 otherwise
  v2         - the velocity corresponding to r2 (km/s)
  coe        - orbital elements [h e RA incl w TA a]
               where h    = angular momentum (km^2/s)
                     e    = eccentricity
                     RA   = right ascension of the ascending node (rad)
                     incl = orbit inclination (rad)
                     w    = argument of perigee (rad)
                     TA   = true anomaly (rad)
                     a    = semimajor axis (km)
  T          - period of elliptic orbit (s)
 
  User M-functions required: gibbs, coe_from_sv
%}
% ----------------------------------------------
 
clear all; clc
deg = pi/180;
global mu
 
%...Data declaration for Example 5.1:
mu = 398600;
r1 = [-294.32 4265.1 5986.7];
r2 = [-1365.5 3637.6 6346.8];
r3 = [-2940.3 2473.7 6555.8];
%...
 
%...Echo the input data to the command window:
fprintf('-----------------------------------------------------')
fprintf('\n Example 5.1: Gibbs Method\n')
fprintf('\n\n Input data:\n')
fprintf('\n  Gravitational parameter (km^3/s^2)  = %g\n', mu)
fprintf('\n  r1 (km) = [%g  %g  %g]', r1(1), r1(2), r1(3))
fprintf('\n  r2 (km) = [%g  %g  %g]', r2(1), r2(2), r2(3))
fprintf('\n  r3 (km) = [%g  %g  %g]', r3(1), r3(2), r3(3))
fprintf('\n\n');
 
%...Algorithm 5.1:
[v2, ierr] = gibbs(r1, r2, r3);
 
%...If the vectors r1, r2, r3, are not coplanar, abort:
if ierr == 1
    fprintf('\n  These vectors are not coplanar.\n\n')
    return
end
 
%...Algorithm 4.2:
coe  = coe_from_sv(r2,v2,mu);
 
h    = coe(1);
e    = coe(2);
RA   = coe(3);
incl = coe(4);
w    = coe(5);
TA   = coe(6);
a    = coe(7);
 
%...Output the results to the command window:
fprintf(' Solution:')
fprintf('\n');
fprintf('\n  v2 (km/s) = [%g  %g  %g]', v2(1), v2(2), v2(3))
fprintf('\n\n  Orbital elements:');
fprintf('\n    Angular momentum (km^2/s)  = %g', h) 
fprintf('\n    Eccentricity               = %g', e)
fprintf('\n    Inclination (deg)          = %g', incl/deg)
fprintf('\n    RA of ascending node (deg) = %g', RA/deg)
fprintf('\n    Argument of perigee (deg)  = %g', w/deg)
fprintf('\n    True anomaly (deg)         = %g', TA/deg) 
fprintf('\n    Semimajor axis (km)        = %g', a)
%...If the orbit is an ellipse, output the period:
if e < 1
    T = 2*pi/sqrt(mu)*coe(7)^1.5;
    fprintf('\n    Period (s)                 = %g', T)
end 
fprintf('\n-----------------------------------------------------\n')
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

