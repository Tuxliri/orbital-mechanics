function ds = EOM_RSW(t,s,mu,a_per_rsw)

% EOM_RSW Equation of motion that returns the derivative of the keplerian
% elements at a given time, including the perturbing accelerations in the
% RSW (Radial-transveral-out of plane) reference frame
% 
% PROTOTYPE:
%   ds = EOM(t,s,mu,a_per)
% 
% INPUT:
%   t[1]        time of evaluation
%   s[6]        state vector containing the keplerian elements
%               [a,e,i,RAAN,omega,f]            [km,-,rad,rad,rad,rad]
%   mu[1]       gravitational parameter of main body
%   a_per_rsw[3]   function handle that returns the vector of perturbing
%               accelerations [ar,as,aw]
% 
% OUTPUT:
%   ds[6]   derivative of the state vector
% 
% CONTRIBUTORS:
%   Davide Iafrate      14-12-2020

% Extract the accelerations in the rsw frame
acc = a_per_rsw(t,s);
ar = acc(1);
as = acc(2);
aw = acc(3);

% Extract the keplerian elements from the state vector
a = s(1);
e = s(2);
i = s(3);
RAAN = s(4);
omega = s(5);
f = s(6);

% Calculate the needed orbital values
p = a*(1 - e^2);
r = p / (1+e*cos(f));
v = sqrt(2*mu/r - mu/a);
h = sqrt(p*mu);

% Compute the derivatives of the state
da = 2*a^2/h * (e*sin(f)*ar + p*as/r);

de = 1/h * (p*sin(f)*ar + ((p+r)*cos(f) + r*e)*as);

di = r*cos(f+omega) *aw/h;

dRAAN = r*sin(f+omega)*aw/(h*sin(i));

domega = 1/(h*e) *( - p*cos(f)*ar + (p+r)*sin(f)*as) ...
    - r*sin(f+omega)*cos(i)*aw/(h*sin(i));

df = h/r^2 + 1/(e*h) * (p*cos(f)*ar - (p+r)*sin(f)*as);

% Vector of derivatives of the state
ds = [da; de; di; dRAAN; domega; df];