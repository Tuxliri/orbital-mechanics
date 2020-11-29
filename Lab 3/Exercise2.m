% Exercise 1

% Set path to default
path(pathdef);
% Add [...] folder to path
addpath(genpath('functions\'));

addpath(genpath('../Common\'));

%% Physical parameters
muE = astroConstants(13);

% Input
a1 = 12500;         % [km]
e1 = 0;             % [-]
i1 = 0;             % [deg]
RAAN1 = 0;          % [deg]
omega1 = 0;         % [deg]
f1 = deg2rad(120);  % [rad]
% T1 = 2*pi*sqrt(muE/
[r1,v1] = kep2car(a1,e1,i1,RAAN1,omega1,f1,muE);


% Input
a2 = 9500;          % [km]
e2 = 0.3;           % [-]
i2 = 0;             % [deg]
RAAN2 = 0;          % [deg]
omega2 = 0;         % [deg]
f2 = deg2rad(250);  % [rad]

[r2,v2] = kep2car(a2,e2,i2,RAAN2,omega2,f2,muE);

TOF = 3300;         % [s]
%% Call Lambert solver

[A,P,E,ERROR,v1L,v2L,TPAR,THETA] =...
    lambertMR(r1,r2,TOF,muE,0,0,0,1);

%% Compute the total cost of the maneuvre

DV = norm(v1-v1L) + norm(v2L-v2)

% (r1, v1,t1)
% (r2, v2, t1+Dt)
tspan = (1:1:TOF);    % Use the period instead [IMPROVEMENT]


[r,v] = propagator(r1,v1L,muE,tspan);
% [rVIRTUAL,vVIRTUAL] =  = propagator(r1,v1L,muE,tspan);
PLOTORBIT(r);
