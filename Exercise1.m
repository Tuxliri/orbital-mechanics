% Exercise 1

% Set path to default
path(pathdef);
% Add [...] folder to path
addpath(genpath('functions\'));

addpath(genpath('../Common\'));

%% Physical parameters
muE = astroConstants(13);

% Input
r1 = [-21800, 37900, 0];    % [km]
r2 = [27300, 27700, 0];     % [km]
Dt = 15*3600 + 6*60 + 40;

%% Call Lambert solver

[A,P,E,ERROR,v1,v2,TPAR,THETA] =...
    lambertMR(r1,r2,Dt,muE,0,0,0,1);

%%
% (r1, v1,t1)
% (r2, v2, t1+Dt)
tspan = (1:1:80000);    % Use the period instead [IMPROVEMENT]

[r,v] = propagator(r1,v1,muE,tspan);

PLOTORBIT(r);

