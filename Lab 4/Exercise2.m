clc
clearvars
close all

%% Set path to default
path(pathdef);
% Add [...] folder to path
addpath(genpath('../functions\'));

addpath(genpath('../Common\'));

%% useful constants
muSun = astroConstants(4);
muE = astroConstants(13);
AU = astroConstants(2);

% Define Data
rE = [0; -1; 0] * AU;
radius_e = astroConstants(23);      % Earth mean radius
Vm = [31.5; 4.69; 0];
Vp = [38.58; 0; 0];

% Compute planet velocity
VEmag = sqrt(muSun/norm(rE));

VE = VEmag .* cross([0;0;1],rE)/norm(rE);

% Compute velocities relative to the planet
vinfM = Vm - VE;                % Entry velocity in SOI
vinfP = Vp - VE;                % Exit velocity in SOI

% Solve powered flyby
h_min = 400;                    % Minimum height above the planet [km]
[dv_p,delta,Delta,rp] = flybypow(vinfM,vinfP,muE,radius_e+h_min);

% Compute useful data
deltadeg = rad2deg(delta);
h_ga = rp - radius_e;
% compute the total delta v

% Plot the two planetocentric hyperbolic arcs