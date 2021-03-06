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
[dv_p,delta,a,e,Delta,rp,vp0] = flybypow(vinfM,vinfP,muE,radius_e+h_min);

% Compute useful data
deltadeg = rad2deg(delta);
h_ga = rp - radius_e;
% compute the total delta v

% Plot the two planetocentric hyperbolic arcs
tspanA = linspace(0,-10000,10000);         % Plotting backward in time
[r0,v0A] = kep2car(a(1),e(1),0, 0,0,0,muE);
[r_A,v_A] = propagator(r0,v0A,muE,tspanA);

tspanB = linspace(0,10000,10000);
[r0,v0B] = kep2car(a(2),e(2),0,0,0,0,muE);
[r_B,v_B] = propagator(r0,v0B,muE,tspanB);


% Plot the arrival hyperbola
figure(1);
hold on
legend on
ARR = plot3(r_A(:,1),r_A(:,2),r_A(:,3), 'LineWidth', 2);
ARR.DisplayName = 'Incoming hyperbola';
grid on
ylabel('y [km]')
xlabel('x [km]')

% Plot the departure hyperbola
DEP = plot3(r_B(:,1),r_B(:,2),r_B(:,3), 'LineWidth', 2);
DEP.DisplayName = 'Outcoming hyperbola';

% Plot the asymptotes
b_dep = a(2)*tan(acos(1/e(2)));
b_arr = a(1)*tan(acos(1/e(1)));

x2 = linspace(r0(1)-a(2),r_A(end,1),10000);

y_dep = -b_dep/a(2) .* (x2 - r0(1) + a(2));

x1 = linspace(r0(1)-a(1),r_A(end,1),10000);
y_arr = b_arr/a(1) .* (x1 - r0(1) + a(1));

ARR_asymptote = plot3(x1,y_arr,zeros(1,length(x1)),'LineWidth',2);
DEP_asymptote = plot3(x2,y_dep,zeros(1,length(x2)),'LineWidth',2);