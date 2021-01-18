% Test for the assignment 1 using the genetic algorithm
clc
clearvars
close all

%% Set path to default
path(pathdef);

% Add [...] folder to path
addpath(genpath('functions\plots\'));
addpath(genpath('functions\'));
addpath(genpath('../Common\'));

%% Problem definition
muSun = astroConstants(4);
AU = astroConstants(2);

Neptune.ID = 8;    % Departure planet
Venus.ID = 2;      % Gravity assist planet
Mercury.ID = 1;    % Arrival planet

Neptune.mu = astroConstants(18);
Venus.mu = astroConstants(12);
Mercury.mu = astroConstants(11);

R_Venus = astroConstants(22);
h_p_min = 300;      % [km ]
%% Time windows

% Parabolic ToF calculation: planets orbits are assumed circular and
% coplanar. The radiuses of the approximated orbits are calculated at the
% earliest date possible.
clear var; close all; clc;

% Give dates
min_dep = [2031 02 01 00 00 00];
max_arr = [2071 02 01 00 00 00];
min_tdep = date2mjd2000(min_dep);
max_tarr = date2mjd2000(max_arr);

% Neptune -> Venus
[~, r_N, ~] = ephemeris(min_tdep, 8);  % Neptune
[~, r_V, ~] = ephemeris(min_tdep, 2);  % Venus
t_p_NV = parabolicTOF(norm(r_N), norm(r_V), muSun)*365; % NV ToF [ days ]

% Venus -> Mercury
[~, r_M, ~] = ephemeris(min_tdep, 1);  % Mercury
t_p_VM = parabolicTOF(norm(r_V), norm(r_M), muSun)*365; % VM ToF [ days ]

% Earliest dates of flyby and arrival
min_tfly = min_tdep + t_p_NV; 
min_fly = mjd20002date(min_tfly); % Earliest flyby
min_tarr = min_tfly + t_p_VM; 
min_arr = mjd20002date(min_tarr); % Earliest arrival

% Latest departure and flyby
max_tfly = max_tarr - t_p_VM; 
max_fly = mjd20002date(max_tfly); % Latest flyby
max_tdep = max_tfly - t_p_NV; 
max_dep = mjd20002date(max_tdep); % Latest departure

% Minimum radius of perigee at flyby, sum of radius of the planet and of the atmosphere
rpmin = R_Venus + h_p_min;   % [ km ]
lb = [min_tdep; min_tfly; min_tarr]; % lower boundary
ub = [max_tdep; max_tfly; max_tarr]; % upper boundary

%% Coarse grid search
departure = linspace(min_tdep, max_tdep, 300);
gravityassist = linspace(min_tfly, max_tfly, 300);
arrival = linspace(min_tarr, max_tarr, 300);

% preallocation
DV = zeros(length(departure),length(gravityassist),length(arrival));

% Parallel pool initialization
if max(size(gcp)) == 0 % parallel pool needed
    parpool % create the parallel pool
end

dep_size =  length(departure);
ga_size = length(gravityassist);
arr_size = length(arrival);

tStart = tic;
parfor i = 1 :dep_size
    for j = 1:ga_size
        tof_1 = gravityassist(j) - departure(i);
        for k = 1:arr_size
            tof_2 = arrival(k) - gravityassist(j);
            %% Set NaN for TOFs<parabolic TOF and rp < rpmin
            if  tof_2 < t_p_VM || tof_1 < t_p_NV 
                DV(i,j,k) = NaN;
                            
            else
                [DV(i,j,k), ~, ~, ~, FLYBY,~,~]= GAtransfer(Neptune,Venus,Mercury,departure(i),gravityassist(j),arrival(k));
                if FLYBY.rp < rpmin
                    DV(i,j,k) = 1000000;       % Set an impossibly high value
                end
             
            end
        end
    end
end

[dv_min_grid,loc] = min(DV(:));
[ii,jj,kk] = ind2sub(size(DV),loc);

x_GRID = [departure(ii); gravityassist(jj); arrival(kk)];

GRID_TIME = toc(tStart);

%% IMPORTANT SETS CONSTRAINTS ON x1<x2<x3 ( as a sys of 3 inequalities)
% x1<x2 , x2<x3, x1<x3

A = [1 -1 0;
    0 1 -1;
    1 0 -1];

b = [-t_p_NV; -t_p_VM; -(t_p_NV + t_p_VM)]; % DO NOT MULTIPLY BY 365

Aeq = [];
beq = [];

objfun = @(x) GAtransfer(Neptune,Venus,Mercury,x(1),x(2),x(3));

nonlincon = @(x) flyby_CONSTR(Neptune,Venus,Mercury,x(1),x(2),x(3),rpmin);

%% Constrained optimization function
options = optimoptions('fmincon','OptimalityTolerance',1e-12,'StepTolerance',1e-12);
x_PERF = fmincon(@(x) objfun(x),x_GRID,A,b,Aeq,beq,lb,ub,nonlincon,options);


%% Compute the positions and velocities of the planets at the optimum times
% as found by the GRID + FMINCON procedure
[DV, DV_dep, DV_arr, DV_ga, FLYBY, TRANSFER1, TRANSFER2]...
    = GAtransfer(Neptune, Venus, Mercury, x_PERF(1), x_PERF(2), x_PERF(3));

%% Genetic algorithm function - NOT USED
% if max(size(gcp)) == 0 % parallel pool needed
%     parpool % create the parallel pool
% end
% 
% numvars = 3;        % Three variables, tdep,tGA and tarr
% 
% GAopts = optimoptions('ga','UseParallel',true,'PopulationSize',500);
% % GAopts = optimoptions('ga','PopulationSize',3e6,'Display','off');
% [x_0, fval] = ga(objfun,numvars,A,b,Aeq,beq,lb,ub,nonlincon,GAopts);
% x_PERF_GA = fmincon(@(x) objfun(x),x_0,A,b,Aeq,beq,lb,ub,nonlincon,options);

%% Data for the assignment first and second arcs

% Transfer dates
dep = mjd20002date(x_PERF(1))
fly = mjd20002date(x_PERF(2))
arr = mjd20002date(x_PERF(3))

% first arc, initial point
[a1, e1, i1, RAAN1, omega1, f1_in] = car2kep(TRANSFER1.r1,TRANSFER1.v1,muSun);
a1 = a1/AU
e1 = e1
i1 = rad2deg(i1)
RAAN1 = rad2deg(RAAN1)
omega1 = rad2deg(omega1)
f1_in = rad2deg(f1_in)
r1_in = TRANSFER1.r1/AU
v1_in = TRANSFER1.v1
r1_in_mod = norm(TRANSFER1.r1)/AU
v1_in_mod = norm(TRANSFER1.v1)

% first arc, final point
[~, ~, ~, ~, ~, f1_fin] = car2kep(TRANSFER1.r2,TRANSFER1.v2,muSun);
f1_fin = rad2deg(f1_fin)
r1_fin = TRANSFER1.r2/AU
v1_fin = TRANSFER1.v2
r1_fin_mod = norm(TRANSFER1.r2)/AU
v1_fin_mod = norm(TRANSFER1.v2)


% second arc, initial point
[a2, e2, i2, RAAN2, omega2, f2_in] = car2kep(TRANSFER2.r1,TRANSFER2.v1,muSun);
a2 = a2/AU
e2 = e2
i2 = rad2deg(i2)
RAAN2 = rad2deg(RAAN2)
omega2 = rad2deg(omega2)
f2_in = rad2deg(f2_in)
r2_in = TRANSFER2.r1/AU
v2_in = TRANSFER2.v1
r2_in_mod = norm(TRANSFER2.r1)/AU
v2_in_mod = norm(TRANSFER2.v1)


% second arc, final point
[~, ~, ~, ~, ~, f2_fin] = car2kep(TRANSFER2.r2,TRANSFER2.v2,astroConstants(4));
f2_fin = rad2deg(f2_fin)
r2_fin = TRANSFER2.r2/AU
v2_fin = TRANSFER2.v2
r2_fin_mod = norm(TRANSFER2.r2)/AU
v2_fin_mod = norm(TRANSFER2.v2)

%% Data for the assignment flyby
clc;

% Venus' SOI radius on flyby date
mass_V = 4.867*10^24; % [ kg ]
mass_S = 1.989*10^30; % [ kg ]
[kep, ~] = uplanet(x_PERF(2), Venus.ID);
[r_V, ~] = kep2car(kep(1), kep(2), kep(3), kep(4), kep(5), kep(6), muSun);
r_V_SOI = norm(r_V)*(mass_V/mass_S)^(2/5)/AU

% Probe escape velocities in planetocentric RF
[~, ~, v_V] = ephemeris(x_PERF(2), Venus.ID);
v_minus = TRANSFER1.v2 - v_V
v_plus = TRANSFER2.v1 - v_V

% Hyperbolas parameters (modified flybypow before I have right results)
a = FLYBY.a
e = FLYBY.e
vp0 = FLYBY.vp0
ap = -Venus.mu/norm(v_plus)^2;
vp_plus = sqrt(Venus.mu*(1/abs(ap) + 2/rpmin))
theta_inf = rad2deg(acos(-1./FLYBY.e))
beta = 180 - theta_inf

%% Mission cost
dv_perc = DV_ga/DV*100

%% Transfer plots
clc; close all;

% Given orbits
figure(1) 
date = dep;
sun = 1;
obj = [1 2 8];
scale = [7000 2500  1000 30]; % DO NOT CHANGE!!!
orbit = 2;
Plot_Celestial_Objects(sun,obj,orbit,date,scale)
title('Given orbits')
legend('Mercury', 'Venus','Neptune')

figure(2)
obj = [1 2];
scale = [1000 500 10]; % DO NOT CHANGE!!!
Plot_Celestial_Objects(sun,obj,orbit,date,scale)
title('Given orbits')
legend('Mercury', 'Venus')

% transfer

figure(3) % transfer plus planets
r1 = [TRANSFER1.r1; TRANSFER2.r1];
r2 = [TRANSFER1.r2; TRANSFER2.r2];
v1 = [TRANSFER1.v1; TRANSFER2.v1];
v2 = [TRANSFER1.v2; TRANSFER2.v2];
disp = 1;
plotTrajectory(r1, r2, v1, v2, muSun, disp)
hold on;
date = dep;
sun = 1;
obj = [1 2 8];
scale = [7000 2500  1000 30]; % DO NOT CHANGE!!!
orbit = 2;
Plot_Celestial_Objects(sun,obj,orbit,date,scale)
title('Complete transfer trajectory')
legend('Orb. 1','Orb. 2','Arc 1','Arc 2','Mercury', 'Venus','Neptune')

figure(4) % transfer plus planets zoommed
disp = 1;
plotTrajectory(r1, r2, v1, v2, muSun, disp)
hold on;
obj = [1 2];
scale = [1000 500 10]; % DO NOT CHANGE!!!
Plot_Celestial_Objects(sun,obj,orbit,date,scale)
title('Complete transfer trajectory')
legend('Orb. 1','Orb. 2','Arc 1','Arc 2','Mercury', 'Venus')
axis equal

figure (5) % transfer only
plotTrajectory(r1, r2, v1, v2, muSun, disp)
title('Complete transfer trajectory')
legend('Orb. 1','Orb. 2','Arc 1','Arc 2')
axis equal

figure (6) % first arc only
plotTrajectory(TRANSFER1.r1, TRANSFER1.r2, TRANSFER1.v1, TRANSFER1.v2, muSun, disp)
title('First arc')
legend('Orb. 1','Arc 1')
axis equal

figure (7) % second arc only
plotTrajectory(TRANSFER2.r1, TRANSFER2.r2, TRANSFER2.v1, TRANSFER2.v2, muSun, disp)
title('Second arc')
legend('Orb. 2','Arc 2')
axis equal
grid on
%% plot hyperbola

% Plot the two planetocentric hyperbolic arcs
tspanA = linspace(0,-10000,10000);         % Plotting backward in time
[r0,v0A] = kep2car(FLYBY.a(1),FLYBY.e(1),0, 0,0,0,Venus.mu);
[r_A, ~] = propagator(r0,v0A,Venus.mu,tspanA);
tspanB = linspace(0,10000,10000);
[r0,v0B] = kep2car(a(2),e(2),0,0,0,0,Venus.mu);
[r_B, ~] = propagator(r0,v0B,Venus.mu,tspanB);

% Plot the arrival hyperbola
figure(8);
title('Venus flyby')
hold on;
legend on;
ARR = plot3(r_A(:,1),r_A(:,2),r_A(:,3), 'LineWidth', 2);
ARR.DisplayName = 'Incoming hyperbola';

% Plot the departure hyperbola
DEP = plot3(r_B(:,1),r_B(:,2),r_B(:,3), 'LineWidth', 2);
DEP.DisplayName = 'Outcoming hyperbola';

% Asymptotes not needed
hold on;
PlotObject(Venus.ID)