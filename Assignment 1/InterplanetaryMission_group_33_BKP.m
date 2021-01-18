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

Neptune.Color = [62, 102, 249]./255;
Venus.Color = [254,235,186]./255;
Mercury.Color = [210, 210, 210]./255;


%% Setting Constraints
min_tdep = date2mjd2000([2031 02 01 00 00 00]);
min_tGA = date2mjd2000([2043 04 20 00 00 00]);
min_tarr = date2mjd2000([2043 04 28 00 00 00]);

max_tdep = date2mjd2000([2058 11 06 00 00 00]);
max_tGA = date2mjd2000([2071 01 23 00 00 00]);
max_tarr = date2mjd2000([2071 02 01 00 00 00]);

h_p_min = 300;      % [km ]
% Minimum radius of perigee at flyby, sum of radius of the planet and of the atmosphere
rpmin = R_Venus + h_p_min;   % [ km ]

lb = [min_tdep; min_tGA; min_tarr]; % lower boundary
ub = [max_tdep; max_tGA; max_tarr]; % upper boundary

% Calculate the parabolic ToF for departure conditions
[~,R_N0,~] = ephemeris(min_tdep,Neptune.ID);
[~,R_V0,~] = ephemeris(min_tdep,Venus.ID);
[~,R_M0,~] = ephemeris(min_tdep,Mercury.ID);

TOF_p_NV = parabolicTOF(norm(R_N0),norm(R_V0),muSun);
TOF_p_VM = parabolicTOF(norm(R_V0),norm(R_M0),muSun);

%% Coarse grid search
departure = linspace(min_tdep,max_tdep,300);
gravityassist = linspace(min_tGA,max_tGA,300);
arrival = linspace(min_tarr,max_tarr,300);

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
            if  tof_2 < TOF_p_VM*365 || tof_1 < TOF_p_NV*365
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

b = [-TOF_p_NV; -TOF_p_VM; -(TOF_p_NV+TOF_p_VM)]*365;

Aeq = [];
beq = [];

objfun = @(x) GAtransfer(Neptune,Venus,Mercury,x(1),x(2),x(3));

nonlincon = @(x) flyby_CONSTR(Neptune,Venus,Mercury,x(1),x(2),x(3),rpmin);

%% Constrained optimization function
options = optimoptions('fmincon','OptimalityTolerance',1e-12,'StepTolerance',1e-12);
x_PERF = fmincon(@(x) objfun(x),x_GRID,A,b,Aeq,beq,lb,ub,nonlincon);

%% Genetic algorithm function
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

%% Compute the positions and velocities of the planets at the optimum times
% as found by the GRID + FMINCON procedure
[DV, DV_dep, DV_arr, DV_ga, FLYBY, TRANSFER1, TRANSFER2]...
    = GAtransfer(Neptune, Venus, Mercury, x_PERF(1), x_PERF(2), x_PERF(3));

[outputArg1] = assignmentplot(Neptune,Venus,Mercury,x_PERF,TRANSFER1,TRANSFER2);
  
% assignmentplot(,x_PERF,TRANSFER1,TRANSFER2)
% Plot_Celestial_Objects(1,[Neptune.ID,Venus.ID,Mercury.ID],1,x_PERF(1))