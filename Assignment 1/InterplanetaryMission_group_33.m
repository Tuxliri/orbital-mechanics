% Test for the assignment 1 using the genetic algorithm
clc
clearvars
close all

%% Set path to default
path(pathdef);
% Add [...] folder to path
addpath(genpath('functions\'));

addpath(genpath('../Common\'));

% Upper and lower bounds for the different mission phases
% x1 = tdep, x2 = tGA, t3 = tarr
%% Problem definition
Neptune.ID = 8;    % Departure planet
Venus.ID = 3;      % Gravity assist planet
Mercury.ID = 1;    % Arrival planet

Neptune.mu = astroConstants(18);
Venus.mu = astroConstants(12);
Mercury.mu = astroConstants(11);

Neptune.Color = [62, 102, 249]./255;
Venus.Color = [254,235,186]./255;
Mercury.Color = [210,210,210]./255;

%% Setting Constraints
min_tdep = date2mjd2000([2031 02 01 00 00 00]);
min_tarr = date2mjd2000([2043 03 01 00 00 00]);
min_tGA = date2mjd2000([2043 02 01 00 00 00]);

max_tdep = date2mjd2000([2058 02 01 00 00 00]);
max_tGA = date2mjd2000([2070 02 01 00 00 00]);
max_tarr = date2mjd2000([2071 02 01 00 00 00]);

% Minimum radius of perigee at flyby, sum of radius of the planet and of the atmosphere [km]
rpmin = 6300;   

lb = [min_tdep; min_tGA; min_tarr];
ub = [max_tdep; max_tGA; max_tarr];

%% Coarse grid search
departure = linspace(min_tdep,max_tdep,100);
gravityassist = linspace(min_tGA,max_tGA,100);
arrival = linspace(min_tarr,max_tarr,100);

% preallocation
DV = zeros(length(departure),length(gravityassist),length(arrival));

for i = 1 : length(departure)
    for j = 1:length(gravityassist)
        tof_1 = gravityassist(j) - departure(i);
        for k = 1:length(arrival)
            tof_2 = arrival(k) - gravityassist(j);
            %% Set NaN for negative TOFs and rp < rpmin
            if  tof_2 < 0 || tof_1 < 0
                DV(i,j,k) = NaN;
                
            
            else
            [DV(i,j,k), ~, ~, FLYBY,~,~]...
                 = GAtransfer(Neptune,Venus,Mercury,departure(i),gravityassist(j),arrival(k));
             
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

%% IMPORTANT SETS CONSTRAINTS ON x1<x2<x3 ( as a sys of 3 inequalities)
% x1<x2 , x2<x3, x1<x3

A = [1 -1 0;
    0 1 -1;
    1 0 -1];

b = [-300; -20; -350];

Aeq = [];
beq = [];

objfun = @(x) GAtransfer(Neptune,Venus,Mercury,x(1),x(2),x(3));

nonlincon = @(x) flyby_CONSTR(Neptune,Venus,Mercury,x(1),x(2),x(3),rpmin);

%% Constrained optimization function
options = optimoptions('fmincon','OptimalityTolerance',1e-12,'StepTolerance',1e-12,'MaxFunctionEvaluations',3000, ...
'MaxIterations',1000    );
x_PERF = fmincon(@(x) objfun(x),x_GRID,A,b,Aeq,beq,lb,ub,nonlincon);

%% Genetic algorithm function
numvars = 3;        % Three variables, tdep,tGA and tarr

GAopts = optimoptions('ga','PopulationSize',3e20,'MaxStallTime',120);
[x_0, fval] = ga(objfun,numvars,A,b,Aeq,beq,lb,ub,nonlincon);

for i = 1:3
    mjd20002date(x_PERF(i))
end

%% Compute the positions and velocities of the planets at the optimum times
% as found by the GRID + FMINCON procedure

[DV, DV_dep, DV_arr, FLYBY,TRANSFER1,TRANSFER2]...
    = GAtransfer(Neptune,Venus,Mercury,x_PERF(1),x_PERF(2),x_PERF(3));


    
assignmentplot(Neptune,Venus,Mercury,x_PERF,TRANSFER1,TRANSFER2)

