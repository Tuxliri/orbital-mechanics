% Test for the assignment 1 using the genetic algorithm
clc
clearvars
close all

%% Set path to default
path(pathdef);
% Add [...] folder to path
addpath(genpath('../functions_assignment\'));

addpath(genpath('../Common\'));

% Upper and lower bounds for the different mission phases
% x1 = tdep, x2 = tGA, t3 = tarr
%% Problem definition
Neptune.ID = 8;    % Departure planet
Venus.ID = 2;      % Gravity assist planet
Mercury.ID = 1;    % Arrival planet

Neptune.mu = astroConstants(18);
Venus.mu = astroConstants(12);
Mercury.mu = astroConstants(11);

%% Setting Constraints
min_tdep = date2mjd2000([2031 02 01 00 00 00]);
min_tarr = date2mjd2000([2050 02 01 00 00 00]);
min_tGA = date2mjd2000([2049 02 01 00 00 00]);

max_tdep = date2mjd2000([2033 02 01 00 00 00]);
max_tarr = date2mjd2000([2060 02 01 00 00 00]);
max_tGA = date2mjd2000([2059 02 01 00 00 00]);

% Launch vehicle and probe maximum deltav
LV_max_vinf = 120;       % Maximum hyperbolic excess velocity of launch vehicle [km/s]
SC_max_vinf = 100;       % Maximum deltav of the spacecraft [km/s]
rpmin = 6300;           % Minimum radius of perigee at flyby [km]

lb = [min_tdep; min_tGA; min_tarr];
ub = [max_tdep; max_tGA; max_tarr];

%%IMPORTANT SETS CONSTRAINTS ON x1<x2<x3 ( as a sys of 3 inequalities)
% x1<x2 , x2<x3, x1<x3

A = [1 -1 0;
    0 1 -1;
    1 0 -1];
b = [-300;-20;-350];
Aeq = [];
beq = [];

objfun = @(x) GAtransfer(Neptune,Venus,Mercury,x(1),x(2),x(3));

nonlincon = @(x) flyby_CONSTR(Neptune,Venus,Mercury,x(1),x(2),x(3), LV_max_vinf, SC_max_vinf,rpmin);

numvars = 3;        % Three variables, tdep,tGA and tarr

% Genetic algorithm function
% GAopts = optimoptions('ga','PopulationSize',1e15);
[x_0, fval] = ga(objfun,numvars,A,b,Aeq,beq,lb,ub,nonlincon);

% Constrained optimization function
options = optimoptions('fmincon','OptimalityTolerance',1e-12,'StepTolerance',1e-12);
x_PERF = fmincon(objfun,x_0,A,b,Aeq,beq,lb,ub,nonlincon);

%% Compute the positions and velocities of the planets at the optimum times
% [kepA, r_A, v_A] = ephemeris(x_OPT(1),Neptune);
% 
% [kepB, r_B, v_B] = ephemeris(x_OPT(2),Venus);   % Orbital radius of the GA planet
% 
% [kepC, r_C, v_C] = ephemeris(x_OPT(3),Mercury);
% 
% % plot them
% 
% NEP = plot3(r_A(1),r_A(2),r_A(3),'MarkerSize',30)
% grid on
% NEP = plot3(r_A(1),r_A(2),r_A(3),'MarkerSize',30)