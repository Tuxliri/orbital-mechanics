clc
clearvars
close all

%% Set path to default
path(pathdef);
% Add [...] folder to path
addpath(genpath('functions\'));

addpath(genpath('../Common\'));

%% useful constants
muSun = astroConstants(4);
DAY2SECS = 24*60*60;

%% Input

% Design an interplanetary transfer with minimum DeltaV_tot between planetA
% and planetB, under the following mission requirements:

MISSION = {'Mercury',1,[2023 11 01; 2025 01 01], [2024 04 01; 2025 03 01], 7.0;
           'Venus',2,[2024 06 01; 2026 11 01], [2024 12 01; 2027 06 01], 3.0;
           'Mars',4,[2025 08 01; 2031 01 01], [2026 01 01; 2032 01 01],3.5;
           'Jupiter',5,[2026 06 01; 2028 06 01], [2028 06 01; 2034 01 01],9.1;
           'Mars Express',4,[2003 04 01; 2003 08 01], [2003 09 01; 2004 03 01],4;
           };
    

% Numbers of the mission to be selected
ID = 4;


%% Create a departure window vector
PROFILE = MISSION(ID,:);

planetA = 3;
planetB = PROFILE{2};
window_A = [PROFILE{3} zeros(2,3)];
window_B = [PROFILE{4} zeros(2,3)];
vinf = PROFILE{5};

%% calculate the total DVs
windowPlanetA = [date2mjd2000(window_A(1,1:6)) date2mjd2000(window_A(2,1:6))];

windowPlanetB = [date2mjd2000(window_B(1,1:6)) date2mjd2000(window_B(2,1:6))];


[DV, T1, T2, DV1, DV2,DEP_ORBIT,ARR_ORBIT,r0TF,v0TF] = tfdesigner(planetA,planetB,windowPlanetA,windowPlanetB,300,300);

%% find the minimum DV

min_DV = min(min(DV));

%% find the refined DV taking the approximated one as initial guess

%% plot the porkchop plot
% make a few levels
levels = linspace(floor(min_DV),floor(min_DV+10),11);  % fix the level for higher or lower DVs

% Convert mjd2000 vectors to date format
dates1 = zeros(1,length(T1));
dates2 = zeros(1,length(T2));

for i = 1:length(T1)
    dates1(i) = datenum(mjd20002date(T1(i)));
end

for i = 1:length(T2)
    dates2(i) = datenum(mjd20002date(T2(i)));
end

% plot the contour
contour(dates1,dates2,DV',levels,'ShowText','on');


% put axis labels to gregorian dates at 45 deg
xtickangle(45);
ytickangle(0);

datetick('x','yyyy mmm dd','keeplimits','keepticks')
datetick('y','yyyy mmm dd','keeplimits','keepticks')
ylabel('Arrival date');
xlabel('Departure date');



% put a nice colourbar

hbc = colorbar;
hbc.Title.String = '$\Delta v$ [km/s]';
hbc.Title.Interpreter = 'latex';
ha = gca;
ha.FontSize = 13;
hbc.Title.FontSize = 15;

grid on


%% Constant tof lines
hold on
caxis([14 20])
[C2,h2] = contour(dates1,dates2, dates2' - dates1,100:300:3000,'k');
clabel(C2,h2);

%% Find the index of the minimum of DVtot
 [i_min, k_min] = find(DV==min_DV);

% Find the departure date of the optimal Lambert transfer orbit
OPT_DEP_MJD2000 = T1(i_min);
OPT_ARR_MJD2000 = T2(k_min);
 
% Find the initial position and velocity of the optimal Lambert transfer orbit
OPT_r0 = reshape(r0TF(i_min,k_min,1:3),1,3);
OPT_v0 = reshape(v0TF(i_min,k_min,1:3),1,3);

%% calculate the transfer arc

% Calculate the period of the transfer orbit, to plot the remaining part
[a,~,~,~,~,~] = car2kep(OPT_r0,OPT_v0,muSun);
T_TF = 2*pi*sqrt(a^3 / muSun);              % [s]

tarc = linspace(OPT_DEP_MJD2000.*DAY2SECS,OPT_ARR_MJD2000.*DAY2SECS,1000);
[r_TF_ARC,v_TF_ARC] = propagator(OPT_r0,OPT_v0,muSun,tarc);
 
%% calculate the transfer orbit


%% calculate planet A orbit
[kepA,~,~] = ephemeris(OPT_DEP_MJD2000,planetA);
T_A = 2*pi*sqrt(kepA(1)^3 / muSun);              % orbital period of planetA [s]
tspanA = linspace(OPT_DEP_MJD2000,OPT_DEP_MJD2000 + T_A/DAY2SECS,1000);

[kep_el_A,R_A,V_A] = ephemeris(tspanA,planetA);

%% calculate planet B orbit
[kepB,~,~] = ephemeris(OPT_DEP_MJD2000,planetB);
T_B = 2*pi*sqrt(kepB(1)^3 / muSun);              % orbital period of planetB [s]

tspanB = linspace(OPT_DEP_MJD2000,OPT_DEP_MJD2000 + T_B/DAY2SECS,1000);

[kep_el_B,R_B,V_B] = ephemeris(tspanB,planetB);

%% Plotting
R = [r_TF_ARC R_A R_B];

figure(2);
hold on
plot3(r_TF_ARC(:,1),r_TF_ARC(:,2),r_TF_ARC(:,3), 'LineWidth', 2)
plot3(R_A(:,1),R_A(:,2),R_A(:,3), 'LineWidth', 2)
plot3(R_B(:,1),R_B(:,2),R_B(:,3), 'LineWidth', 2)

axis equal
grid on

view(30,45)
xlabel('x')
ylabel('y')
zlabel('z')