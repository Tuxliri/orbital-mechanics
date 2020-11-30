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
DAY2SECS = 24*60*60;

%% Input

% Design an interplanetary transfer with minimum DeltaV_tot between planetA
% and planetB, under the following mission requirements:

MISSION = {'Mercury',1,[2023 11 01; 2025 01 01], [2024 04 01; 2025 03 01], 7.0,[210,210,210];
           'Venus',2,[2024 06 01; 2026 11 01], [2024 12 01; 2027 06 01], 3.0,[255,235,186];
           'Mars',4,[2025 08 01; 2031 01 01], [2026 01 01; 2032 01 01],3.5,[217, 83, 25];
           'Jupiter',5,[2026 06 01; 2028 06 01], [2028 06 01; 2034 01 01],9.1,[201,144,57];
           'Saturn',6,[2027 09 01; 2029 10 01], [2030 04 01; 2036 03 01],11.5,[200,130,60];
           'Uranus',7,[2027 01 01; 2029 01 01], [2031 04 01; 2045 12 01],12.1,[213, 251, 252];
           'Neptune',8,[2025 01 01; 2026 10 01], [2036 01 01; 2055 06 01],12.5,[62, 102, 249];
           'Mars Express',4,[2003 04 01; 2003 08 01], [2003 09 01; 2004 03 01],4,[217, 83, 25];
           };
    

% Numbers of the mission to be selected
ID = 5;


%% Create a departure window vector
PROFILE = MISSION(ID,:);

planetA = 3;
planetB = PROFILE{2};
window_A = [PROFILE{3} zeros(2,3)];
window_B = [PROFILE{4} zeros(2,3)];
vinf = PROFILE{5};
A_color = [0, 114, 189]./254;
if length(PROFILE) > 5
    B_color = PROFILE{6}./254;
end

%% calculate the total DVs

windowPlanetA = [date2mjd2000(window_A(1,1:6)) date2mjd2000(window_A(2,1:6))];

windowPlanetB = [date2mjd2000(window_B(1,1:6)) date2mjd2000(window_B(2,1:6))];


[DV, T1, T2, DV1, DV2,DEP_ORBIT,ARR_ORBIT,r0TF,v0TF] = tfdesigner(planetA,planetB,windowPlanetA,windowPlanetB,100,100);

%% find the initial guess minimum DV
% constrained problem: the required DV1 must be less than the maximum
% excess velocity from the launcher vinf
min_DV0 = min(DV(DV1 < vinf));

% index of the minimum of DVtot
[i_min, k_min] = find(DV==min_DV0);

% Find the departure date of the optimal Lambert transfer orbit
OPT_DEP_MJD2000 = T1(i_min);
OPT_ARR_MJD2000 = T2(k_min);

%% find the minmum DV dates in MJD2000 format through numerical optimization
options = optimoptions('fminunc','OptimalityTolerance',1e-12,'StepTolerance',1e-12);

% Departure and arrival dates vector
% x = fminunc(@(x) designer_opt(planetA,planetB,x(1),x(2)), [OPT_DEP_MJD2000 OPT_ARR_MJD2000],options);
% 
% OPT_DEP_MJD2000 = x(1);
% OPT_ARR_MJD2000 = x(2);

% Find the orbital paremeters of the optimal transfer orbit

[min_DV, ~, ~, ~, ~,~, ~,OPT_r0,OPT_v0]...
    = tfdesigner(planetA,planetB,OPT_DEP_MJD2000,OPT_ARR_MJD2000,1,1);

OPT_r0 = reshape(OPT_r0,1,3);
OPT_v0 = reshape(OPT_v0,1,3);

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
caxis([levels(1) levels(end)])
[C2,h2] = contour(dates1,dates2, dates2' - dates1,100:300:3000,'k');
clabel(C2,h2);

%% calculate the transfer arc

% Calculate the period of the transfer orbit, to plot the remaining part
[a,~,~,~,~,~] = car2kep(OPT_r0,OPT_v0,muSun);
T_TF = 2*pi*sqrt(a^3 / muSun);              % [s]

tarc = linspace(OPT_DEP_MJD2000.*DAY2SECS,OPT_ARR_MJD2000.*DAY2SECS,1000);
[r_TF_ARC,v_TF_ARC] = propagator(OPT_r0,OPT_v0,muSun,tarc);
 
%% calculate the transfer orbit
torbit = linspace(OPT_ARR_MJD2000.*DAY2SECS,OPT_DEP_MJD2000.*DAY2SECS + T_TF,1000);
[r_TF_ORBIT,v_TF_ORBIT] = propagator(r_TF_ARC(end,:),v_TF_ARC(end,:),muSun,torbit);

%% calculate planet A orbit

[kep_el_A,R_A,V_A,R_A_dep,R_A_arr] = planet_orbit(planetA,OPT_DEP_MJD2000,OPT_ARR_MJD2000,muSun);

%% calculate planet B orbit

[kep_el_B,R_B,V_B,R_B_dep,R_B_arr] = planet_orbit(planetB,OPT_DEP_MJD2000,OPT_ARR_MJD2000,muSun);

%% Plotting

figure(3);
hold on
legend on
TF = plot3(r_TF_ARC(:,1),r_TF_ARC(:,2),r_TF_ARC(:,3), 'LineWidth', 2);
TF.Color = '#77AC30';
TF.DisplayName = 'Transfer Arc';

TF_ORBIT = plot3(r_TF_ORBIT(:,1),r_TF_ORBIT(:,2),r_TF_ORBIT(:,3),'--', 'LineWidth', 2);
TF_ORBIT.DisplayName = 'Transfer orbit';
TF_ORBIT.Color = TF.Color;

PA_ORBIT = plot3(R_A(:,1),R_A(:,2),R_A(:,3), 'LineWidth', 2) ;  
PA_ORBIT.DisplayName = 'departure planet orbit';
PA_ORBIT.Color = A_color;

PB_ORBIT = plot3(R_B(:,1),R_B(:,2),R_B(:,3), 'LineWidth', 2);   
PB_ORBIT.DisplayName = 'arrival planet orbit';
PB_ORBIT.Color = B_color;

% plot initial and final position of departure planet

ARR_A= plot3(R_A_arr(1),R_A_arr(2),R_A_arr(3));
ARR_A.Marker = 'o';
ARR_A.MarkerFaceColor = PA_ORBIT.Color;
ARR_A.MarkerSize = 10;
ARR_A.HandleVisibility = 'off';

DEP_A= plot3(R_A_dep(1),R_A_dep(2),R_A_dep(3));
DEP_A.Marker = 'o';
DEP_A.MarkerFaceColor = PA_ORBIT.Color;
DEP_A.MarkerSize = 10;
DEP_A.HandleVisibility = 'off';

% plot initial and final position of arrival planet

ARR_B = plot3(R_B_arr(1),R_B_arr(2),R_B_arr(3),'-o','MarkerFaceColor','#BC2731','MarkerSize',10);
ARR_B.HandleVisibility = 'off';

DEP_B = plot3(R_B_dep(1),R_B_dep(2),R_B_dep(3),'-o','MarkerFaceColor','#BC2731','MarkerSize',10);
DEP_B.HandleVisibility = 'off';

% plot transfer windows
DEP_WIN=plot3(DEP_ORBIT(:,1),DEP_ORBIT(:,2),DEP_ORBIT(:,3),'Color',[PA_ORBIT.Color 0.5],'LineWidth', 10);
DEP_WIN.DisplayName = 'Departure Window';

ARR_WIN=plot3(ARR_ORBIT(:,1),ARR_ORBIT(:,2),ARR_ORBIT(:,3),'Color',[PB_ORBIT.Color 0.5],'LineWidth', 10);
ARR_WIN.DisplayName = 'Arrival Window';


axis equal
grid on



view(30,30)
xlabel('x')
ylabel('y')
zlabel('z')