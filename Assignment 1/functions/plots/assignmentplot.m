function [TF] = assignmentplot(planetA,planetB,planetC,datevec,TF1,TF2)

OPT_ARR_MJD2000 = datevec(3);
OPT_GA = datevec(2);
OPT_DEP_MJD2000 = datevec(1);

DAY2SECS = 24*3600;
%% calculate the transfer arc
muSun = astroConstants(4);
% Calculate the period of the first transfer orbit, to plot the remaining part
[a1,~,~,~,~,~] = car2kep(TF1.r1,TF1.v1,muSun);
T_TF1 = 2*pi*sqrt(a1^3 / muSun);              % [s]

% Calculate the period of the first transfer orbit, to plot the remaining part
[a2,~,~,~,~,~] = car2kep(TF2.r1,TF2.v2,muSun);
T_TF2 = 2*pi*sqrt(a2^3 / muSun);              % [s]

tarc1 = linspace(OPT_DEP_MJD2000,OPT_GA,10000);
tarc2 = linspace(OPT_GA,OPT_ARR_MJD2000,1000);

[r_TF_ARC1,~] = propagator(TF1.r1,TF1.v1,muSun,tarc1.*DAY2SECS);
[r_TF_ARC2,~] = propagator(TF2.r1,TF2.v1,muSun,tarc2.*DAY2SECS);
%% calculate the transfer orbit1
torbit1 = linspace(OPT_ARR_MJD2000.*DAY2SECS,OPT_GA.*DAY2SECS + T_TF1,1000);
[r_TF_ORBIT1,~] = propagator(TF1.r1,TF1.v1,muSun,torbit1);

%% calculate the transfer orbit2
torbit2 = linspace(OPT_GA.*DAY2SECS,OPT_DEP_MJD2000.*DAY2SECS + T_TF2,1000);
[r_TF_ORBIT2,~] = propagator(TF2.r1,TF2.v1,muSun,torbit2);

%% calculate planet A orbit

[kep_el_A,R_A,V_A,R_A_dep,R_A_GA,R_A_arr] = planet_orbit(planetA.ID,OPT_DEP_MJD2000,OPT_GA,OPT_ARR_MJD2000,muSun);

%% calculate planet B orbit

[kep_el_B,R_B,V_B,R_B_dep,R_B_GA,R_B_arr] = planet_orbit(planetB.ID,OPT_DEP_MJD2000,OPT_GA,OPT_ARR_MJD2000,muSun);

%% calculate planet C orbit

[kep_el_C,R_C,V_C,R_C_dep,R_C_GA,R_C_arr] = planet_orbit(planetC.ID,OPT_DEP_MJD2000,OPT_GA,OPT_ARR_MJD2000,muSun);

%% Plotting

figure(3);
hold on
legend on
TF = plot3(r_TF_ARC1(:,1),r_TF_ARC1(:,2),r_TF_ARC1(:,3), 'LineWidth', 2);
TF.Color = 'r';
TF.DisplayName = 'Transfer Arc 1';

hold on
legend on
TF = plot3(r_TF_ARC2(:,1),r_TF_ARC2(:,2),r_TF_ARC2(:,3), 'LineWidth', 2);
TF.Color = 'b';
TF.DisplayName = 'Transfer Arc 2';

PA_ORBIT = plot3(R_A(:,1),R_A(:,2),R_A(:,3), 'LineWidth', 2) ;  
PA_ORBIT.DisplayName = 'departure planet orbit';
PA_ORBIT.Color = planetA.Color;

PB_ORBIT = plot3(R_B(:,1),R_B(:,2),R_B(:,3), 'LineWidth', 2);   
PB_ORBIT.DisplayName = 'arrival planet orbit';
PB_ORBIT.Color = planetB.Color;

PC_ORBIT = plot3(R_C(:,1),R_C(:,2),R_C(:,3), 'LineWidth', 2);   
PC_ORBIT.DisplayName = 'arrival planet orbit';
PC_ORBIT.Color = planetC.Color;

% plot initial position of departure planet

DEP_A= plot3(R_A_dep(1),R_A_dep(2),R_A_dep(3));
DEP_A.Marker = 'o';
DEP_A.MarkerFaceColor = planetA.Color;
DEP_A.MarkerSize = 10;
DEP_A.HandleVisibility = 'off';

% plot FLYBY position of planet B

ARR_B = plot3(R_B_GA(1),R_B_GA(2),R_B_GA(3));
ARR_B.Marker = 'o';
ARR_B.MarkerFaceColor = planetB.Color;
ARR_B.MarkerSize = 10;
ARR_B.HandleVisibility = 'off';

% plot final position of arrival planet

ARR_C = plot3(R_C_arr(1),R_C_arr(2),R_C_arr(3));
ARR_C.Marker = 'o';
ARR_C.MarkerFaceColor = planetC.Color;
ARR_C.MarkerSize = 10;
ARR_C.HandleVisibility = 'off';

% plot transfer windows

axis equal
grid on



view(30,30)
xlabel('x')
ylabel('y')
zlabel('z')
end
