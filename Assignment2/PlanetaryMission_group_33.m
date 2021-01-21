% Assignment 2 : Planetary Explorer
clc
clear all
close all

%% Set path to default
path(pathdef);
% Add [...] folder to path
addpath(genpath('functions/'));
addpath(genpath('../Common/'));

% Conversion constants
DAY2SECS=24*3600;

% DATA

const = astroConstants([23 13 9]);
R_E = const(1);					% radius of the Earth                   [ km ]
muE = const(2);					% gravitational parameter of Earth      [ km^3/s^2 ]
J2 = const(3);					% second zonal armonic of earth         [ - ]
omega_e = deg2rad(15.04)/3600;  % angular velocity of Earth's rotation  [ rad/s ]
gw_longitude0 = 0;              % longitude of greenwhich   at time t0  [ rad ]
Psr = 4.56e-6;                  % Solar radiation pressure at 1AU       [N/m^2]

NDAYS = 1200;                      % Number of days to propagate for perturbations [days]

%% Assigned orbit parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOOSE THE DESIRED SATELLITE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1 - ASSIGNED SATELLITE
% 2 - Validation of SRP perturbation
% 3 - "POLAR" SATELLITE Real TLEs comparison
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
selection = 1;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

switch selection
    case 1  % ASSIGNED SATELLITE
        a = 40718;		% [km]	semi-major axis
        e = 0.6177; 	% [-]	eccentricity
        i = deg2rad(78.2195);	% [deg]	inclination
        RAAN = deg2rad(0);
        omega = deg2rad(40);
        f0 = deg2rad(0);
        hp = 15566.491;	% [km]	height of perigee
        t0 = 0;
        k = 1;
        m = 1;
        periods = 10;               % number of periods to plot
        
        date0 = [2021 01 01 00 00 00];  % Initial date at time t=0
        
        % SRP perturbation data
        Cr = 1.2;		% [-] reflectivity coefficient
        Am = 4;			% [m^2/kg] area-to-mass ratio
        
    case 2
        % Validation of SRP perturbation, data from Curtis Example 10.9
        % here we validate our model using only SRP as a perturbation
        a = 10085.44;
        e = 0.025422;
        i = deg2rad(88.3924);
        RAAN = deg2rad(45.3812);
        omega = deg2rad(227.493);
        f0 = deg2rad(343.427);
        date0 = jd2date(2438400.5);  % Initial date at time t=0
        k = 1;
        m = 1;
        periods = 4;               % number of periods to plot
        
        J2 = 0;                      % Overwrite J2 coefficient to test SRP
        t0 = 0;
        
        % model
        
        % SRP perturbation data
        Cr = 2;		% [-] reflectivity coefficient
        Am = 2;			% [m^2/kg] area-to-mass ratio
        
    case 3 % Real TLEs comparison
        
        a = 3.487516398368765E+04;		% [km]	semi-major axis
        e = 6.931588472978080E-01; 	% [-]	eccentricity
        i = deg2rad(84.72791562472406);	% [deg]	inclination
        RAAN = deg2rad(1.123750117977185E+01);
        omega = deg2rad(2.244379313913353E+02);
        f0 = deg2rad(1.226135757537926E+02);
        
		%ID:23802 POLAR
        t0 = 0;
        k = 1;
        m = 1;
        periods = 10;               % number of periods to plot
        
        date0 = [2000 01 02 00 00 00];  % Initial date at time t=0
        
        % SRP perturbation data
        mass = 1297;    % [kg]
        area = 2.4*1.8; % [m^2]
        Cr = 1.2;		% [-] reflectivity coefficient
        Am = area/mass;			% [m^2/kg] area-to-mass ratio
        NDAYS = 1035; 
        
end

%% Calculate the semi-major axis of the repeating ground track

a_repeating = repeating_ground_track(m,k,muE,omega_e);

%% Calculate the semi-major axis of the perturbed repeating ground track

a_secular = repeating_ground_track(m, k, muE, omega_e, J2, R_E, e, i);

%% Create the state vectors from orbital parameters
% calculate the orbital periods
T = 2*pi*sqrt(a^3/muE);                         % Original
T_repeating = 2*pi*sqrt(a_repeating^3/muE);     % Repeating ground track
T_secular = 2*pi*sqrt(a_secular^3/muE);         % Repeating secular ground track for the perturbed secular orbit

%% create state vectors

state_vec = [a e i RAAN omega f0];                   % Original
state_vec_repeating = [a_repeating e i RAAN omega f0];           % Repeating ground track
state_vec_secular = [a_secular e i RAAN omega f0];   % Repeating ground track for the perturbed secular orbit


%% Create time vectors
for j = 1 : 3
    if j == 1
        % one period
        tEnd1=T;
        tEnd2=T_repeating;
        tEnd3=T_secular;
    elseif j == 2
        % one day
        tEnd1=24*3600;
        tEnd2=tEnd1;
        tEnd3=tEnd1;
    else
        % ten days
        tEnd1=10*24*3600;
        tEnd2=tEnd1;
        tEnd3=tEnd1;
    end
    
    t = (0:1:tEnd1);                    	% Original
    t_repeating = (0:1:tEnd2);	% Repeating ground track
    t_secular = (0:1:tEnd3);    	% Repeating secular ground track for the perturbed secular orbit
    
    
    %% Compute RA, declination, lon and latitude
    % Unperturbed original
    [alpha, delta, lon, lat] = groundtrack(state_vec, gw_longitude0, t, omega_e, muE, t0);
    
    % J2 perturbed original
    [alpha_j2, delta_j2, lon_j2, lat_j2] = groundtrack(state_vec, gw_longitude0, t, omega_e, muE, t0, J2, R_E);
    
    % Repeating unperturbed ground track
    [alpha_repeating, delta_repeating, lon_repeating, lat_repeating] = groundtrack(state_vec_repeating, gw_longitude0, t_repeating, omega_e, muE, t0);
    
    % Repeating perturbed ground track
    [alpha_secular, delta_secular, lon_secular, lat_secular] = groundtrack(state_vec_secular,gw_longitude0,t_secular,omega_e, muE,t0, J2, R_E);
    
    %% Plotting of the groundtracks
    % unperturbed J2 perturbed original for 1 orbit, 1 day and 10 days
    figure
    
    % Unperturbed original
    plot_groundtrack(lon,lat,'y');
    set(gca,'FontSize',20)
    % J2 perturbed original
    plot_groundtrack(lon_j2,lat_j2,'r');
    set(gca,'FontSize',20)
    legend ('Unperturbed original', 'Start', 'End','J2 perturbed original', 'Start', 'End','Orientation','horizontal','Location','northoutside' );
    
    switch j
        case 1
            title('GT of the unpert. 2BP and of the 2BP pert. by J2 - 1 orbit')
            
            % unperturbed and perturbed repeating GT
            figure()
            
            % Repeating ground track
            plot_groundtrack(lon_repeating,lat_repeating, 'y' );
            set(gca,'FontSize',20)
            % J2 Perturbed ground track
            plot_groundtrack(lon_secular, lat_secular, 'r');
            set(gca,'FontSize',20)
            legend ('Unperturbed repeating GT', 'Start', 'End','J2 perturbed repeating GT', 'Start', 'End','Orientation','horizontal','Location','northoutside','FontSize',20);
            title('repeating GT of the unpert. 2BP and of the 2BP pert. by J2 - 1 orbit','FontSize',20)    
        case 2
            title('GT of the unpert. 2BP and of the 2BP pert. by J2 - 1 day','FontSize',20)
        case 3
            title('GT of the unpert. 2BP and of the 2BP pert. by J2 - 10 days','FontSize',20)
    end
end

%%

% Initial conditions
kep0 = [a e i RAAN omega f0];
tspan = t0:100:NDAYS*DAY2SECS;

%% Get the perturbed orbital elements
[t_gauss,kep_gauss] = ORBITPROPAGATOR(t0,kep0,tspan,date0,J2,Cr,Psr,Am);

tStart = tic;

% Initial cartesian elements
[r0, v0] =  kep2car(a,e,i,RAAN,omega,f0,muE);

% Find out how much time the gaussian propagation took
tGAUSS = toc(tStart);

%% Solve with the cartesian elements
kepB = zeros(length(tspan),6);

tStart = tic;

[ri, vi] = propagator(r0,v0,muE,tspan,J2,R_E,date0,Cr,Psr,Am);

% Convert cartesian to keplerian
for j = 1:length(ri(:,1))
    
    [kepB(j,1),kepB(j,2),kepB(j,3),kepB(j,4),kepB(j,5),kepB(j,6)] = car2kep(ri(j,:),vi(j,:),muE);
end

% Find out how much time the cartesian propagation took
tCART = toc(tStart);

%% Filtering lower frequencies

% Cut-off period
Tfilter = 3*T;

% Number of points for the filtering window
nwindow = nearest( Tfilter / (sum(diff(t_gauss)) / (numel(t_gauss)-1) ) );

% Filter elements ( no unwrapping)
kep_filtered = movmean( kep_gauss,nwindow,1);

%% Compare the two solutions through plotting
% wrapping


% Add a little gap to cartesian elements to avoid numerical problems

kepB(:,3:5) = unwrap(kepB(:,3:5),[],1);
kepB(5:end,6) = unwrap(kepB(5:end,6),[],1);

kep_gauss(:,3:6) = unwrap(kep_gauss(:,3:6),[],1);

kep_filtered(:,3:6) = unwrap(kep_filtered(:,3:6),[],1);

kepB(:,3:6) = rad2deg(kepB(:,3:6));
kep_gauss(:,3:6) = rad2deg(kep_gauss(:,3:6));
kep_filtered(:,3:6) = rad2deg(kep_filtered(:,3:6));

tspan = tspan./DAY2SECS;
LINEWIDTH = 2;

%% Plotting
% Semi-major axis
figure(2)

subplot(2,3,1)
plot(tspan,kep_gauss(:,1)-a,tspan,kepB(:,1)-a,tspan,kep_filtered(:,1)-a);
legend('Gauss equations','Cartesian','Secular (filtered)')
grid on
xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('$\mathbf{a-a_0 [Km]}$','Interpreter', 'latex','Fontsize', 14)

% Eccentricity
subplot(2,3,2)
plot(tspan,kep_gauss(:,2)-e,tspan,kepB(:,2)-e,tspan,kep_filtered(:,2)-e);
legend('Gauss equations','Cartesian','Secular (filtered)')
grid on
xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('$\mathbf{e-e_0 [-]}$','Interpreter', 'latex','Fontsize', 14)


% inclination
subplot(2,3,3)
plot(tspan,kep_gauss(:,3)-rad2deg(i),tspan,kepB(:,3)-rad2deg(i),tspan,kep_filtered(:,3)-rad2deg(i));
legend('Gauss equations','Cartesian','Secular (filtered)')
grid on
xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('$\mathbf{i-i_0 [deg]}$','Interpreter', 'latex','Fontsize', 14)

% RAAN
subplot(2,3,4)
plot(tspan,kep_gauss(:,4)-rad2deg(RAAN),tspan,kepB(:,4)-rad2deg(RAAN),tspan,kep_filtered(:,4)-rad2deg(RAAN));
legend('Gauss equations','Cartesian','Secular (filtered)')
grid on
xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('$\mathbf{\Omega -\Omega_0 [deg]}$','Interpreter', 'latex','Fontsize', 14)

% omega
subplot(2,3,5)
plot(tspan,kep_gauss(:,5)-rad2deg(omega),tspan,kepB(:,5)-rad2deg(omega),tspan,kep_filtered(:,5)-rad2deg(omega));
legend('Gauss equations','Cartesian','Secular (filtered)')
grid on
xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('$\mathbf{\omega -\omega_0  [deg]}$','Interpreter', 'latex','Fontsize', 14)


% f
subplot(2,3,6)
plot(tspan,kep_gauss(:,6)-deg2rad(f0),tspan,kepB(:,6)-deg2rad(f0),tspan,kep_filtered(:,6)-deg2rad(f0));
legend('Gauss equations','Cartesian','Secular (filtered)')
grid on
xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('$\mathbf{f-f_0  [deg]}$','Interpreter', 'latex','Fontsize', 14)

%% Compare the difference between the two solutions through plotting
figure(3)

% Semi-major axis
subplot(2,3,1)
semilogy(tspan,abs((kep_gauss(:,1) - kepB(:,1)))./kep0(1));
grid on
xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('${|a_{Car} - a_{Gauss}| / a0 [-]}$','Interpreter', 'latex')

% Eccentricity
subplot(2,3,2)
semilogy(tspan,abs(kep_gauss(:,2) - kepB(:,2)));
grid on
xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('${|e_{Car} - e_{Gauss}| [-]}$','Interpreter', 'latex')


% inclination
subplot(2,3,3)
semilogy(tspan,abs(kep_gauss(:,3) - kepB(:,3)) ./(2*pi));
grid on
xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('${|i_{Car} - i_{Gauss}| / 2\pi [-]}$','Interpreter', 'latex')

% RAAN
subplot(2,3,4)
semilogy(tspan,abs(kep_gauss(:,4) - kepB(:,4)) ./(2*pi));
grid on
xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('${|\Omega_{Car} - \Omega_{Gauss}| / 2\pi [-]}$','Interpreter', 'latex')

% omega
subplot(2,3,5)
semilogy(tspan,abs(kep_gauss(:,5) - kepB(:,5)) ./(2*pi));
grid on
xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('${|\omega_{Car} - \omega_{Gauss}| / 2\pi [-]}$','Interpreter', 'latex')


% omega
subplot(2,3,6)
semilogy(tspan,abs(kep_gauss(:,6) - kepB(:,6)) ./ kep0(5));
grid on
xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('${|f_{Car} - f_{Gauss}| / f_{Gauss} [-]}$','Interpreter', 'latex')

%% Plot REAL ORBITAL Elements
if selection == 3
    load('Polar_Real_Time_Elements_hourly.mat')

    % Creating a table containing the real time orbital elements
    TABLE = table(REALTLES);
    %Converting table to array
    POLAR = table2array(TABLE);

    % Creating a vector of julian days for the the TLEs
    DATES = POLAR(:,1);
    INITIALDATE = jd2date(DATES(1));
    FINALDATE = jd2date(DATES(end));

    % Assigning the six orbital elements to KEP matrix

    % Extracting semi-major axis vector
    KEP(:,1) = POLAR(:,11);
    % Extracting eccentricity vector
    KEP(:,2) = POLAR(:,2);
    % Extracting inclination vector
    KEP(:,3) = POLAR(:,4);
    % Extracting RAAN vector
    KEP(:,4) = POLAR(:,5);
    % Extracting argument of perigee vector
    KEP(:,5) = POLAR(:,6);
    % Extracting true anomaly vector
    KEP(:,6) = POLAR(:,10);

    KEP(:,6) = deg2rad(KEP(:,6));
    KEP(3:end,6) = unwrap(KEP(3:end,6),[],1);
    KEP(:,6) = rad2deg(KEP(:,6));

    DATES = DATES-DATES(1);

    figure(4)

    % Semi-major axis
    subplot(3,2,1)
    plot(DATES,KEP(:,1),'b.',tspan,kepB(:,1),'r-',tspan,kep_gauss(:,1),'k-','LineWidth',2);
    legend('TLEs','Cartesian','Gauss equations');
    grid on
    xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
    ylabel('$\mathbf{a [Km]}$','Interpreter', 'latex','Fontsize', 14)
    % saveas

    % Eccentricity
    % jd = juliandate(DATES);
    subplot(3,2,2)
    plot(DATES,KEP(:,2),'b.',tspan,kepB(:,2),'r-',tspan,kep_gauss(:,2),'k-','LineWidth',2);
    legend('TLEs','Cartesian','Gauss equations');
    grid on
    xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
    ylabel('$\mathbf{e [-]}$','Interpreter', 'latex','Fontsize', 14)


    % inclination
    subplot(3,2,3)
    plot(DATES,KEP(:,3),'b.',tspan,kepB(:,3),'r-',tspan,kep_gauss(:,3),'k-','LineWidth',2);
    legend('TLEs','Cartesian','Gauss equations');
    grid on
    xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
    ylabel('$\mathbf{i [deg]}$','Interpreter', 'latex','Fontsize', 14)

    % RAAN
    subplot(3,2,4)
    plot(DATES,KEP(:,4),'b.',tspan,kepB(:,4),'r-',tspan,kep_gauss(:,4),'k-','LineWidth',2);
    legend('TLEs','Cartesian','Gauss equations');
    grid on
    xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
    ylabel('$\mathbf{\Omega  [deg]}$','Interpreter', 'latex','Fontsize', 14)

    % omega
    subplot(3,2,5)
    plot(DATES ,KEP(:,5),'b.',tspan,kepB(:,5),'r-',tspan,kep_gauss(:,5),'k-','LineWidth',2);
    legend('TLEs','Cartesian','Gauss equations');
    grid on
    xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
    ylabel('$\mathbf{\omega  [deg]}$','Interpreter', 'latex','Fontsize', 14)


    % f
    subplot(3,2,6)
    plot(DATES,KEP(:,6),'b.',tspan,kepB(:,6),'r-',tspan,kep_gauss(:,6),'k-','LineWidth',2);
    legend('TLEs','Cartesian','Gauss equations');
    grid on
    xlabel('${time [days]}$','Interpreter', 'latex','Fontsize', 14)
    ylabel('$\mathbf{f  [deg]}$','Interpreter', 'latex','Fontsize', 14)
end
%% Plot the Orbit evolution
for j=1:2
    %number of orbits
    switch j
        case 1
            norb=1;
        case 2
            norb=1000;
    end
    pointsperorb=1000;
    figure
    axis equal
    rin=r0;
    vin=v0;
    MAP = colormap('jet');
    lenmap = length(MAP(:,1));
    colscale=(norb/lenmap);
    for q=1:norb
        col=ceil(q/colscale);
        tvect=linspace((q-1)*T, q*T, pointsperorb);
        [r, v] = propagator(rin,vin,muE,tvect,J2,R_E,date0,Cr,Psr,Am);
        rin=r(end,:);
        vin=v(end,:);
        pos=[r(:,1) r(:,2) r(:,3)];
        switch j
            case 1
                plot3(pos(:,1),pos(:,2),pos(:,3),'Color',MAP(col,:),'LineWidth',3);
                set(gca,'FontSize',20)
            case 2
                plot3(pos(:,1),pos(:,2),pos(:,3),'Color',MAP(col,:));
                set(gca,'FontSize',20)
        end
        hold on
    end
    scatter3(r0(1), r0(2), r0(3),'LineWidth',3);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    switch j
        case 1
            title('Orbit representation','FontSize',20);
        case 2
            title('Orbit evolution','FontSize',20);
            D = colorbar('EastOutside');
            D.Ticks = linspace(1,lenmap,5);
            D.TickLabels = {'1','250','500','750','1000'};
            D.Label.String = 'number of orbits';
    end
    grid on
    C = imread('EarthTexture.jpg');
    [x, y, z] = ellipsoid(0,0,0, 6378.135, 6378.135,6356.750,1E2);
    surf(x,y,z,circshift(flip(C),[0,0]), 'FaceColor','texturemap','EdgeColor','none');
    axis equal;
    hold on;
end