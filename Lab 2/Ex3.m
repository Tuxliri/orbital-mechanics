% Exercise 3: Computation of ground tracks
% 


clc
clear all
close all


%% Select the type of orbit
% 

orbit = 2;

switch orbit
    case 1
        a = 8350;
        e = 0.19760;
        i = deg2rad(60);
        RAAN = deg2rad(270);
        omega = deg2rad(45);
        f0 = deg2rad(230);
        periods = 1;               % number of periods to plot
        k =12;
        m = 1;
        
    case 2      % Molnyia orbit
        a = 26600;
        e = 0.74;
        i = deg2rad(63.4);
        RAAN = deg2rad(50);
        omega = deg2rad(280);
        f0 = deg2rad(0);
        periods = 1;               % number of periods to plot
        k = 1;
        m = 1;
        
    case 3      % Circular LEO orbit
        a = 7178.137;
        e = 0;
        i = [0 30 98];
        
        i = i(3);   % select which inclination you want
        
        switch i    % select the number of revolutions k and earth rotations m
                    % based on the inclination chosen for the orbit
            
            case 0
                k=20;
                m=2;
                
            case 30
                k=29;
                m=2;
                
            case 98
                k=15;
                m=1;
        end
                    
                
        i = deg2rad(i);
        RAAN = deg2rad(0);
        omega = deg2rad(40);
        f0 = deg2rad(0);
        periods = 30;               % number of periods to plot
        
        
    case 4      % GEO orbit
        a = 42164;
        e = 0;
        i = 0;
        RAAN = deg2rad(0);
        omega = deg2rad(0);
        f0 = deg2rad(20);
        periods = 500;               % number of periods to plot
        m = 1;
        k = 1;
        
end

% DATA

mu_E = 398600;                  % gravitational parameter of Earth      [ km^3/s^2 ]
omega_e = deg2rad(15.04)/3600;  % angular velocity of Earth's rotation  [ rad/s ]
R_E = 6378.137;                 % radius of the Earth                   [ km ]
gw_longitude0 = 0;              % longitude of greenwhich               [ rad ]
J2 = 0.00108263;                % second zonal armonic of earth         [ - ]

%% Calculate the semi-major axis of the repeating ground track

a_new = repeating_ground_track(m,k,mu_E,omega_e); 

%% Calculate the semi-major axis of the perturbed repeating ground track

a_secular = repeating_ground_track_J2(m, k, mu_E, omega_e, J2, R_E, e, i);

%% Create the state vectors from orbital parameters

% calculate the orbital periods
T = 2*pi*sqrt(a^3/mu_E);                    % Original
T_new = 2*pi*sqrt(a_new^3/mu_E);            % Repeating ground track
T_secular = 2*pi*sqrt(a_secular^3/mu_E);    % Repeating secular ground track for the perturbed secular orbit

%% Create time vectors

t = (0:1:periods*T);                    % Original
t_new = (0:1:periods*T_new);            % Repeating ground track
t_secular = (0:1:periods*T_secular);    % Repeating secular ground track for the perturbed secular orbit

%% create state vectors

state_vec = [a e i RAAN omega f0];                   % Original
state_vec_new = [a_new e i RAAN omega f0];           % Repeating ground track
state_vec_secular = [a_secular e i RAAN omega f0];   % Repeating ground track for the perturbed secular orbit

%% Compute RA, declination, lon and latitude
% Original
[alpha, delta, lon, lat] = groundtrack(state_vec, gw_longitude0, t, omega_e, mu_E, 0);   

% Repeating ground track
[alpha_new, delta_new, lon_new, lat_new] = groundtrack(state_vec_new, gw_longitude0, t_new, omega_e, mu_E, 0); 

% Repeating perturbed ground track
[alpha_secular, delta_secular, lon_secular, lat_secular] = groundtrack_J2(state_vec_secular,gw_longitude0,t_secular,omega_e, mu_E,0, J2, R_E);

%% Plotting of the groundtracks

% Original
plot_groundtrack(lon,lat,'#77AC30');

 % Repeating ground track
plot_groundtrack(lon_new,lat_new, 'r' );

% Perturbed secular ground track
% plot_groundtrack(lon_secular, lat_secular, 'y');

legend ('Original', 'Start', 'End','Repeating Ground track', 'Start', 'End', 'Repeating Secular Ground track', 'Start', 'End', 'Orientation','horizontal','Location','northoutside' );
