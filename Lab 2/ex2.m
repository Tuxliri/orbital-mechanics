% Exercise 1: Computation of ground tracks
% 
% Plot the ground track for the following orbits:
% 1.
% ð?‘Ž = 8350 km, ð?‘’ = 0.1976, ð?‘– = 60 deg, Î© = 270 deg,
% ð?œ” = 45 deg, ð?‘“0 = 230 deg (taken from [1], Example 4.12)

clc
clear all
close all


%% Select the type of orbit
% 

orbit = 4;

switch orbit
    case 1
        a = 8350;
        e = 0.19760;
        i = deg2rad(60);
        RAAN = deg2rad(270);
        omega = deg2rad(45);
        f0 = deg2rad(230);
        periods = 15;               % number of periods to plot
        k =12;
        m = 1;
        
    case 2      % Molnyia orbit
        a = 26600;
        e = 0.74;
        i = deg2rad(63.4);
        RAAN = deg2rad(50);
        omega = deg2rad(280);
        f0 = deg2rad(0);
        periods = 30;               % number of periods to plot
        k = 2;
        m = 1;
        
    case 3      % Circular LEO orbit
        a = 7178.137;
        e = 0;
        i = [0 30 98];
        
        i = i(2);   % select which inclination you want
        
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
        omega = deg2rad(9.18854);
        f0 = deg2rad(0);
        periods = 5;               % number of periods to plot
        m = 1;
        k = 1;
  
end

% DATA

mu_earth = 398600;              % gravitational parameter of Earth      [ km^3/s^2 ]
omega_e = deg2rad(15.04)/3600;  % angular velocity of Earth's rotation  [ rad/s ]
R_earth = 6378.137;             % radius of the Earth                   [ km ]
gw_longitude0 = 0;              % longitude of greenwhich [ rad ]

%% Calculate the semi-major axis of the repeating ground track

a_new = repeating_ground_track(m,k,mu_earth,omega_e); 

%% Create the state vectors from orbital parameters

% calculate the orbital periods
T = 2*pi*sqrt(a^3/mu_earth);            % Original
T_new = 2*pi*sqrt(a_new^3/mu_earth);    % Repeating ground track


t = (0:1:periods*T);            % Original
t_new = (0:1:periods*T_new);    % Repeating ground track

state_vec = [a e i RAAN omega f0];                  % Original
state_vec_new = [a_new e i RAAN omega f0];           % Repeating ground track

 % Original
[alpha, delta, lon, lat] = groundtrack(state_vec, gw_longitude0, t, omega_e, mu_earth, 0);                         
 % Repeating ground track
[alpha_new, delta_new, lon_new, lat_new] = groundtrack(state_vec_new, gw_longitude0, t_new, omega_e, mu_earth, 0); 

 % Original
plot_groundtrack(lon,lat,'#77AC30');

 % Repeating ground track
plot_groundtrack(lon_new,lat_new, 'r' );

legend ('Ground track', 'Start', 'End','Repeating Ground track', 'Start', 'End', 'Orientation','horizontal','Location','northoutside' );
