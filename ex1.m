% Exercise 1: Computation of ground tracks
% 
% Plot the ground track for the following orbits:
% 1.
% 𝑎 = 8350 km, 𝑒 = 0.1976, 𝑖 = 60 deg, Ω = 270 deg,
% 𝜔 = 45 deg, 𝑓0 = 230 deg (taken from [1], Example 4.12)

clc
clear all
close all


%% Select the type of orbit
% 

orbit = 1;

switch orbit
    case 1
        a = 8350;
        e = 0.19760;
        i = deg2rad(60);
        RAAN = deg2rad(270);
        omega = deg2rad(45);
        f0 = deg2rad(230);
        periods = 3.25;               % number of periods to plot
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
        
        i = i(1);   % select which inclination you want
        
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
                    
                
        i = deg2rad(0);
        RAAN = deg2rad(0);
        omega = deg2rad(40);
        f0 = deg2rad(0);
        periods = 5;               % number of periods to plot
        
        
    case 4      % GEO orbit
        a = 7178.137;
        e = 0;
        i = deg2rad(0);
        RAAN = deg2rad(0);
        omega = deg2rad(40);
        f0 = deg2rad(0);
        periods = 5;               % number of periods to plot
end

% DATA

mu_earth = 398600;              % gravitational parameter of Earth      [ km^3/s^2 ]
omega_e = deg2rad(15.04)/3600;  % angular velocity of Earth's rotation  [ rad/s ]
R_earth = 6378.137;             % radius of the Earth                   [ km ]
gw_longitude0 = 0;              % longitude of greenwhich [ rad ]

%% Create the state vector from orbital parameters

% calculate the orbital period
T = 2*pi*sqrt(a^3/mu_earth);


t = (0:1:periods*T);


    state_vec = [a e i RAAN omega f0];
    
    [alpha, delta, lon, lat] = groundtrack(state_vec, gw_longitude0, t, omega_e, mu_earth, 0); 

    plot_groundtrack(lon,lat);
