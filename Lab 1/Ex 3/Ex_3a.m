%% Exercise 3a - Implement a solver for Kepler’s equation

clc
clear all
close all

%% Problem's parameter
a = 7000;           % Semi latus rectum         [ km ]
mu = 398600;        % Earth's gravitational parameter [ km^3/s^2 ]
e = 0.95;            % Eccentricity                    [ - ]
t0 = 0;             % Initial time                    [ s ]
f0 =0;              % Initial true anomaly            [ deg ]

% Calculate period

T = 2*pi*sqrt(a^3/mu);      % Period                            [ s ]

t = (t0:10:2*T);

f = zeros(1,length(t));
for i=1:length(t)
    f(i) = fsolve_kepler_eq_solver(t(i),e,a,mu,0,0);
end

plot(t,f,'LineWidth', 2);
title('True Anomaly f')
ylabel('f [deg]')
xlabel('t [s]')
grid on