%% Exercise 3b - Plot f(t) for different orbits

clc
clear all
close all

%% Problem's parameter
a = 7000;                       % Semi latus rectum               [ km ]
mu = 398600;                    % Earth's gravitational parameter [ km^3/s^2 ]
e = [0 0.2 0.4 0.6 0.8 0.95];   % Eccentricities                  [ - ]
t0 = 0;                         % Initial time                    [ s ]
f0 = 0;                         % Initial true anomaly            [ deg ]
k = 2;                          % Number of periods to plot       [ - ]
N = 150;                        % Number of points to calculate   [ - ]

f = zeros(length(e),N);
for i=1:length(e)
    f(i,:) = true_anomaly(e(i), a, mu, N, k, t0,f0);
end


%% Plotting 2d diagram of true anomaly

figure()

for i= 1:length(e)
    plot(linspace(t0,length(f(1,:)),N),f(i,:),'LineWidth', 1, 'DisplayName', strcat('e= ',num2str(e(i))))
    hold on
end

legend
title('True Anomaly f')
ylabel('f [deg]')
xlabel('t [s]')
grid on

%% Plotting 3d diagram of true anomaly

figure()
[X,Y] = meshgrid(linspace(t0,length(f(1,:)),N),e);
surf(X,Y,f)
xlabel('t [s]')
ylabel('e [-]')
zlabel('f [deg]')