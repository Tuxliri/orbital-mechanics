close all

%% Inputs
% Physical parameters

k = 10;        % spring constant                  [N/m]
c = 0.3;       %viscous damping coefficient      [N*s/m]
m = 1;         % mass of the block                [Kg]

% Oscillator parameters

omega0 = sqrt(k/m);
gamma = c/2*m;

%% Check if oscillator is underdamped

gamma^2 < omega0^2

% Initial conditions

y0 = [15; 0];

% Set time span

tspan = linspace(0, 10, 1000);

% Set options

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the analytical calculation

Y = an_harmonic_oscill(tspan, y0, omega0, gamma);

% Perform the numerical integration

[ T, Y1 ] = ode113( @(t,y) ode_harm_oscill(t,y,omega0,gamma), tspan, y0, options);

% Plot the results
%% Analytical solution
% Position Plot

figure()
subplot(2,1,1)
plot(tspan, Y(1,:), '-')
xlabel('Time [T]');
ylabel('Position [L]');
title('Position');

% Velocity 
subplot(2,1,2)
plot(tspan, Y(2,:), '-')
xlabel('Time [T]');
ylabel('Velocity [L/T]');
title('Velocity');

sgtitle('Analytical Solution')

%% Numerical Solution

figure()
subplot(2,1,1)
plot(tspan, Y(1,:), '-')
xlabel('Time [T]');
ylabel('Position [L]');
title('Position');

% Velocity 
subplot(2,1,2)
plot(tspan, Y1(:,2), '-')
xlabel('Time [T]');
ylabel('Velocity [L/T]');
title('Velocity');

sgtitle('Numerical Solution')