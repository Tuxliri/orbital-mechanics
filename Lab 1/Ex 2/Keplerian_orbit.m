%% Exercise 2a -  Integrate numerically a Keplerian orbit (two-body problem)
clc
clear all
close all

%% Set options for ODE solver
opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','on');

%% Set initial conditions

mu = 398600;                % Earth's gravitational parameter   [ km^3/s^2 ]
r0 = [ -7128.137, 0, 0 ];   % Initial radius vector             [ km ]
v0 = [ 0, -9.781, 0 ];      % Initial velocity vector           [ km/s ]

%% Calculate theoretical constants
% Specific energy (constant)                                    [ km^2/s^2 ]
epsilon0 = 0.5*norm(v0)^2 - mu/norm(r0);

a = -mu/(2*epsilon0);        % Semi-major axis                   [ km ]

y0 = [ r0(1) v0(1) r0(2) v0(2) r0(3) v0(3) ];       % Initial state vector

% Define time of integration
T = 2*pi*sqrt(a^3/mu);      % Period                            [ s ]
year = 31536000;            % 1 year period                     [ s ]
tspan = (0:100:T);

%% Compute the integration of the ode

[t, y] = ode113(@(t,y) twobodyode(t,y,mu), tspan, y0, opts);

%% Calculate numerically the velocities, radiuses and specific energy

v = [y(:,2) y(:,4) y(:,6)];                   % Velocities       [ km/s ]
r = [y(:,1) y(:,3) y(:,5)];                   % Radiuses         [ km ]

epsilon = zeros(1,length(v));
for i=1:length(v)
    epsilon(i) = 0.5*(norm(v(i,:)))^2 - mu/norm(r(i,:));  % Specific Energy [km^2/s^2]
end

%% Calculate angular momentum vector h and eccentricity vector and check that
%  they remain constant in magnitude and direction

h = zeros(length(r),3);
e = zeros(length(r),3);
e_abs = zeros(1,length(v));
h_abs = zeros(1,length(v));
theta = zeros(1,length(v));
for i=1:length(r)
    h(i,:) = cross(r(i,:),v(i,:));           % Angular momentum (constant)       [ km^2/s ]
    e(i,:) = cross(v(i,:) , h(i,:))./mu - r(i,:)./norm(r(i,:)); % Eccentricity vector
    e_abs(i) = norm(e(i,:));
    h_abs(i) = norm(h(i,:));
    theta(i) = acos(dot(r(i,:),e(i,:))/(norm(r(i,:))*norm(e(i,:))));  % True anomaly    [-]
end

%% Radial and transversal components of the velocity

v_r = zeros(1,length(v));
v_theta = zeros(1,length(v));

for i=1:length(e)
    v_r(i) = dot(v(i,:),r(i,:))/norm(r(i,:));
    v_theta(i) = mu*(1 + norm(e(i,:))*cos(theta(i)))/norm(h(i,:));
end

%% Perpendicularity check e - h
% Vectors are perpendicular if their dot product is 0

perpendicularity = zeros(1,length(v));

for i = 1:length(e)
    perpendicularity(i) = dot(e(i,:),h(i,:));
end

%% Plot specific energy
figure('Name','Specific Energy','NumberTitle','off');
plot(t,epsilon,'Color', 'b','LineWidth',2)
grid on
title('Specific Energy')
ylabel('$\mathbf{\epsilon [km^2/s^2]}$','Interpreter','latex')
xlabel('t [s]')

%% Plot h - e perpendicularity

figure('Name','Orbit Properties','NumberTitle','off');
subplot(2,2,1)
plot(t, perpendicularity,'Color', 'b','LineWidth',2)
grid on
title('h - e dot product')
ylabel('$\mathbf{h}\cdot \mathbf{e}\;\bf[km^2/s]$','Interpreter','latex','Fontsize', 14)
xlabel('t [s]')

%% Plot Angular momentum h
subplot(2,2,2)
plot(t, h_abs, '-', t, h(:,1), '--', t, h(:,2),':', t, h(:,3), '-.', 'LineWidth',2)
grid on
title('Angular momentum')
ylabel('||h||, hx, hy, hz [km^2/s]')
xlabel('t [s]')
legend('||h||', 'hx', 'hy', 'hz')

%% Plot eccentricity e

subplot(2,2,3)
plot(t, e_abs, '-', t, e(:,1), '--', t, e(:,2),':', t, e(:,3), '-.', 'LineWidth',2)
grid on
title('Eccentricity vector')
ylabel('||e||, ex, ey, ez [-]')
xlabel('t [s]')
legend('||e||', 'ex', 'ey', 'ez')

%% Plot radial and transversal velocity

subplot(2,2,4)
plot(t, v_theta, 'LineWidth', 2)
grid on
hold on
plot(t,v_r, 'LineWidth', 2)
title('Radial and transversal velocity')
ylabel('$\mathbf{v_{r},v_{\theta }\;[km/s]}$','Interpreter', 'latex','Fontsize', 14)
xlabel('$ t $','Interpreter', 'latex', 'Fontsize', 14)
legend('$v_{\theta}$','$v_{r}$','Interpreter', 'latex', 'Fontsize', 12)

%% Orbit 3D plot

figure('Name','Orbit Plot','NumberTitle','off');
plot3(y(:,1),y(:,3),y(:,5), 'Color', 'r', 'LineWidth', 2)
axis equal
axis([-inf inf -inf inf -20000 20000])
xlabel('$\mathbf{R_{x} [Km]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('$\mathbf{R_{y} [Km]}$','Interpreter', 'latex','Fontsize', 14)
zlabel('$\mathbf{R_{z} [Km]}$','Interpreter', 'latex','Fontsize', 14)
grid on
hold on
plot3(0,0,0, '-o','Color','b','MarkerSize',20,'MarkerFaceColor','#D9FFFF')
