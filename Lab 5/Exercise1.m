clc
clearvars
close all

%% Set path to default
path(pathdef);
% Add [...] folder to path
addpath(genpath('../functions_LABS\'));

addpath(genpath('../Common\'));

%% useful constants
const = astroConstants([23 13 9]);
R_e = const(1);
muE = const(2);
J2 = const(3);

% Initial conditions
kep0 = [7571 0.01 deg2rad(87.9) pi pi 0];
t0 = 0;
Torb = 2*pi/sqrt(muE/kep0(1)^3);
tspan = t0:1:100*Torb;

% Get the perturbed orbital elements
[t_gauss,kep_gauss] = ORBITPROPAGATOR(t0,kep0,tspan);

% Extract initial orbital elements
a0 = kep0(1);
e0 = kep0(2);
i0 = kep0(3);
RAAN0 = kep0(4);
omega0 = kep0(5);
f0 = kep0(6);

% Initial cartesian elements
[r0, v0] =  kep2car(a0,e0,i0,RAAN0,omega0,f0,muE);
y0 = [r0(1); v0(1); r0(2); v0(2); r0(3); v0(3)];

%% Solve with the cartesian elements
opts = odeset('Reltol',1e-16,'Stats','on');
[tB, y] = ode113(@(t,y) twobodyode_j2(t,y,muE,J2,R_e), tspan,y0,opts);

% Convert cartesian to keplerian
kepB = zeros(length(tB),6);
for i = 1:length(tB)
    ri = [y(i,1),y(i,3),y(i,5)];
    vi = [y(i,2),y(i,4),y(i,6)];
    [kepB(i,1),kepB(i,2),kepB(i,3),kepB(i,4),kepB(i,5),kepB(i,6)] = car2kep(ri,vi,muE);
end

%% Filtering lower frequencies

% Cut-off period
Tfilter = 3*Torb;

% Number of points for the filtering window
nwindow = nearest( Tfilter / (sum(diff(t_gauss)) / (numel(t_gauss)-1) ) );

% Filter elements ( no unwrapping)
kep_filtered = movmean( kep_gauss,nwindow,1);

%% Compare the two solutions through plotting
figure(1)

% Semi-major axis
subplot(2,3,1)
plot(tspan./Torb,kep_gauss(:,1),tspan./Torb,kepB(:,1),tspan./Torb,kep_filtered(:,1));
legend('Gauss equations','Cartesian','Secular (filtered)')
grid on
xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('$\mathbf{a [Km]}$','Interpreter', 'latex','Fontsize', 14)

% Eccentricity
subplot(2,3,2)
plot(tspan./Torb,kep_gauss(:,2),tspan./Torb,kepB(:,2),tspan./Torb,kep_filtered(:,2));
legend('Gauss equations','Cartesian','Secular (filtered)')
grid on
xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('$\mathbf{e [-]}$','Interpreter', 'latex','Fontsize', 14)


% inclination
subplot(2,3,3)
plot(tspan./Torb,rad2deg(kep_gauss(:,3)),tspan./Torb,rad2deg(kepB(:,3)),tspan./Torb,rad2deg(kep_filtered(:,3)));
legend('Gauss equations','Cartesian','Secular (filtered)')
grid on
xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('$\mathbf{i [deg]}$','Interpreter', 'latex','Fontsize', 14)

% RAAN
subplot(2,3,4)
plot(tspan./Torb,rad2deg(kep_gauss(:,4)),tspan./Torb,rad2deg(kepB(:,4)),tspan./Torb,rad2deg(kep_filtered(:,4)));
legend('Gauss equations','Cartesian','Secular (filtered)')
grid on
xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('$\mathbf{\Omega  [deg]}$','Interpreter', 'latex','Fontsize', 14)

% omega
subplot(2,3,5)
plot(tspan./Torb,rad2deg(kep_gauss(:,5)),tspan./Torb,rad2deg(kepB(:,5)),tspan./Torb,rad2deg(kep_filtered(:,5)));
legend('Gauss equations','Cartesian','Secular (filtered)')
grid on
xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('$\mathbf{\omega  [deg]}$','Interpreter', 'latex','Fontsize', 14)


% f
subplot(2,3,6)
plot(tspan./Torb,rad2deg(kep_gauss(:,6)),tspan./Torb,unwrap(rad2deg(kepB(:,6))),tspan./Torb,rad2deg(kep_filtered(:,6)));
legend('Gauss equations','Cartesian','Secular (filtered)')
grid on
xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('$\mathbf{f  [deg]}$','Interpreter', 'latex','Fontsize', 14)

%% Compare the difference between the two solutions through plotting
figure(2)

% Semi-major axis
subplot(2,3,1)
semilogy(tspan./Torb,abs((kep_gauss(:,1) - kepB(:,1)))./kep0(1));
grid on
xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('${|a_{Car} - a_{Gauss}| / a0 [-]}$','Interpreter', 'latex')

% Eccentricity
subplot(2,3,2)
semilogy(tspan./Torb,abs(kep_gauss(:,2) - kepB(:,2)));
grid on
xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('${|e_{Car} - e_{Gauss}| [-]}$','Interpreter', 'latex')


% inclination
subplot(2,3,3)
semilogy(tspan./Torb,abs(kep_gauss(:,3) - kepB(:,3)) ./(2*pi));
grid on
xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('${|i_{Car} - i_{Gauss}| / 2\pi [-]}$','Interpreter', 'latex')

% RAAN
subplot(2,3,4)
semilogy(tspan./Torb,abs(kep_gauss(:,4) - kepB(:,4)) ./(2*pi));
grid on
xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('${|\Omega_{Car} - \Omega_{Gauss}| / 2\pi [-]}$','Interpreter', 'latex')

% omega
subplot(2,3,5)
semilogy(tspan./Torb,abs(kep_gauss(:,5) - kepB(:,5)) ./(2*pi));
grid on
xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('${|\omega_{Car} - \omega_{Gauss}| / 2\pi [-]}$','Interpreter', 'latex')


% omega
subplot(2,3,6)
semilogy(tspan./Torb,abs(kep_gauss(:,6) - unwrap(kepB(:,6))) ./ kep0(5));
grid on
xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
ylabel('${|f_{Car} - f_{Gauss}| / f_{Gauss} [-]}$','Interpreter', 'latex')
