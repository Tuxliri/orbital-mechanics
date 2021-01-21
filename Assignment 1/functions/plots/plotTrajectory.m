function plotTrajectory(r1,r2,v1,v2,mu,orbit)
%
% DESCRIPTION:
%   Given the initial and final position and velocity for each arc that
%   composes the trajectory, this algorithm displays the trajectory and the
%   orbits it is composed of.
%
% PROTOTYPE:
%   plotTrajectory(r1, r2, v1, v2, mu, orbit)
%
% INPUT:
%   r1    [ n x 1 ] position at the beginning of each arc that composes the
%                   trajectory [ km ]
%   r2    [ n x 1 ] position at the end of each arc that composes the
%                   trajectory [ km ]
%   v1    [ n x 1 ] velocity at the beginning of each arc that composes the
%                   trajectory [ km/s ]
%   v2    [ n x 1 ] velocity at the end of each arc that composes the
%                   trajectory [ km/s ]
%   mu    [ 1 ]     planetary constant of the celestial body the probe is
%                   orbiting   [ km^3/s^2 ]
%   orbit [ 1 ]     0 displays the trajectory only, if 1 displays orbits as
%                   well       [ - ]
%
% OUTPUT:
%   Graphical representation of probe's trajectory
%
% CALLED FUNCTIONS:
%   car2kep
%   kep2car
%
% CONTRIBUTORS:
%   Alkady Marwan
%   Bossi Nunez Pedro
%   Davide Demartini
%   Davide Iafrate
%
% VERSIONS:
%   18-01-2021: First version
%

% Determine how many arcs the trajectory is composed of
[m, n] = size(r1);

% Empty vectors for keplerian elements
a = zeros(1,m); e = zeros(1,m); i = zeros(1,m); RAAN = zeros(1,m); om = zeros(1,m); f1 = zeros(1,m); f2 = zeros(1,m);

% Determine keplerian elements + initial & final position for each arc
for j = 1:m
    [a(1,j), e(1,j), i(1,j), RAAN(1,j), om(1,j), f1(1,j)] = car2kep(r1(j,:), v1(j,:), mu);
    [~, ~, ~, ~, ~, f2(1,j)] = car2kep(r2(j,:), v2(j,:), mu);
end

% Display the entire orbit
if orbit == 1
    for j = 1:m
        color = '--';
        plotOrbit(a(j), e(j), i(j), RAAN(j), om(j), mu, color)
        hold on;
    end
end


% Determine position for an array of f (for each arc)
for j = 1:m % Select arc
    x = []; y  = []; z = [];
    if (f2(j) - f1(j)) < 0
        f2(j) = f2(j) + 2*pi;
    end
    for  f = [f1(j):0.001:f2(j)] % Create array of f
        [r_vec, ~] = kep2car(a(j), e(j), i(j), RAAN(j), om(j), f, mu);
        x = [x r_vec(1)];
        y = [y r_vec(2)];
        z = [z r_vec(3)];
    end
    plot3(x(1,:), y(1,:), z(1,:),'LineWidth',2)
    hold on;
end
end