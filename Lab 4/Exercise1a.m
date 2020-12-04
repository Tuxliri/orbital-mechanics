clc
clearvars
close all

%% Set path to default
path(pathdef);
% Add [...] folder to path
addpath(genpath('../functions\'));

addpath(genpath('../Common\'));

%% useful constants
muSun = astroConstants(4);
muE = astroConstants(13);
AU = astroConstants(2);

%% Define Data
rE = [1; 0; 0] * AU;
radius_e = astroConstants(23);      % Earth mean radius

% Compute planet velocity
VEmag = sqrt(muSun/norm(rE));

VE = VEmag .* cross([0;0;1],rE)/norm(rE);

Delta = 9200;       % Asymptote distance [km]

% Select the appropriate case for the asymptote position
asymptote = 'under';

switch asymptote
    
    case 'front'    % leading-side flyby
        
        % In the case of the asymptote in front of the planet, the axis of
        % rotation of vinf u=[0; 0; -1] is parallel to z but in the
        % opposite direction
                       
        R = @(x) [cos(-x) -sin(-x) 0
                  sin(-x)  cos(-x) 0
                  0         0    1];
        
    case 'behind'   % trailing-side flyby
        
        % In the case of the asymptote behind the planet, the axis of
        % rotation of vinf u=[0; 0; 1] is parallel to z
 
        
        R = @(x) [cos(x) -sin(x) 0
                  sin(x)  cos(x) 0
                  0         0    1];
        
    case 'under'
        
         % In the case of the asymptote under the planet, the axis of
         % rotation of vinf u=[0; -1; 0] is parallel to y and in the
         % opposite direction
         
         R = @(x) [cos(-x) 0 sin(-x)
                   0       1      0
                   -sin(-x) 0 cos(-x)];
       
end

%% Exercise 1a
vinfM = [15.1; 0; 0]; % Incoming velocity in planetocentric frame [km/s]

[rp,a,e,delta] = flybyunpow(vinfM,muE,Delta);

vinfP = R(delta)*vinfM;
dv = vinfP - vinfM;

V_minus = VE + vinfM;
V_plus = VE + vinfP;

%% Propagate the heliocentric trajectories
tspanA = linspace(0,-20000000,10000);
[r_A, v_A] = propagator(rE,V_minus,muSun,tspanA);

tspanB = linspace(0,20000000,10000);
[r_B, v_B] = propagator(rE,V_plus,muSun,tspanB);

%% Plotting

figure(1)
hold on
legend on
grid on
axis equal

SUN = plot3(0,0,0,'-o', 'MarkerSize',10);
SUN.MarkerFaceColor = '#F7CB45'; 
SUN.MarkerEdgeColor = SUN.MarkerFaceColor;
SUN.HandleVisibility = 'off';

rE = rE./AU;
EARTH = plot3(rE(1),rE(2),rE(3),'-o', 'MarkerSize',5);
EARTH.MarkerFaceColor = '#285FF4'; 
EARTH.MarkerEdgeColor = EARTH.MarkerFaceColor;
EARTH.HandleVisibility = 'off';

r_A =r_A./AU;
r_B = r_B./AU;

ARR_arc = plot3(r_A(:,1),r_A(:,2),r_A(:,3),'LineWidth',2);
ARR_arc.DisplayName = 'Before flyby';

DEP_arc = plot3(r_B(:,1),r_B(:,2),r_B(:,3),'LineWidth',2);
DEP_arc.DisplayName = 'After flyby';

xlabel('$x \left [ AU \right ]$','Interpreter','latex')
ylabel('$y \left [ AU \right ]$','Interpreter','latex')
zlabel('$z \left [ AU \right ]$','Interpreter','latex')