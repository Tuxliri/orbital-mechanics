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

Delta = (9200:1000:13200);       % Asymptote distance [km]

% Select the appropriate case for the asymptote position
asymptote = 'front';

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

%% Exercise 1b
vinfM = [15.1; 0; 0]; % Incoming velocity in planetocentric frame [km/s]

%% Plotting

figure(1)
hold on
legend on
grid on
axis equal


%% Propagate the heliocentric trajectories
V_minus = VE + vinfM;

tspanA = linspace(0,-15000000,10000);
[r_A, v_A] = propagator(rE,V_minus,muSun,tspanA);

r_A =r_A./AU;   % Scale to AU
ARR_arc = plot3(r_A(:,1),r_A(:,2),r_A(:,3),'LineWidth',2);
ARR_arc.DisplayName = 'Before flyby';

for i=1:length(Delta)
[rp(i),a(i),e(i),delta(i)] = flybyunpow(vinfM,muE,Delta(i));

vinfP(:,i) = R(delta(i))*vinfM;
dv = vinfP(:,i) - vinfM;


V_plus(:,i) = VE + vinfP(:,i);


tspanB = linspace(0,15000000,10000);
[r_B, v_B] = propagator(rE,V_plus(:,i),muSun,tspanB);


r_B = r_B./AU;
str=sprintf('Delta = %d', Delta(i));
plt(i)=plot3(r_B(:,1),r_B(:,2),r_B(:,3),'LineWidth',2);
plr(i).displayName = str;
end


SUN = plot3(0,0,0,'-o', 'MarkerSize',10);
SUN.MarkerFaceColor = '#F7CB45'; 
SUN.MarkerEdgeColor = SUN.MarkerFaceColor;
SUN.HandleVisibility = 'off';

rE = rE./AU;
EARTH = plot3(rE(1),rE(2),rE(3),'-o', 'MarkerSize',5);
EARTH.MarkerFaceColor = '#285FF4'; 
EARTH.MarkerEdgeColor = EARTH.MarkerFaceColor;
EARTH.HandleVisibility = 'off';

xlabel('$x \left [ AU \right ]$','Interpreter','latex')
ylabel('$y \left [ AU \right ]$','Interpreter','latex')
zlabel('$z \left [ AU \right ]$','Interpreter','latex')
