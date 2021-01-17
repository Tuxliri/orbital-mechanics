function plot_groundtrack(lon,lat,col)
% Plotting routine for ground tracks
%  
% PROTOTYPE:
%   a = plot_groundtrack(lon,lat,col)
% 
% INPUT:
%   lon [ n ]      Vector of longitudes      [ rad ]
%   lat [ n ]      Vector of latitudes       [ rad ]
%   col [ 1 ]      Color for the plot        [ - ]
%
% OUTPUT:
%
% CONTRIBUTORS:
%   Davide Iafrate      
%   Alkady Marwan
%   Pedro Bossi Núñez
%   Davide Demartini
%
% VERSIONS:
%   2020-10-16: First version
%   2020-12-15: Second version
%   2021-01-16: Third version

% Input check
if nargin < 3
    col = 'r'
end

% Constrain latitude and longitude to be in [-180°,180°]
lon_wrapped = wrapTo180(rad2deg(lon));
lat_wrapped = wrapTo180(rad2deg(lat));

% Fill the "HOLE"
for j = 2:length(lon_wrapped)
    if lon_wrapped(j)*lon_wrapped(j-1) < 0
        lat_wrapped(j) = NaN;
    end
end

% Plot the ground track
plot(lon_wrapped,lat_wrapped,'Color', col, 'linewidth', 2)
hold on

% Highlight starting point (circle)
plot(lon_wrapped(1),lat_wrapped(1), 'Marker','o', 'Color', col, 'linewidth', 2)

% Highlight final point (square)
plot(lon_wrapped(end),lat_wrapped(end), 'Marker', 's', 'Color', col, 'linewidth', 2)

% Plot world map
axis([-180 180 -90 90])
xticks(-180:30:180);
yticks(-90:30:90);
I = imread('EarthTexture.jpg'); 
h = image(xlim,-ylim,I); 
uistack(h,'bottom')
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend ('Ground track', 'Start', 'End', 'Orientation','horizontal','Location','northoutside' )
end