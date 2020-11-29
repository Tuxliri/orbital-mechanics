function plot_groundtrack(lon,lat,col)
% Plotting routine for ground tracks
%  
% PROTOTYPE:
%   a = repeating_ground_track(m, k, mu, omega)
% 
% INPUT:
%   lon[n]      vector of longitudes      [rad]
%   lat[n]      vector of latitudes       [rad]
%   col[1]    color for the plot        [-]
%
% OUTPUT:
%   ---

% TRY TO FIX THE "HOLE" AT 0 LONGITUDE (you should set NaN when lon == 0 and lon
% 
% 
% eps = 1.e-3;
% 
% for j = 1:length(lon)
%     if lon(j) < eps || abs(lon(j) - 2*pi) < 1.e-3
%         lat(j) = NaN;
%         lon(j) = NaN;
%     end
% end

% wrap latitude and longitude to [-180 180] degrees

lon_wrapped = wrapTo180(rad2deg(lon));
lat_wrapped = wrapTo180(rad2deg(lat));

% FIX THE "HOLE 
for j = 2:length(lon_wrapped)
    if lon_wrapped(j)*lon_wrapped(j-1) < 0
        lat_wrapped(j) = NaN;
    end
end

plot(lon_wrapped,lat_wrapped,'Color', col, 'linewidth', 2)
hold on

% display a circle on the starting point
plot(lon_wrapped(1),lat_wrapped(1), 'Marker','o', 'Color', col, 'linewidth', 2)
hold on

% display a square on the final point

plot(lon_wrapped(end),lat_wrapped(end), 'Marker', 's', 'Color', col, 'linewidth', 2)
hold on

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

