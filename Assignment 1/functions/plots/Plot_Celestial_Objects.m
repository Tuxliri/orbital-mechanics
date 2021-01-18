function Plot_Celestial_Objects(sun,obj,orbit,date,scale)
%
% DESCRIPTION:
%   Given a date, this algorithm plots the position of whished planets
%   around the sun whith their respective orbits. See PlotObject for more
%   info.
%
% PROTOTYPE:
%   Plot_Celestial_Objects(sun, obj, orbit, date, scale)
%
% INPUT:
%   sun   [ 1 ]: diplsy Sun
%           0 Don't display
%           1 Display
%   obj   [ n x 1 ]: vector containing the numbers of the celestial objects
%       that must be displayed.
%           1 Mercury
%           2 Venus
%           3 Earth
%           4 Mars
%           5 Jupiter
%           6 Saturn
%           7 Uranus
%           8 Neptune
%   orbit [ 1 ]: numerical value used to select what to dysplay.
%           0 Celestial object only
%           1 Celestial object's orbit only
%           2 Celestial object and it's orbit
%   date  [6 x 1]: date in the Gregorian calendar, as a 6-elements vector
%       [year, month, day, hour, minute, second]. For dates before
%       1582, the resulting date components are valid only in the Gregorian
%       proleptic calendar. This is based on the Gregorian calendar but 
%       extended to cover dates before its introduction. Date must be after
%       12:00 noon, 24 November -4713. 
%   scale [ 9 x 1 ]: how big should a planet be displayed (experiment with 
%       different data). The 1st position is for Mercury, the 9th is for
%       the Sun
%
% OUTPUT:
%   Graphical representation of selected celestial objects at given date
%
% CALLED FUNCTIONS:
%   date2mjd2000
%   uplanet
%   PlotObject
%   astroConstants
%   plotOrbit
%   kep2car
%
% CONTRIBUTORS:
%   Alkady Marwan
%   Bossi Nunez Pedro
%   Davide Demartini
%   Davide Iafrate
%
% VERSIONS:
%   20-12-2020: First version
%   18-01-2021: Second Version
%

% Convert Gregorian (current) date into MJD200
mjd = date2mjd2000(date);

% Extract ephemeresis information for given Celestial Objects
a = zeros(length(obj),1);
e = zeros(length(obj),1);
i = zeros(length(obj),1);
Om = zeros(length(obj),1);
om = zeros(length(obj),1);
theta = zeros(length(obj),1);
for j = 1:length(obj)
    [kep, ~] = uplanet(mjd, obj(j));
    a(j) = kep(1);
    e(j) = kep(2);
    i(j) = kep(3);
    Om(j) = kep(4);
    om(j) = kep(5);
    theta(j) = kep(6);
end

% Calculate CO positions
r = zeros(length(obj),3); % r [ n x 3]
mu = astroConstants(4);
for j = 1:length(a)
    [r_, ~] = kep2car(a(j), e(j), i(j), Om(j), om(j), theta(j), mu); % r [ n x 3 ]
    r(j,1) = r_(1);
    r(j,2) = r_(2);
    r(j,3) = r_(3);
end

% Display COs or orbits or both
if orbit == 0
    for j = 1:length(obj)
        PlotObject(obj(j), r(j,:), scale(j));
        hold on;
    end
elseif orbit == 1
    for j = 1:length(obj)
        plotOrbit(obj(j), date);
        hold on;
    end
elseif orbit == 2 
    % The two for cycles must be separated in order to have meaningful
    % legend when plotting
    for j = 1:length(obj)
        plotOrbit(obj(j), date);
        hold on;
    end
    for i = 1:length(obj)
        PlotObject(obj(i), r(i,:), scale(i));
        hold on;
    end
else
    error('Invalid orbit value');
end

% Display Sun
if sun == 1
    if length(obj) < 9
        j = length(obj);
    end
    PlotObject(0, [0 0 0], scale(j+1));
end
end