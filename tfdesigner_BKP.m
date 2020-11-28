function [DV, T1, T2, DV1, DV2,dep_orbit, arr_orbit,rtf0,vtf0]...
    = tfdesigner(planetA,planetB,windowA,windowB,ndep,narr)
% Get the different DV needed for the transfer between two celestial bodies
% and the respective times of departure and arrival
% PROTOTYPE:
%    [DV, T1, T2, DV1, DV2, dep_orbit, arr_orbit, TFarc, mu] = tfdesigner(planetA,planetB,windowA,windowB)
%
%  INPUT :
%	planetA[1]    Integer number identifying the departure celestial body (< 11)
%                   1:   Mercury
%                   2:   Venus
%                   3:   Earth
%                   4:   Mars
%                   5:   Jupiter
%                   6:   Saturn
%                   7:   Uranus
%                   8:   Neptune
%                   9:   Pluto
%                   10:  Sun
%
%	planetB[1]    Integer number identifying the departure celestial body (< 11)
%   windowA[2x6]    Vector with transfer window opening and closing times in
%                 DATE from departure body; first row is opening of
%                 transfer window, second row is closing of transfer window
%   windowB[2x6]    Vector with transfer window opening and closing times in
%                 DATE to arrival body; first row is opening of
%                 transfer window, second row is closing of transfer window
%   ndep [1]     number of points for which to get the position vector of
%                the departure orbit
%   narr[1]      number of points for which to calculate the position
%                vector of the arrival orbit
%
%
% OUTPUT:
%   DV[N,M]     total delta V for transfer   [km/s]
%   T1[N]       departure window times MJD2000 [days]
%   T2[M]       arrvial window times MJD2000 [s]
%   DV1[N,M]    cost of departure maneuvre   [km/s]
%   DV2[N,M]    cost of arrival maneuvre     [km/s]
%   dep_orbit[ndep,3] position vector of departure orbit [km]
%   arr_orbit[narr,3] position vector of arrival orbit [km]
%
%   TF_arc     position vector of transfer arc[km]
%   mu          planetary constant of central body [km^3/s^2]
%
% CONTRIBUTORS:
%
%   Davide Iafrate
%
% VERSIONS
%   2020-11-23: First version

narginchk(6,6);

DAY2SECS = 24*3600;
mu = astroConstants(4);

% Create the vector of dates for which to get the keplerian elements of the
% planetary bodies

A_dep_start = date2mjd2000(windowA(1,1:6));
A_dep_end = date2mjd2000(windowA(2,1:6));

B_arr_start = date2mjd2000(windowB(1,1:6));
B_arr_end = date2mjd2000(windowB(2,1:6));

dep_window = linspace(A_dep_start,A_dep_end,ndep);
arr_window = linspace(B_arr_start,B_arr_end,narr);

% get keplerian elements, positions and velocities of planet A and B

[kepA,RA,VA] = ephemeris(dep_window,planetA);
[kepB,RB,VB] = ephemeris(arr_window,planetB);
dep_orbit = [RA; VA];
arr_orbit = [RB; VB];

%%
DV = zeros(ndep,narr);
DV1 = zeros(ndep,narr);
DV2 = zeros(ndep,narr);

rtf0 = zeros(ndep,narr,3);
vtf0 = zeros(ndep,narr,3);
for i=1:ndep
    for k=1:narr
        %% Planet A - departure planet
        
        % positions and velocities along the whole departure window
        ri = RA(i,1:3);
        vi = VA(i,1:3);
        
        t1 = dep_window(i)*DAY2SECS;    % convert to seconds
        
        if arr_window(k) < dep_window(i)
            DV(i,k) = NaN;
            DV1(i,k) = NaN;
            DV2(i,k)= NaN;
            
        else
            %% Planet B - departure planet
            
            t2 = arr_window(k)*DAY2SECS;    % convert to seconds
            
            % final positions and velocities along the whole
            % arrival window
            
            rf = RB(k,1:3);
            vf = VB(k,1:3);
            
            %% Compute the total DeltaV for each departure time T1 and T2
            % NB tfixed_design function REQUIRES the date in seconds
            
            [DV(i,k), ~, ~,DV1(i,k),DV2(i,k), rtf0(i,k,1:3), vtf0(i,k,1:3)] = tfixed_design(ri,vi,rf,vf,t1,t2,mu);
            
        end
    end
end

T1 = dep_window;
T2 = arr_window;