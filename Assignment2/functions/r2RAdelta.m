function [delta,alpha] = r2RAdelta(r)
% r2RAdelta function that converts the position vector into the
%   rightascention and delta  
%   
% PROTOTYPE:
%    [delta, alpha] = r2RAdelta(r)
% 
% INPUT:
%   r     [ 3 x 1 ]    Initial radius vector             [ km ]
%
% OUTPUT:
%   alpha [ n x 1 ]    Right ascension vector            [ rad ]
%   delta [ n x 1 ]    Declination vector                [ rad ]
%
% CONTRIBUTORS:
%   Davide Demartini
%   Davide Iafrate      
%   Alkady Marwan
%   Pedro Bossi Núñez
%
% VERSIONS
%   15-12-2020: First version

% Set void vectors
r_mag = zeros(1, length(r(:,3)));
delta = zeros(1, length(r(:,3)));
alpha = zeros(1, length(r(:,3)));

% Declination
for i = 1:length(r(:,3))
    r_mag(i) = norm(r(i,:));
    delta(i) = asin(r(i,3)/r_mag(i));
end

% Right ascension
for i = 1:length(r(:,3))
    if r(i,2) > 0
        alpha(i) = acos((r(i,1)/r_mag(i))/cos(delta(i)));
    else
        alpha(i) = 2*pi - acos((r(i,1)/r_mag(i))/cos(delta(i)));
    end
end
end