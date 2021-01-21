function PlotObject(obj,r,scale)
% Celestial object representation
%
% DESCRIPTION:
%   Given a date, this algorithm plots the position of whished planets
%   around the sun whith their respective orbits. When planets are
%   represented, their dimention with respect to their orbit is enlarged to
%   allow visibility.
%
% PROTOTYPE:
%   PlotObject(obj, r, scale)
%
% INPUT:
%   obj  [ 1 ]: Object selection. [ - ]
%          0 Sun
%          1 Mercury
%          2 Venus
%          3 Earth
%          4 Mars
%          5 Jupiter
%          6 Saturn
%          7 Uranus
%          8 Neptune
%   r     [ 3 x 1 ]: Position vector [ km ]
%   scale [ 1 ]: Factor used to enlarge planets
%
% OUTPUT:
%   Graphical representation of selected celestial object
%
% CALLED FUNCTIONS:
%
% AUTHORS:
%   Davide Demartini
%
% VERSIONS:
%   20-12-2020: First version
%
% FUTURE UPDATES:
%   Input check, every planet is a sphere and not a spheroid, insert a
%   factor that accounts how big the represented orbits are (to enlarge
%   planets plot).

% Input check (basic)
if nargin <= 2
    r = [0 0 0]; % Position the reference frame at the center of the planet
    scale = 1; % Display planet with real dimentions
end 

% Celestial object representation
if obj == 0 % Sun representation
    C = imread('TextureSun.jpg');
    theta = 0;
    radius = 696340*scale;
    [x, y, z] = ellipsoid(0, 0, 0, radius, radius, radius, 1E3);
    surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]),'FaceColor','texturemap','Edgecolor','none');
elseif obj == 1 % Mercury representation
    C = imread('TextureMercury.jpg');
    theta = 0;
    radius = 2439.7*scale;
    [x, y, z] = ellipsoid(r(1), r(2), r(3), radius, radius, radius, 1E2);
    surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]),'FaceColor','texturemap','Edgecolor','none');
elseif obj == 2 % Venus representation
    C = imread('TextureVenus.jpg');
    theta = 0;
    radius = 6051.8*scale;
    [x, y, z] = ellipsoid(r(1), r(2), r(3), radius, radius, radius, 1E2);
    surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]),'FaceColor','texturemap','Edgecolor','none');
elseif obj == 3 % Earth representation
    C=imread('TextureEarth.jpg');
    theta=0;
    radius = 6378.135*scale;
    [x, y, z] = ellipsoid(r(1), r(2), r(3), radius, radius, radius, 1E2);
    surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]),'FaceColor','texturemap','Edgecolor','none');
elseif obj == 4 % Mars representation
    C = imread('TextureMars.jpg');
    theta = 0;
    radius = 3389.5*scale;
    [x, y, z] = ellipsoid(r(1), r(2), r(3), radius, radius, radius, 1E2);
    surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]),'FaceColor','texturemap','Edgecolor','none');
elseif obj == 5 % Jupiter representation
    C = imread('TextureJupiter.jpg');
    theta = 0; 
    radius = 69911*scale;
    [x, y, z] = ellipsoid(r(1), r(2), r(3), radius, radius, radius, 1E2);
    surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]),'FaceColor','texturemap','Edgecolor','none');
elseif obj == 6 % Saturn representation
    C = imread('TextureSaturn.jpg');
    theta = 0; 
    radius = 58232*scale;
    [x, y, z] = ellipsoid(r(1), r(2), r(3), radius, radius, radius, 1E2);
    surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]),'FaceColor','texturemap','Edgecolor','none');
elseif obj == 7 % Uranus representation
    C = imread('TextureUranus.jpg');
    theta = 0; 
    radius = 25362*scale;
    [x, y, z] = ellipsoid(r(1), r(2), r(3), radius, radius, radius, 1E2);
    surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]),'FaceColor','texturemap','Edgecolor','none');
elseif obj == 8 % Neptune representation
    C = imread('TextureNeptune.jpg');
    theta = 0;
    radius = 24622*scale;
    [x, y, z] = ellipsoid(r(1), r(2), r(3), radius, radius, radius, 1E2);
    surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]),'FaceColor','texturemap','Edgecolor','none');
else
    error('Invalid CO selection');
end
axis equal;
end