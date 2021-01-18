function [lon,lat] = RAdelta2latlon(alpha,delta,gw_lon0,omega_r,t,t0)
% r2RAdelta function that converts the position vector into the
%   rightascention and delta  
%   
% PROTOTYPE:
%   [lon, lat] = RAdelta2latlon(alpha, delta, gw_lon0, omega, t, t0)
% 
% INPUT:
%   alpha   [ n x 1 ]    Right ascension vector            [ rad ]
%   delta   [ n x 1 ]    Declination vector                [ rad ]
%   gw_lon0 [ 1 ]        Greenwhich longitude at t0        [ rad ]
%   omega_r [ 1 ]        Angular velocity of the Earth     [ rad/s ]
%   t       [ n ]        Vector of times at which the ground track will be 
%                         computed                         [ s ]
%   t0      [ 1 ]        Initial time                      [ s ]
%
% OUTPUT:
%   gw_lon  [ n x 3 ]    Greenwhich longitude at t         [ rad ]
%   lon     [ n x 1 ]    Object longitude at t             [ rad ]
%   lat     [ n x 1 ]    Object latitude at t              [ rad ]
%
% CONTRIBUTORS:
%   Davide Iafrate      
%   Alkady Marwan
%   Pedro Bossi Núñez
%   Davide Demartini
%
% VERSIONS
%   15-12-2020: First version
%   16-01-2021: Second version

gw_lon = gw_lon0 + omega_r*(t - t0);
lon = wrapTo2Pi(alpha - gw_lon);
lat = delta;
end