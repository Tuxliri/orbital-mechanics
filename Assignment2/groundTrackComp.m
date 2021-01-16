function [alpha,delta,lon,lat] = groundTrackComp...
    (r_vect,v_vect,theta_g_0,timevect,mu,om_e)
%
%% ORBITAL MECHANICS - groundTrackComp.m
%
% PROTOTYPE:
% [alpha,delta,lon,lat] = groundTrackComp...
%    (r_vect,v_vect,theta_g_0,timevect,mu,om_e)
%
% DESCRIPTION:
%  This program computes the groundtrack of the satellite
%
% INPUT:
%	r_vect       position vector in geocentric coordinates  [Km]
%   v_vect       velocity vector in geocentric coordinates  [Km/s]
%   theta_g_0    longitude of Greenwich at initial time     [deg]
%   timevect     vector of times for ground track computing [s]
%   mu           planetary constant                         [Km^3/s^2]
%   om_e         Earth's angular velocity                   [rad/s]
%
% OUTPUT:
%   alpha        right ascension in ECEI frame              [deg]
%   delta        declination in ECEI frame                  [deg]
%   lon          longitude wrt rotating Earth               [deg]
%   lat          latitude wrt rotating Earth                [deg]
%
% -------------------------------------
%  Copyright (C) 2020 Pedro Bossi Nunez  <pedro.bossi@pm.me>
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% CODE:
%
% set options
options=odeset('RelTol',1e-13,'Abstol',1e-14);
%
% Variable initialisation
tlen=length(timevect);
alpha=ones(tlen,1);
delta=ones(tlen,1);
lon=ones(tlen,1);
theta_g=ones(tlen,1);
y0=[r_vect v_vect]';
%
% Orbit propagation
[t,Y]= ode113(@(t,y) ode_2body(t, y, mu),timevect,y0,options);
t0=t(1);
pos=[Y(:,1) Y(:,2) Y(:,3)];
pos=[pos(:,1) pos(:,2) pos(:,3)];
vel=[Y(:,4) Y(:,5) Y(:,6)];
%
% Conversion to RA and declination
for i=(1:tlen)
    x=pos(i,1);
    y=pos(i,2);
    z=pos(i,3);
    r=norm(pos(i,:),2);
    delta(i)=asin(z/r);
    alpha(i)=acos((x/r)/cos(delta(i)));
    if sign(y/r)>0
        alpha(i)=alpha(i);
    else
        alpha(i)=2*pi-alpha(i);
    end
end
%
% Conversion to longitude and latitude
for j=(1:tlen)
    theta_g(j)=theta_g_0+om_e*(t(j)-t0);
    theta_g(j)=mod(theta_g(j),2*pi);
    lon(j)=alpha(j)-theta_g(j);
    if lon(j)>pi
        lon(j)=lon(j)-2*pi;
    end
    if lon(j)<-pi
        lon(j)=lon(j)+2*pi;
    end
end
lat=delta.*180./pi;
lon=lon.*180./pi;
end