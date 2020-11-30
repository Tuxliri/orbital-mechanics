function [c,ceq] = launcherCONSTR(windowA,windowB,planetA,planetB,vinf)
%LAUNCHERCONSTR function that builds the matrices for the constrained
%minmum optimization search
% The constraint is that DV1 < vinf
% INPUT:
%   x[] vector containing the departure and arrival windows in MJD2000
%       format
%   planetA[1]      departure planet ID (1-8)
%   planetB[1]      arrival planet ID (1-8)
%   vinf[1] maximum excess velocity from launcher [km/s]
% OUTPUT:
%   c[1] DV1 - vinf [km/s] nonlinear inequality matrix
%   ceq[]   nonlinear equality matrix
[~, ~, ~, DV1, ~,~, ~,~,~]...
    = tfdesigner(planetA,planetB,windowA,windowB,1,1);
c = DV1-vinf;
ceq = [];
end

