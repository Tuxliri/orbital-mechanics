function DV = designer_opt(planetA,planetB,windowA,windowsB)
%DESIGNER_OPTIMIZATION This is a wrapper to find the minimum DV using the
%fminunc or fmincon optimization problems
% INPUT:
%   planetA[1]      departure planet ID (1-8)
%   planetB[1]      arrival planet ID (1-8)
%   x[2x2]          vector with dates in the following 

[DV, ~, ~, ~, ~,~, ~,~,~]...
    = tfdesigner(planetA,planetB,windowA,windowsB,1,1);
end

