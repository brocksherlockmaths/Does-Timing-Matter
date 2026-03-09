function [squareDistance] = orthonormalSquareError(tdat,ydat,model,p)
%orthonormalSquareError Calculates the sum square distances from data
%points to the closest point on the model
%   model is a function handle, first argument is time (independent
%   variable), and second argument is the model parameters as a vector
%   We want to use this function in another minimisation to find the
%   parameters that minimise the distance
%   tdat and ydat must be the same size so that they are coordinate pairs

% Setup square distance
squareDistance = 0;

% Loop through each data point
for i = 1:numel(ydat)
    % Ned to minimise the square distance between point and function
    % i.e., minimise D^2 = (tdat - t)^2 + (ydat - y)^2

    % Function of the distance to minimise
    distFun = @(t) (tdat(i) - t).^2 + (ydat(i) - model(t,p)).^2;
    % Minimise this function with respect to t
    % Minimum should be close to the tdat for the data, so use as
    % startpoint
    [~, minDist] = fminsearch(distFun,tdat(i));
    % Add minDist (the shortest distance) to the square distance
    squareDistance = squareDistance + minDist;

end