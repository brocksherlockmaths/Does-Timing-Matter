function [RMSE] = rootMeanSquareError(SSE, p, n)
%ROOTMEANSQUAREERROR Summary of this function goes here
%   n is number of points, p is number of fitted coefficients

RMSE = sqrt(SSE / (n-p));
end

