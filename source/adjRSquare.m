function [adjrsquare] = adjRSquare(SSE, SST, p, n)
%ADJRSQUARE Summary of this function goes here
%   R2 is the calculated r square value,
%   n is the number of observations
%   p is the number of regression coefficients

adjrsquare = 1 - ((n - 1) / (n - p)) * SSE/SST;
end

