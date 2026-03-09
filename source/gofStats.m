function [gof] = gofStats(expData,SSE,p,n)
%GOFSTATS Summary of this function goes here
%   n is number of points
%   p is the number of parameters

[rsquare, SST] = RSquare(SSE,expData);
[adjrsquare] = adjRSquare(SSE, SST, p, n);
rmse = rootMeanSquareError(SSE,p,n);


gof.SSE = SSE;
gof.SST = SST;
gof.rsquare = rsquare;
gof.adjrsquare = adjrsquare;
gof.rmse = rmse;
end

