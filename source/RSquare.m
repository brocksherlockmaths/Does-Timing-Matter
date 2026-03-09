function [rsquare, SST, ybar] = RSquare(SSE,expData)
%RSQUARE Summary of this function goes here
%   Takes the

ybar = mean(expData,'all');
SST = sum((expData - ybar).^2,'all');

rsquare = 1 - SSE / SST;
end

