function [SST] = sumOfSquareTotal(expData)
%SUMOFSQUARETOTAL Summary of this function goes here
%   Detailed explanation goes here

ybar = mean(expData,'all');
SST = sum((expData - ybar).^2,'all');

end

