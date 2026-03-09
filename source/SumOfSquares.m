function [SSE] = SumOfSquares(fitValues,realValues)
%CALCULATERESIDUALS Summary of this function goes here
% Columns of expData need to be repeats and the rows the time points
% fitValues need to be a column vector with the elements corresponding to a
% value at each of the timepoints of expData
% The adipocyte files are set up like this, just need to ensure fit values
% are column vector

fitValues = fitValues(:); % make a column vector

SSE = (fitValues - realValues);
SSE = sum(SSE.^2,'all','omitmissing');