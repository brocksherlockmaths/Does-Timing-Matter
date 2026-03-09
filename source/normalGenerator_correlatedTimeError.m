function [data] = normalGenerator_correlatedTimeError(theta,t,sigma,model,epsilon,numSamples)
%UNTITLED Summary of this function goes here
%   Measurements at each timepoint is correlated
% Each measurement has some lognormally distributed error that occurs after
% the previous measurement


t = repmat(t,1,numSamples);
timeError = normrnd(0,epsilon,size(t));

t = t + timeError;

mu = model(theta,t);


% replicate mu and sigma to get multiple samples at each timepoint
mu = repmat(mu,1,numSamples);
sigma = repmat(sigma,1,numSamples);
data = normrnd(mu,sigma); % rows are time points columns are samples

end