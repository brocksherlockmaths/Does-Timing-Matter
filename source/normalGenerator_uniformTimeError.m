function [data] = normalGenerator_uniformTimeError(theta,t,sigma,model,interval,numSamples)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


t = repmat(t,1,numSamples);
timeError = unifrnd(interval(1),interval(2),size(t));

t = t + timeError;

mu = model(theta,t);


% replicate mu and sigma to get multiple samples at each timepoint
mu = repmat(mu,1,numSamples);
sigma = repmat(sigma,1,numSamples);
data = normrnd(mu,sigma); % rows are time points columns are samples

end