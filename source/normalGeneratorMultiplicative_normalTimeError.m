function [data] = normalGeneratorMultiplicative_normalTimeError(theta,t,sigma,model,epsilon,numSamples)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


t = repmat(t,1,numSamples);
timeError = normrnd(0,epsilon,size(t));

t = t + timeError;

data = model(theta,t);


% replicate mu and sigma to get multiple samples at each timepoint
mu = zeros(1,numSamples);
sigma = repmat(sigma,1,numSamples);
data = data + data.*normrnd(mu,sigma); % rows are time points columns are samples

end