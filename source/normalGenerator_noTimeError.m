function [data] = normalGenerator_noTimeError(theta,t,sigma,model,numSamples)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mu = model(theta,t);
% mu = log(muX.^2 ./ sqrt(muX.^2 + sigma.^2));
% sigma = sqrt(log(1 + sigma.^2 ./ (muX.^2)));

% replicate mu and sigma to get multiple samples at each timepoint
mu = repmat(mu,1,numSamples);
sigma = repmat(sigma,1,numSamples);
data = normrnd(mu,sigma); % rows are time points columns are samples

end