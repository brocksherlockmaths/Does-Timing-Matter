function [data] = logNGenerator_biasedUniform(theta,t,sigma,model,epsilon,numSamples)
%UNTITLED Summary of this function goes here
%   epsilon is the upper end of intervalfor uniform biased noise

% Each measurement has a uniwue time biased error
% Need t to be replicated for each number of samples to generate
% independent noise on each measurement
t = repmat(t,1,numSamples);
% add the noise
if epsilon <= 0 % biased to early measurements
    timeError = unifrnd(epsilon,0,size(t));
else % biased late measurments
    timeError = unifrnd(0,epsilon, size(t));
end
t = t + timeError;

muX = model(theta,t);
mu = log(muX.^2 ./ sqrt(muX.^2 + sigma.^2));
sigma = sqrt(log(1 + sigma.^2 ./ (muX.^2)));

% replicate mu and sigma to get multiple samples at each timepoint
%mu = repmat(mu,1,numSamples);
%sigma = repmat(sigma,1,numSamples);
data = lognrnd(mu,sigma); % rows are time points columns are samples

end