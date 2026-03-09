function [orderedTime, randomOrderTime] = generate_correlatedTimes(t,muX,sigma,numSamples,distribution)
%UNTITLED Summary of this function goes here
%   Measurements at each timepoint is correlated
% Each measurement has some lognormally distributed error that occurs after
% the previous measurement
% One of the measurements is always taken at the correct time,
% this starts the clock on how long it takes to collect the rest of the
% measurements
% Using lognormal distribution as time between successive measurements is 
%
% orderedTime output has the correlated times always in the same order,
% i.e., every day the measurements are taken on mouse 1, then mouse 2, etc.
% randomOrderTime output has the correlated times in a random order, i.e.,
% experimentalist randomly selects an order of mcie to measure each day

t = t(:)'; % make sure time is row vector
t = repmat(t,numSamples,1);

% Get correlated errors
if strcmp(distribution,'lognormal')
    % lognormal parameters
    mu = log(muX.^2 ./ sqrt(muX.^2 + sigma.^2));
    sigma = sqrt(log(1 + sigma.^2 ./ (muX.^2)));
    errors = lognrnd(mu,sigma,size(t));
elseif strcmp(distribution,'normal')
    errors = normrnd(muX,sigma,size(t));
end
errors = cumsum(errors,1);

% Apply errors to ordered times
orderedTime = t + errors;

% Apply to random ordering
randomOrderError = errors;
for i = 1:size(t,2)
    randomOrderError(:,i) = randomOrderError(randperm(numSamples),i);
end
randomOrderTime = t + randomOrderError;

randomOrderTime = randomOrderTime';
orderedTime = orderedTime';

end