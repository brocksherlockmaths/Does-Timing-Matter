function [logLikelihood] = linearLogLikelihood(tdat,ydat,theta)
%linearLogLikelihood computes the log likelihood function for a linear
%model y(t) = a + b*t where a and b are normally distributed parameters
%   Detailed explanation goes here
% theta(1) = mu_a
% theta(2) = sigma_a
% theta(3) = mu_b
% theta(4) = sigma_b

mu_a = theta(1);
sigma_a = theta(2);
mu_b = theta(3);
sigma_b = theta(4);

% Make data and time paramters column vectors
ydat = ydat(:);
tdat = tdat(:);

% Compute the log likelihood
% logLikelihood = sum(-log(sigma_a) - log(sqrt(2*pi)) - (ydat - mu_a).^2 ./ 2 ./ (sigma_a^2) + ...
%     log(1 + sigma_a / sigma_b .* exp((ydat-mu_a).^2 ./ 2 ./ (sigma_a^2) - ...
%     (ydat-mu_b).^2 .* tdat ./ 2 ./ (sigma_b^2))));

% logLikelihood = sum(log(1/(sigma_a*sqrt(2*pi)) * exp(-(ydat-mu_a).^2 ./2 ./ (sigma_a^2)) + ...
%     1/(sigma_b*sqrt(2*pi)) * exp(-(ydat-mu_b).^2 ./2 ./ (sigma_b^2))));

logLikelihood = sum(-(ydat-mu_a - mu_b.*tdat).^2 ./ (sigma_a^2 + sigma_b^2*tdat.^2)./2 - log(sqrt(sigma_a^2 + sigma_b^2*tdat.^2)) - log(sqrt(2*pi)),'all');

end