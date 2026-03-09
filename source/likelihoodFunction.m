function [likelihood] = likelihoodFunction(tdat,ydat,sigma,model,p)
%Calculates the log likelihood function given a set of parameters and
%experimental data for the 2 compartment model
%   tdat is a column vector with the times of data measurements
%   ydat is a matrix where each row corresponds to measurements taken at
%   each time given in tdat
%   model is a function that takes parameters and computes the model output
%   p is the parameters to be fit

y = model(tdat,p); % calculate the model at each time point

% Use sigma = -1 to specify calculating sigma from the data
if sigma == -1
    errors = y - ydat;
    sigma = std(errors,0,2);
    sigma = repmat(sigma,1,size(ydat,2));
else
    sigma = repmat(sigma,size(ydat,1),size(ydat,2));
end


% Need y and sigma to be same dimension as ydat
y = repmat(y,1,size(ydat,2));


%log(normpdf(ydat,y,sigma))
%logLikelihood = -sum(log(normpdf(ydat,y,sigma)), 'all');
coeff = prod(1./sqrt(2*pi.*sigma.^2),"all");
mu = sum((ydat - y).^2 ./ (2*sigma.^2), 'all');

likelihood = coeff * exp(-mu);

% logLikelihood = 0;
% 
% 
% errors = y(:) - ydat; % matrix where each row is errors at a timepoint
% %variance = var(errors,0,2,'omitmissing'); % variance of errors at each timepoint
% logLikelihood = logLikelihood + sum(errors.^2 ./ (2 * variance) + 1/2 * log(2*pi*variance),'all','omitmissing'); % variance is a column vector consisting of the variance at each time point
% 
% 
% logLikelihood = -logLikelihood;
end