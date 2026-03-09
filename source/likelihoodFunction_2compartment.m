function [likelihood] = likelihoodFunction_2compartment(expData,experiments,insLevels,kex,ken,Total)
%Calculates the log likelihood function given a set of parameters and
%experimental data for the 2 compartment model
%   Given a set of parameters, calculates the likelihood function

[output, ~, ~]...
    = compartment2Combined_totalInsulin(insLevels,kex,ken,Total);


assignin('base','output',output);

% In first instance, assume sigma = 1. sigma is std, sigma^2 = variance
variance = 1;

likelihood = 1;

for i = 1:length(experiments)
    % Loop through all insulin levels for this experiment
    % Get second level struct fieldnames, these should corresond to insulin
    % level
    currExp = experiments{i};
    insLevel = fieldnames(expData.(currExp));
    for j = 1:length(insLevel)
        currIns = insLevel{j};
        if ~isempty(expData.(currExp).(currIns))
            % Calculate the SSE
            modelData = output(strcmp(output.Expt,currExp),:);
            modelData = modelData(strcmp(modelData.InsulinTag,currIns),:);

            % need to interpolate model output to matching timepoints of data
            currTimes = expData.times.(currExp);
            y=interp1(modelData.Time,modelData.Data,currTimes);

            % loop each timepoint
            % for i = 1:length(currTimes)
            %     errors = (y(i) - expData.(currExp).(currIns)(i,:)); % matrix where each row is errors at a timepoint
            %     likelihood = likelihood - sum(errors.^2 ./ (2 * variance) + 1/2 * log(2*pi*variance),'all','omitmissing'); % variance is a column vector consisting of the variance at each time point
            % end
            errors = y(:) - expData.(currExp).(currIns); % matrix where each row is errors at a timepoint
            % variance = var(errors,0,2,'omitmissing') % variance of errors at each timepoint
            likelihood = likelihood * prod(1 ./ sqrt(2 * pi * variance) .*...
                exp(-(errors.^2) ./ (2*variance)),'all','omitmissing'); % variance is a column vector consisting of the variance at each time point
        end
    end
end

likelihood = -likelihood;
end