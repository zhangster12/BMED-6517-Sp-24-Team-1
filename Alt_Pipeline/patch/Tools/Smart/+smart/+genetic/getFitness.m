function [fitness, normalized, ranks] = getFitness(inputArray, P, F)

% -------------------------------------------------------------------------
% This function calculates the fitness of the organism represented by the
% input array, defined as the average standard deviation across all vectors
% of the input array. Input parameters are as follows:
% - inputArray  {NxM}   Columnwise array of cells forming M vectors
% - P                   Propensity to select features with higher fitness
% - F                   Method of calculating fitness
%                       - 'avgstd': Mean standard deviation
%                       - 'maxstd': Maximum standard deviation
% -------------------------------------------------------------------------

% Initialize placeholder for return value
fitness = zeros(length(inputArray), 1);

for i = 1:length(inputArray)
    
    % Vectorize input array
    % Arguments: inputArray, Plot
    vectors = smart.general.plotCell(inputArray{i}, false);
    
    % Handle NaNs by filling missing values
    vectors = fillmissing(vectors, 'nearest');
    
    % Detrend vectors
    dVectors = detrend(vectors, 'linear');
   
    % Get fitness scores
    if strcmp(F, 'avgstd')
        fitness(i) = 1/mean(std(dVectors));
    elseif strcmp(F, 'maxstd')
        fitness(i) = 1/max(std(dVectors));
    elseif strcmp(F, 'scgmean')
        % This is a custom fitness score that may be edited by the user
        fitness(i) = 1/mean(std(dVectors(:,1:4)));
    elseif strcmp(F, 'scgmax')
        % This is a custom fitness score that may be edited by the user
        fitness(i) = 1/max(std(dVectors(:,1:4)));
    end
    
end

% Adjust for propensity
fitness = fitness.^P;

% Normalize fitness
normalized = fitness/norm(fitness, 1);

% Obtain rank
[~, ranks] = sort(normalized, 'descend');

end

