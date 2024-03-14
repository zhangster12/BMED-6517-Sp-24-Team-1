function indices = getOutliers(featureArray, statsArray, tol, numFeatures)

% -------------------------------------------------------------------------
% This function returns the outlier indices in for each feature vector in a
% feature cell array, comparing the feature value to a provided vector of
% expected feature values and corresponding standard deviation. The input
% arguments are as follows:
% - featureArray:   {N x 1} Cell array containing feature vectors
% - statsArray:     {M x 2} Cell array containing expected feature values
%                           as vectors in column 1 and the standard 
%                           deivation in column 2
% - tol:                    Tolerance for outliers (in standard deviations)
% - numFeatures             Number of features to analyze
% -------------------------------------------------------------------------

% Ensure numFeatures des not surpass the number of features
if numFeatures > length(featureArray); numFeatures = length(featureArray); end

% Set placeholder for return value
indices = cell(numFeatures, 1);

% Perform analysis for each feature vector
for i = 1:numFeatures
    
    % Extract vector
    temp = featureArray{i};
    
    % Initialize counter
    counter = 1;
    
    % Remove outliers
    for j = 1:length(temp)
        
        % Don't perform the following if the value is NaN
        if isnan(temp(j)); continue; end
        
        % Calculate bounds
        lowerBound = statsArray{i, 1}(j) - tol*statsArray{i, 2};
        upperBound = statsArray{i, 1}(j) + tol*statsArray{i, 2};
        
        % Return the index if the value is an outlier
        if temp(j) < lowerBound || temp(j) > upperBound
            indices{i}(counter) = j; counter = counter + 1;
        end
        
    end
    
end

end

