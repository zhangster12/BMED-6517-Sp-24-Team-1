function features = detrendFeatures(featureArray)

% -------------------------------------------------------------------------
% This function detrends feature vectors organized in a cell array,
% returning a cell array of detrended features. The input arguments are:
% - featureArray    {N x 1} Cell array of feature vectors
% -------------------------------------------------------------------------

% Set return value placeholder
features = cell(size(featureArray));

% For each features
for i = 1:length(featureArray)
    
    % Extract feature vector
    temp = featureArray{i};
    
    % Detrend vector
    temp = detrend(temp(~isnan(temp)));
    
    % Initialize counter
    counter = 1;
    
    % Populate return value placeholder
    % For each feature value...
    for j = 1:length(featureArray{i})
        
        % If the value is NaN, write NaN to the return value
        if isnan(featureArray{i}(j))
            features{i}(j) = NaN;
        else
            % Else, write the detrended value
            features{i}(j) = temp(counter);
            % Increment counter
            counter = counter + 1;
        end
        
    end
    
end

end

