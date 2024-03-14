function outputArray = plotCell(inputArray, Plot)

% -------------------------------------------------------------------------
% This function accepts an input cell array with columns of cells
% representing a segmented time series, reconstructing the time series and
% optionally plotting the results. Input arguments are as follows:
% - inputArray  {NxM}   Input cell array with columnwise data vectors
% - Plot        Bool    Plot results
% -------------------------------------------------------------------------

% Initialize return value placeholder
outputArray = [];

% For each column...
for i = 1:size(inputArray, 2)
    
    % Obtain the resulting data vectors
    temp = cat(1, inputArray{:,i});
    
    % Append the results to the output array
    outputArray = [outputArray temp];
    
end

% Plot results
if Plot
    figure; hold on; grid on;   % Initialize figure
    for i = 1:size(inputArray, 2)
        area(outputArray(:,i), 'FaceAlpha', 0.3)
    end
end

end

