function indices = labelFeatures(signalSegments, featureType, numFeatures, ...
    maxSegments, signalEnsemble)

% -------------------------------------------------------------------------
% This function is for manually labeling features in time-series data given
% an array in which the columns are the signal segments. The user may
% specify the type of feature to extract, and the number of features per
% segment. The arguments are as follows:
% signalSegments    [M x N]     Array containing time-series column vectors
% featureType       {String}    Type of feature to extract (cell array)
%                           - "Peaks"       : Extract peaks
%                           - "Valleys"     : Extract valleys
%                           - "Inflections" : Extract inflections
% numFeatures       Number of features to extract per segment
% maxSegments       Maximum number of segments to analzye
%                           - "max"         : Analyze all entries
%                           - Int           : Set limit
% signalEnsemble    Ensemble averaged signal
% -------------------------------------------------------------------------

% Obtain the total number of signal segments
numSegments = size(signalSegments, 2);

% Modify by maximum number of segments indicated
if maxSegments ~= "max"
    if maxSegments < numSegments
        numSegments = maxSegments;
    end
end

% Initialize placeholder for return value
indices = zeros(numSegments, numFeatures);

% Determine requested feature types
peaks = false; valleys = false; % Set default flags
% If the feature type has been specified, update the flag
if ~isempty(find([featureType{:}] == "Peaks", 1)); peaks = true; end
if ~isempty(find([featureType{:}] == "Valleys", 1)); valleys = true; end
if ~isempty(find([featureType{:}] == "Inflections", 1)); inflections = true; end

% For each segment, plot the location of each feature and prompt the user
% to input which feature(s) is/are the desired feature(s)

% For each segment
for i = 1:numSegments
    
    % Get peak features for each segment
    [peakIdx, valleyIdx, infIdx] = ...
        cardio.general.getPeaks(signalSegments(:,i));
    
    % Set possible feature value placeholder
    possibleFeatures = [];
    
    % Populate possible feature vector
    if peaks; possibleFeatures = [possibleFeatures peakIdx]; end        % Peaks
    if valleys; possibleFeatures = [possibleFeatures valleyIdx]; end    % Valleys
    if inflections; possibleFeatures = [possibleFeatures infIdx]; end   % Inflections
    
    % Plot the segment with peaks labeled
    plot(signalSegments(:,i))
    hold on; grid on;
    plot(signalEnsemble)
    for j = 1:length(possibleFeatures)
        plot(possibleFeatures(j), signalSegments(possibleFeatures(j),i), 'ok')
    end
    xlabel("Sample"); ylabel("Amplitude");
    title("Possible features for Signal " + string(i))
    hold off

    % Accept user input
    [x, y] = ginput(numFeatures);

    % Find the possible feature points closest to each selected feature
    for j = 1:numFeatures
        distances = zeros(1, length(possibleFeatures));
        point1 = [x(j) y(j)];
        for k = 1:length(possibleFeatures)
            point2 = [possibleFeatures(k) signalSegments(possibleFeatures(k), i)];
            distances(k) = pdist2(point1, point2);
        end
        [~, minIdx] = min(distances);
        indices(i, j) = possibleFeatures(minIdx);
    end

end

% Close all figures
close all

end