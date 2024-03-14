function label = labelSegmentFeatures(signalSegment, featureType, numFeatures, ...
    segmentNumber, signalEnsemble)

% -------------------------------------------------------------------------
% This function is for manually labeling features in time-series data given
% an individual signal segment as a vector rather than array. The user may
% specify the type of feature to extract, and the number of features per
% segment. The arguments are as follows:
% signalSegment     [N x 1]     Vector containing signal segment
% featureType       {String}    Type of feature to extract (cell array)
%                           - "Peaks"       : Extract peaks
%                           - "Valleys"     : Extract valleys
%                           - "Inflections" : Extract inflections
% numFeatures       Number of features to extract per segment
% segmentNumber     Identifier for current segment
% signalEnsemble    Ensemble averaged signal
% -------------------------------------------------------------------------

% Initialize placeholder for return value
label = zeros(1, numFeatures);

% Determine requested feature types
peaks = false; valleys = false; inflections = false; % Set default flags
% If the feature type has been specified, update the flag
if ~isempty(find([featureType{:}] == "Peaks", 1)); peaks = true; end
if ~isempty(find([featureType{:}] == "Valleys", 1)); valleys = true; end
if ~isempty(find([featureType{:}] == "Inflections", 1)); inflections = true; end

% For each segment, plot the location of each feature and prompt the user
% to input which feature(s) is/are the desired feature(s)

% Get peak features for each segment
[peakIdx, valleyIdx, infIdx] = ...
    cardio.general.getPeaks(signalSegment);

% Set possible feature value placeholder
possibleFeatures = [];

% Populate possible feature vector
if peaks; possibleFeatures = [possibleFeatures; peakIdx]; end        % Peaks
if valleys; possibleFeatures = [possibleFeatures; valleyIdx]; end    % Valleys
if inflections; possibleFeatures = [possibleFeatures; infIdx]; end   % Inflections

% Plot the segment with peaks labeled
plot(signalSegment)
hold on; grid on;
plot(signalEnsemble)
for j = 1:length(possibleFeatures)
    plot(possibleFeatures(j), signalSegment(possibleFeatures(j)), 'ok')
end
xlabel("Sample"); ylabel("Amplitude");
title("Possible features for Signal " + string(segmentNumber))
hold off

% Accept user input
[x, y] = ginput(numFeatures);

% Find the possible feature points closest to each selected feature
for j = 1:numFeatures
    distances = zeros(1, length(possibleFeatures));
    point1 = [x(j) y(j)];
    for k = 1:length(possibleFeatures)
        point2 = [possibleFeatures(k) signalSegment(possibleFeatures(k))];
        distances(k) = pdist2(point1, point2);
    end
    [~, minIdx] = min(distances);
    label(j) = possibleFeatures(minIdx);
end

% Close all figures
close all

end

