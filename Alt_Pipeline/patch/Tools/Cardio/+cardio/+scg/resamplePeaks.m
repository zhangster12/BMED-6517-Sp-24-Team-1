function [resampledIndex] = resamplePeaks(originalIndex, peakArray,...
    infArray, expectedValue)

% -------------------------------------------------------------------------
% This function resamples peaks based on statistics about where the
% desired peak is most likely located. The function checks whether the
% extracted peak is within X standard deviations of expected; if not, it
% checks to see whether any points are closer to the mean; if so, it
% returns the point closest to the mean.
%
% The function accepts the following arguments:
% - originalIndex:  Index of the original selected peak
% - peakArray:      Vector containing peak indices
% - infArray:       Vector containing inflection point indices
% - expectedValue:  Expected value of the distribution for the desired peak
% -------------------------------------------------------------------------

% Define boolean flag indicating whether peak has been updated
updateFLAG = false;

% Get distance to the mean
distance2Mean = abs(originalIndex - expectedValue);

% Are there any points closer to the mean?
% For each peak, get the distance to the mean.
for j = 1:length(peakArray)
    % If the distance to the mean is less than distance2Mean,
    % the wrong point was selected. Reselect point.
    if abs(expectedValue - peakArray(j)) < distance2Mean
        distance2Mean = abs(expectedValue - peakArray(j));  % Reset distance2Median
        originalIndex = peakArray(j);
        updateFLAG = true;  % Set updateFLAG
    end
end

% If, after looping through all the peaks, the updateFLAG is still
% false, check to see whether there are inflection points closer to the
% mean, indicating that the true peak is possibly obscured by noise.
if ~updateFLAG
    % For each inflection point, get the distance to the mean.
    for j = 1:length(infArray)
        % If the distance to the mean is less than distance2Median, the
        % wrong point was selected. Reselect point.
        if abs(expectedValue - infArray(j)) < distance2Mean
            distance2Mean = abs(expectedValue - infArray(j));    % Reset distance2Median
            originalIndex = infArray(j);
        end
    end
end

% Return resampled index
resampledIndex = originalIndex;

end