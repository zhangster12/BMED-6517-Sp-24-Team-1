function output = normalizeByLevel(obj, feature, varargin)

% -------------------------------------------------------------------------
% This function normalizes a feature vector with respect to its mean value
% in a specified level.
%
% Arguments (required)
% - feature     Feature     Feature to normalize
%
% Arguments (optional)
% - subject                 Subject for which to compute output
% - baseline    Level       Baseline level for which to compute mean value
%                           (default Level.absBasleine1)
% - return      Level       Level to compare against baseline
%                           (default Level.allAbsolute)
% - ratio       FLAG        Return the output as a ratio (default % change)
% - verbose     FLAG        Print output to console?
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        switch varargin{arg}
            case 'subject'; subject = varargin{arg + 1};
            case 'baseline'; baselineLevel = varargin{arg + 1};
            case 'return'; returnLevel = varargin{arg + 1};
            case 'ratio'; ratio = true;
            case 'verbose'; verbose = true;
        end
    end
end

% Set defaults for optional input arguments
if ~exist('subject', 'var'); subject = 1; end
if ~exist('baselineLevel', 'var'); baselineLevel = Level.absBaseline1; end
if ~exist('returnLevel', 'var'); returnLevel = Level.allAbsolute; end
if ~exist('ratio', 'var'); ratio = false; end
if ~exist('verbose', 'var'); verbose = true; end

% Print progress, if applicable
if verbose; disp("-> Extracting data for baseline level..."); end

% Extract the feature vector at the baseline level
if verbose
    [baselineData, ~] = obj.trainingData('subjects', subject, 'levels', baselineLevel, 'features', feature, 'verbose');
else
    [baselineData, ~] = obj.trainingData('subjects', subject, 'levels', baselineLevel, 'features', feature);
end

% Compute the mean value for the baseline data
baselineMean = mean(baselineData);

% Print progress, if applicable
if verbose; disp("-> Extracting data for return level..."); end

% Extract the feature vector for the return level
if verbose
    [returnData, ~] = obj.trainingData('subjects', subject, 'levels', returnLevel, 'features', feature, 'verbose');
else
    [returnData, ~] = obj.trainingData('subjects', subject, 'levels', returnLevel, 'features', feature);
end

% Normalize the return data by the baseline data
if ratio
    output = returnData./baselineMean;
else
    output = (returnData - baselineMean)./baselineMean;
end
    
end

