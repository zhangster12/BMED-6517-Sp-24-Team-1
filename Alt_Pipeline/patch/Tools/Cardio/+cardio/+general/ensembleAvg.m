function averaged = ensembleAvg(signals, varargin)

% -------------------------------------------------------------------------
% This function performs a rolling-window ensemble average.
%
% Arguments (req'd)
% - signal      [MxN]       N signals of length M
%
% Arguments (opt'l)
% - 'windowLength'  Int     Window length for rolling window
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'windowLength'); windowLength = varargin{arg + 1};
        end
    end
end

% Set defaults
if ~exist('windowLength', 'var'); windowLength = 5; end

% Get the number of signals
numSignals = size(signals, 2);

% Determine upper and lower window bound
bound = floor((windowLength-1)/2);

% Set placeholder for return value
averaged = signals;

% Each signal slice should be replaced with the proper ensemble average
for i = 1:numSignals
    
    % Determine the upper and lower bounds
    lowerBound = max(1, i - bound); upperBound = min(numSignals, i + bound);
    
    % Compute the ensemble average for slices within the bounds
    averaged(:, i) = mean(signals(:, lowerBound:upperBound), 2);
    
end