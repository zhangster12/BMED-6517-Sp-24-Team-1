function features = rollingFeatures(signalSegments, numFeatures, M, tol, Plot, varargin)

% -------------------------------------------------------------------------
% SUMMARY
% This function extracts features from an SCG signal using by implementing
% processSCG.m as a rolling window. processSCG.m is run on a specified
% number of signal segments, and the optimal Gaussian Mixture Model (GMM)
% configuration is saved from the first run. For the second run, the window
% of SCG segments is shifted and the script is run again, using the saved
% GMM as a starting point for feature extraction. This process is repeated
% until all signal semgents have been analyzed.
%
% ARGUMENTS (REQ'D)
% - signalSegments  [NxM]   M vectors of SCG segments of length N
% - numFeatures             Number of time-domain features to extract
%                           - 1: first minimum
%                           - 2: first maximum
%                           - 3: second minimum ...
%                           - 6: third maximum
% - M                       Smoothing factor for exponential moving average
% - tol                     Tolerance for outliers (in standard deviations)
% - Plot            Bool    Plot results?
%
% ARGUMENTS (OPT'L)
% - K                       Maximum cluster size for mixed gaussian model
% - iterations              Number of iterations for mixed gaussian model
% - order                   Order number for lines of best fit
% - ao              [Nx1]   AO points for comparison
% - removed         {2x1}   Vectors of removed indices for AO point and signal
%                   {1}:	Removed indices for AO point
%                   {2}:    Removed indices for SCG signal
% - rawOnly         Bool    Raw indices only?
% - inflect         FLAG    Consider inflection points for resampling
% - windowLength            Length (in number of signal segments) of window
% - windowShift             Shift (in number of signal segments) of window
% - iterationDecay          Iteration drop after first pass
% - verbose         FLAG    Display progress?
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Part 0: Setup
% -------------------------------------------------------------------------

warning('off', 'all');

% Extract optional arguments (if necessary)
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'windowLength'); windowLength = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'windowShift'); windowShift = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'iterationDecay'); iterationDecay = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'iterations'); iterations = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        end
    end
end

% Set defaults
if ~exist('windowLength', 'var'); windowLength = size(signalSegments, 2); end
if ~exist('windowShift', 'var'); windowShift = 1; end
if ~exist('iterationDecay', 'var'); iterationDecay = 1; end
if ~exist('iterations', 'var'); iterations = 100; end
if ~exist('verbose', 'var'); verbose = false; end

% Set placeholder for return values
features = zeros(numFeatures, size(signalSegments, 2));

% -------------------------------------------------------------------------
% Part 1: Extract Features
% -------------------------------------------------------------------------

% Print progress, if indicated
if verbose; bar = waitbar(0, "Processing Initial Window"); end

% Extract the first window
window = signalSegments(:, 1:windowLength);

% Perform the first pass
SCG = cardio.scg.processSCG(window, numFeatures, M, tol, Plot, varargin{:});

% Write features to array
for f = 1:numFeatures; features(f, 1:windowLength) = SCG.features{f}(:); end

c = 1; % Initialize the counter

% Perform the following while there are still segments to analyze
while c*windowShift <= size(signalSegments, 2)
    
    % Print progress, if indicated
    if verbose; waitbar(c*windowShift/size(signalSegments, 2), bar, "Progress: " + string(c*windowShift) + ...
            " of " + string(size(signalSegments, 2))); end
    
    % Extract the new window
    endPoint = min(windowLength + c*windowShift, size(signalSegments, 2));
    window = signalSegments(:, c*windowShift:endPoint);
    
    % Update variable input arguments with GMM start point
    newVarargin = varargin;
    newVarargin{end + 1} = 'startModel'; newVarargin{end + 1} = SCG.startModel;
    newVarargin{end + 1} = 'iterations'; newVarargin{end + 1} = ceil(iterations/iterationDecay);
    
    % Process SCG
    SCG = cardio.scg.processSCG(window, numFeatures, M, tol, Plot, newVarargin{:});
    
    % Write features to array
    for f = 1:numFeatures
        features(f, c*windowShift:endPoint) = SCG.features{f}(:);
    end
    
    c = c + 1;  % Increment the counter
    
end

% Close waitbar if necessary
if verbose; close(bar); end

end