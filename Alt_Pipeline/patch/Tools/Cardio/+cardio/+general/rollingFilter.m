function [filtered, varargout] = rollingFilter(signal, varargin)

% -------------------------------------------------------------------------
% This function filters a continuous signal using a filter applied to a
% rolling window. The filter criterion is set to the mean +/- the specified
% number of standard deviations. The removed indices may optionally be
% interpolated.
%
% Arguments (required)
% - signal      [Nx1]   Signal vector
%
% Arguments (optional)
% - tol                 Tolerance for outliers in stddevs (default: 2)
% - winRad              Radius of rolling window (default: 10)
% - overlap             Percent overlap [0, 1] of rolling window (default: 0.5)
% - no_interp	FLAG    Do not linearly interpolate outliers? (default: true)
% - verbose     FLAG    Display progress? (default: false)
%
% Outputs
% - varargout{1}    Indices of outliers
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'tol'); tol = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'winRad'); winRad = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'overlap'); overlap = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'interp'); interp = true;
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        end
    end
end

% Set defaults for optional inputs
if ~exist('tol', 'var'); tol = 2; end
if ~exist('winRad', 'var'); winRad = 100; end
if ~exist('overlap', 'var'); overlap = 0.5; end
if ~exist('interp', 'var'); interp = true; end
if ~exist('verbose', 'var'); verbose = false; end

% Initialize waitbar, if indicated
if verbose; h = waitbar(0, "Computing Rolling Window..."); end

% Obtain rolling window of signal vector (ends padded with nans)
[windowed, idx] = cardio.general.slideWindow(signal, winRad, overlap, 'nan');

% Get the number of windows
numWindows = size(windowed, 2);

% Set placeholder for outlier indices
outlierIdx = [];

% For each window...
for i = 1:numWindows
    
    % Update waitbar, if indicated
    if verbose; waitbar(i/numWindows, h, "Filtering Windows..."); end
    
    % Detrend the window
    windowed(:, i) = detrend(windowed(:, i));
    
    % Find the mean and standard deviation of the window
    winMean = nanmean(windowed(:, i)); winStd = nanstd(windowed(:, i));
    
    % Find outlier points
    outliers = find(windowed(:, i) < winMean - tol*winStd | windowed(:, i) > winMean + tol*winStd);
    
    % If there are outliers, perform the following
    if ~isempty(outliers)
    
        % Get the actual position of each outlier
        windowIdx = idx(i, 1):idx(i, 2);    % True indices for each point
        outliers = windowIdx(outliers);

        % Add outlier points to outlier index placeholder
        outlierIdx = [outlierIdx; outliers(:)];
    
    end
    
end; if verbose; close(h); end

% Limit outlier indices to unique values
outlierIdx = unique(outlierIdx);

% Set output placeholder
filtered = signal;

% If necessary, interpolate the outlier points
if interp; filtered(outlierIdx) = nan; ...
        filtered = fillmissing(filtered, 'linear'); end

% Format outputs
if nargout > 1; varargout{1} = outlierIdx; end

end

