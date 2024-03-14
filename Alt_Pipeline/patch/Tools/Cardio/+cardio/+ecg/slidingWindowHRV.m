function [hrv, varargout] = slidingWindowHRV(indices, Fs, varargin)

% -------------------------------------------------------------------------
% This function computes the heart rate variability over an extended
% recording interval using a sliding window. This function computes heart 
% rate variability using one of three methods:
% 1. R-R interval RMS difference ('difference')
% 2. Power Spectral Density - High vs. Low Frequency ('spectral')
% 3. Poincare Method ('poincare')
%
% Arguments (required)
% - indices     [Nx1]   R-peak time indices (ms)
% - Fs                  Sampling frequency (Hz)
%
% Arguments (optional)
% - difference  FLAG    Use window method for HRV (default)
% - spectral    FLAG    Use spectral method for HRV
% - poincare    FLAG    Use Poincare method for HRV
% - winLength           Window length (samples)
% - overlap             Percent overlap [0, 1]
% - verbose     FLAG    Display waitbar with progress indicator
% - ends        FLAG    Preserve ends of signal (with nan padding)
% - tol                 Tolerance for outliers (in stardard deviations)
%
% Outputs
% - hrv         [Nx1]   Heart rate variability
% - LF/Cov1     [Nx1]   Low frequency component (spectral) or first
%                       eigenvalue (poincare)
% - HF/Cov2     [Nx1]   High frequency component (spectral) or second
%                       eigenvalue (poincare)
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'Fs'); Fs = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'difference'); difference = true;
        elseif strcmp(varargin{arg}, 'spectral'); spectral = true;
        elseif strcmp(varargin{arg}, 'poincare'); poincare = true;
        elseif strcmp(varargin{arg}, 'winLength'); winLength = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'overlap'); overlap = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        elseif strcmp(varargin{arg}, 'ends'); ends = true;
        elseif strcmp(varargin{arg}, 'tol'); tol = varargin{arg + 1};
        end
    end
end

% Set defaults for optional input arguments
if ~exist('winLength', 'var'); winLength = round(length(indices)/10); end
if ~exist('overlap', 'var'); overlap = 0.5; end
if ~exist('difference', 'var'); difference = false; end
if ~exist('spectral', 'var'); spectral = false; end
if ~exist('poincare', 'var'); poincare = false; end
if ~exist('verbose', 'var'); verbose = false; end
if ~difference && ~spectral && ~poincare; difference = true; end
if ~exist('ends', 'var'); ends = false; end
if ~exist('tol', 'var'); tol = inf; end

% If the spectral method is chosen but Fs or ECG is not provided, return an
% error message to the console
if spectral && ~exist('Fs', 'var')
    disp("-> Error in slidingWindowHRV(): Fs must be provided when using spectral method")
    hrv = nan; return
end

% Convert the window length into a window radius
winLength = round(winLength/2);

% Create a sliding window of the indices
if ~ends
    windowedIndices = cardio.general.slideWindow(indices, winLength, overlap);
else
    windowedIndices = cardio.general.slideWindow(indices, winLength, overlap, 'preserveEnds', 'nan');
end

% Create a placeholder for return value
hrv = zeros(size(windowedIndices, 2), 1);
HF = zeros(size(windowedIndices, 2), 1);
LF = zeros(size(windowedIndices, 2), 1);
Cov1 = zeros(size(windowedIndices, 2), 1);
Cov2 = zeros(size(windowedIndices, 2), 1);

% If verbose, display a waitbar
if verbose; f = waitbar(0, "Initializing..."); end

% For each window, compute the HRV
for i = 1:size(windowedIndices, 2)
    
    % If verbose, update the waitbar
    if verbose; waitbar(i/size(windowedIndices, 2), f, "Processing Window " ...
            + string(i) + " of " + string(size(windowedIndices, 2))); end
    
    % Remove nans, if 'ends' was selected
    if ends
        % Remove nans from windowedIndices
        temp = windowedIndices(:, i); window = temp(~isnan(temp));
    else
        window = windowedIndices(:, i);
    end
    
    % Format command based on HRV method selected
    if difference
        hrv(i) = cardio.ecg.computeHRV(window, Fs, 'difference', 'tol', tol);
    elseif spectral
        [hrv(i), LF(i), HF(i)] = cardio.ecg.computeHRV(window, Fs, 'spectral', 'tol', tol);
    elseif poincare
        [hrv(i), Cov1(i), Cov2(i)] = cardio.ecg.computeHRV(window, Fs, 'poincare', 'tol', tol);
    end
    
end

% Specify variable output arguments
if difference
    if nargout > 1; varargout{1} = []; end
    if nargout > 2; varargout{2} = []; end
elseif spectral
    if nargout > 1; varargout{1} = LF; end
    if nargout > 2; varargout{2} = HF; end
elseif poincare
    if nargout > 1; varargout{1} = Cov1; end
    if nargout > 2; varargout{2} = Cov2; end
end

% If verbose, close the waitbar
if verbose; close(f); end

end

