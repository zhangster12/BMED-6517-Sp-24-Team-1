function [windowed, varargout] = slideWindow(signal, windowRadius, overlap, varargin)

% -------------------------------------------------------------------------
% This function divides an input signal with a sliding window.
%
% Arguments (required)
% - signal          [Nx1]   Signal of length N
% - windowRadius            Radius of each window in samples: length/2 - 1
% - overlap                 Percent overlap of each window [0, 1]
%
% Arguments (optional)
% - nan             FLAG    Pad with NaNs (default zeros)
% - preserveEnds    FLAG    Preserve the ends of the signal by returning
%                           shortened windows centered on the end points
%
% Outputs
% - varargout{1}    [Mx2]   Start/stop indices for each of M windows
% -------------------------------------------------------------------------

% Parse optional inputs
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'nan'); nans = true;
        elseif strcmp(varargin{arg}, 'preserveEnds'); preserveEnds = true;
        end
    end
end

% Set defaults for optional inputs
if ~exist('nans', 'var'); nans = false; end
if ~exist('preserveEnds', 'var'); preserveEnds = false; end

% Ensure that the signal is a column vector
signal = signal(:);

% Initialize flags, placeholders, and counters
windowed = []; FLAG = true; idx = 1; counter = 0; 
windowLength = 2*windowRadius + 1; indices = [];

% Compute the index jump for the specified percent overlap
jump = round((1 - overlap)*windowLength);

% -------------------------------------------------------------------------
% Standard method (do not preserve end points)
% -------------------------------------------------------------------------

if ~preserveEnds

    % Perform the following for each window
    while FLAG

        % Increment the counter
        counter = counter + 1;
        
        % Determine whether the end of the signal has been reached
        if idx + (windowLength - 1) > length(signal)

            % If so, compute the buffer size and append to the final window
            bufferLength = windowLength - (length(signal) - idx) - 1;
            if ~nans; buffer = zeros(bufferLength, 1); else; buffer = nan*ones(bufferLength, 1); end
            finalWindow = [signal(idx:end); buffer];
            windowed = [windowed finalWindow];

            % Append window indices to index placeholder
            indices = [indices; idx length(signal)];
            
            % Throw the flag
            FLAG = false;

        else

            % Get the end index
            endIdx = idx + windowLength - 1;
            
            % Append the window to the return value
            windowed = [windowed signal(idx:endIdx)];
            
            % Append window indices to index placeholder
            indices = [indices; idx endIdx];

        end

        % Increment the index to obtain the desired overlap
        idx = idx + jump;

    end

end

% -------------------------------------------------------------------------
% Method of end preservation
% -------------------------------------------------------------------------

if preserveEnds
    
    for idx = 1:jump:length(signal)
        
        % Define the beginning and end points of the window, with "idx" as
        % the centroid of the window
        beginPoint = idx - windowRadius; endPoint = idx + windowRadius;
        
        % Limit the beginning and end points by the signal limits
        beginPoint = max(1, beginPoint);            % Beginning must be > 1
        endPoint = min(length(signal), endPoint);   % Ending must be < signal length
        
        % Append beginning and ending points to index placeholder
        indices = [indices; beginPoint endPoint];
        
        % Extract the window
        window = signal(beginPoint:endPoint);
        
        % If the window size is less than 2*radius + 1, append nans or 0's
        residual = (2*windowRadius + 1) - length(window);
        if residual > 0 && nans
            window = [window; nan*ones(residual, 1)];
        elseif residual > 0
            window = [window; zeros(residual, 1)];
        end
        
        % Append the window to the return value
        windowed = [windowed window];
        
    end
    
end

% Prepare output arguments
if nargout > 1; varargout{1} = indices; end

end

