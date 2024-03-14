function varargout = selectTemplate(signals, varargin)

% -------------------------------------------------------------------------
% This function selects the most representative template in the signal set
% such that the selected signal minimizes the average Wasserstein distance
% between itself and the rest of the signal set after time-warping to
% adjust for natural morphological variability.
%
% Arguments (requried)
% - signals     [NxM]   M signals of length N
%
% Arguments (optional)
% - distance            Distance metric ('euclidian', 'dtw', or 'dtfm')
%                       (default: Wasserstein)
% - warp                Warp method ('dtw', 'dtfm', or 'none')
% - align       FLAG    Apply cross-correlation for signal alignment
% - Fs                  Sampling frequency
% - maxDist             If using DTFM warping, specify the max distance (in
%                       samples) between candidate features
% - verbose     FLAG    Print progress?
%
% Returns
% (1) rank        [Mx1]   Index of template with rank m (in) [1, M]
% (2) score       [Mx1]   Score of template with rank m (in) [1, M]
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'distance'); metric = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'warp'); warp = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'align'); align = true;
        elseif strcmp(varargin{arg}, 'Fs'); Fs = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'maxDist'); maxDist = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        end
    end
end

% Set defaults for optional input arguments
if ~exist('metric', 'var'); metric = 'wasserstein'; end
if ~exist('warp', 'var'); warp = 'dtw'; end
if ~exist('align', 'var'); align = false; end
if ~exist('Fs', 'var'); Fs = 2000; end
if ~exist('maxDist', 'var'); maxDist = 50; end
if ~exist('verbose', 'var'); verbose = false; end

% Get the number of signals
numSig = size(signals, 2);

% Normalize signals
signals = normalize(signals);

% Initialize distance matrix
distance = zeros(numSig);

% Initialize waitbar
if verbose; f = waitbar(0, "Selecting Template..."); end

% Initialize counter
counter = 0; total = numSig^2;

% For each template signal s_i
for i = 1:numSig
    
    % For each signal s_j
    for j = 1:numSig
        
        % Print progress, if indicated
        if verbose; waitbar(counter/total, f, "Processing Template " + string(i) + " of " + string(numSig)); end
        counter = counter + 1;
        
        % Skip if i and j are equal
        if i == j; continue; end
        
        % Extract the two signals
        s_i = signals(:, i); s_j = signals(:, j);
        
        % Align the signals, if indicated
        if align
            [s_i, s_j, ~, ~] = cardio.general.alignSignals(s_i, s_j, Fs, Fs, false);
        end
        
        % Warp the signals, if indicated
        switch warp
            case 'dtw'
                
                % Warp signals with DTW
                [dtwDistance, path_i, path_j] = dtw(s_i, s_j);
                s_i = s_i(path_i); s_j = s_j(path_j);
                
            case 'dtfm'
                
                % Warp signals with DTFM
                [dtfmDistance, temp_i, temp_j] = ...
                    cardio.sqi.peakMatch(s_j, s_i, 'maxDist', maxDist);
                s_i = temp_i{1}; s_j = temp_j{1};
                
        end
        
        % Compute the post- warping/alignment distance
        switch metric
            case 'dtw'
                distance(i, j) = dtwDistance;
            case 'dtfm'
                distance(i, j) = dtfmDistance;
            case 'euclidian'
                distance(i, j) = norm(s_j - s_i);
            otherwise
                [~, distance(i, j)] = cardio.sqi.wasserstein(s_j, s_i);
        end
        
    end
    
end

% Close the waitbar
if verbose; close(f); end

% Compute the average distance for each template
avgDistance = mean(distance, 2);

% Rank the signals
[score, rank] = sort(avgDistance, 'ascend');

% Return output arguments
if nargout > 0; varargout{1} = rank; end
if nargout > 1; varargout{2} = score; end

end