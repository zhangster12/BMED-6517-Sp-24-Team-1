function varargout = separateBeats(signals, varargin)

% -------------------------------------------------------------------------
% This function performs ensemble averaging and beat-separation of the 
% signal based on R-peaks of the corresponding ECG recorded concurrently.
% Input parameters are the following:
%
% Arguments (REQ'D)
% - signals     [MxN]   N signal vectors of length M to slice
%
% Arguments (OPT'L)
% - indices     [Dbl]   Indices at which to split signal
% - ecg         [Dbl]   Concurrent ECG signal with which to split signal
% - threshold           Threshold for R-peak detection
% - minDist             Minimum distance between ECG R peaks
% - samples             Number of samples per signal slice
% - backward    FLAG    Separate beats from back?
% - plot        FLAG    Plot ECG signal and sliced signals?
% - nanpad      BOOL    NaN-pad each beat based on its RR interval?
% - offset              Offset (in +/- samples)
%
% Output Arguments
% If ECG provided:
% --> varargout{1}      Sliced ECG signals
% --> varargout{2}      Sliced channel signals
% If indices provided:
% --> varargout         Sliced channel signals
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for i = 1:length(varargin)
        if strcmp(varargin{i}, 'indices'); indices = varargin{i + 1};
        elseif strcmp(varargin{i}, 'ecg'); ecg = varargin{i + 1};
        elseif strcmp(varargin{i}, 'threshold'); threshold = varargin{i + 1};
        elseif strcmp(varargin{i}, 'minDist'); minDist = varargin{i + 1};
        elseif strcmp(varargin{i}, 'samples'); samples = varargin{i + 1};
        elseif strcmp(varargin{i}, 'backward'); backward = true;
        elseif strcmp(varargin{i}, 'plot'); Plot = true;
        elseif strcmp(varargin{i}, 'nanpad'); nanpad = varargin{i + 1};
        elseif strcmp(varargin{i}, 'offset'); offset = varargin{i + 1};
        end
    end
end

% Set defaults
if ~exist('threshold', 'var'); threshold = 1; end
if ~exist('minDist', 'var'); minDist = 1; end
if ~exist('backward', 'var'); backward = false; end
if ~exist('Plot', 'var'); Plot = false; end
if ~exist('samples', 'var'); samples = size(signals, 1); end
if ~exist('indices', 'var'); indices = []; end
if ~exist('nanpad', 'var'); nanpad = false; end
if ~exist('offset', 'var'); offset = 0; end

% Determine whether to process an ECG signal
if exist('ecg', 'var'); processECG = true; else; processECG = false; end

% Determine how many signals are present and the signal length
numSignals = size(signals, 2); sigLength = size(signals, 1);

% Find peaks from filtered ECG if indices are not provided
if isempty(indices)
    [~, peaks] = findpeaks(ecg, 'MinPeakHeight', threshold, 'MinPeakDistance', minDist); 
else; peaks = indices;  % If indices are provided, separate at indices
end

% Figure: Labeled ECG R-Peaks
if Plot
    % Initialize subfigures and plot signals
    figure; subplot(numSignals + 1, 1, 1); hold on; grid on
    plot(sigLength, ecg); plot(peaks,ecg(peaks), 'rv', 'MarkerFaceColor', 'g')
    % Format subfigure
    title('ECG'); xlabel('Sample'); ylabel('ECG'); legend('ECG','R-peaks')
end

% Define empty matrices for holding slices
if processECG; sliced_ecg = zeros(samples, length(peaks)); end
sliced_signals = zeros(numSignals, samples, length(peaks));

% Convert peaks matrix to integers
peaks = int64(peaks);

% Populate the matrices
% For each peak, the "slice" includes a fixed #samples before/after that peak
for slice = 1:length(peaks)-1
    if processECG
        if ~backward
            if offset == 0
                sliced_ecg(:, slice) = ecg(peaks(slice):(peaks(slice) + samples - 1));
            else
                startIdx = max(1, peaks(slice) + offset);
                endIdx = min((startIdx + samples - 1), length(ecg));
                if endIdx == length(ecg); startIdx = endIdx - samples + 1; end
                sliced_ecg(:, slice) = ecg(startIdx:endIdx);
            end
        else
            if offset == 0
                sliced_ecg(:, slice) = ecg(peaks(slice+1) - samples + 1:peaks(slice+1));
            else
                endIdx = min(peaks(slice + 1) + offset, length(ecg));
                startIdx = max(1, endIdx - samples + 1);
                if startIdx == 1; endIdx = startIdx + samples - 1; end
                sliced_ecg(:, slice) = ecg(startIdx:endIdx);
            end
        end
        if nanpad
            rrint = peaks(slice + 1) - peaks(slice);
            sliced_ecg((rrint + 1):end, slice) = nan;
        end
    end
    for channel = 1:numSignals
        if ~backward
            if offset ~= 0
                startIdx = max(1, peaks(slice) + offset);
                endIdx = min((startIdx + samples - 1), length(signals(:, channel)));
                if endIdx == length(signals(:,channel)); startIdx = endIdx - samples + 1; end
                sliced_signals(channel, :, slice) = signals(startIdx:endIdx, channel);
                % sliced_signals(channel, :, slice) = signals(peaks(slice):(peaks(slice) + samples - 1), channel);
            else
                sliced_signals(channel, :, slice) = signals(peaks(slice):(peaks(slice) + samples - 1), channel);
            end
        else
            if offset ~= 0
                endIdx = min(peaks(slice + 1) + offset, length(signals(:, channel)));
                startIdx = max(1, endIdx - samples + 1);
                if startIdx == 1; endIdx = startIdx + samples - 1; end
                sliced_signals(channel, :, slice) = signals(startIdx:endIdx);
                % sliced_signals(channel, :, slice) = signals(peaks(slice+1) - samples + 1:peaks(slice+1), channel);
            else
                sliced_signals(channel, :, slice) = signals(peaks(slice+1) - samples + 1:peaks(slice+1), channel);
            end
        end
        if nanpad
            rrint = peaks(slice + 1) - peaks(slice);
            sliced_signals(channel, (rrint + 1):end, slice) = nan;
        end
    end
end

% Figure: Plot of sliced channels
if Plot
    for channel = 1:numSignals
        subplot(numSignals + 1, 1, channel + 1); hold on
        plot(sliced_signals(channel, :, :))
        xlabel('Sample'); ylabel('Amplitude')
        title("Channel " + string(channel)); grid on
    end
end

% Squeeze output if there is only one channel
if numSignals == 1; sliced_signals = squeeze(sliced_signals); end

% Format outputs
if processECG
    varargout{1} = sliced_ecg; varargout{2} = sliced_signals; varargout{3} = peaks;
else; varargout{1} = sliced_signals;
end

end