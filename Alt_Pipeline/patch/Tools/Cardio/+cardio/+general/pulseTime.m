function varargout = pulseTime(signals, varargin)

% -------------------------------------------------------------------------
% This function infers the arrival time of pulsatile signals that have been
% segmented using a concurrent ECG. Peak-finding methods include first- or
% second-derivative peaks before the global maximum, global maximum, or
% minimum before the global maximum. A re-sampling algorithm may be applied
% to the extracted features to eliminate discontinuities.
%
% Arguments (required)
% - signals     [MxN]   Array of N signal vectors of length M
%
% Arguments (optional)
% - 'diff1peak'     FLAG    Pull first derivative peak
% - 'diff2peak'     FLAG    Pull second derivative peak
% - 'globalMax'     FLAG    Pull global maximum
% - 'globalMin'     FLAG    Pull global minimum
% - 'threshold'             Pull index at which signal passes threshold
%                           (provide threshold)
% - '%increase'             Pull index at which signal has percent increase
%                           (provide percent increase)
% - 'simpleIncrease'        Pull index at which signal has percent increase
%                           (without limiting by global maximum)
% - 'simplediff2peak'       Pull second derivative peak (without limiting
%                           by global maximum)
% - 'simplediff1peak'       Pull first derivative peak (without limiting by
%                           global maximum)
% - '%decrease'             Pull index at which signal has percent decrease
%                           (provide percent decrease)
% - 'tangents'              Intersecting tangents method
% - 'memory'                Number of prior samples to incorporate in
%                           re-sampling algorithm
% - 'plot'          FLAG    Visualize Results?
% - 'min'                   Minimum index for pulse detection
% - 'max'                   Maximum index for pulse detection
% - 'smooth'                Smoothing factor for differentiation
% - 'forback'       FLAG    Run algorithm forwards and backwards (must have
%                           memory > 0)
% - 'annotate'      FLAG    Annotate first #memory peaks for buffer
% - 'init_buffer'   [1xM]   Initial buffer vector of length #memory
% FOR FILTERING:
% - 'Fs'                                    Sampling frequency
% - 'fstop1'/'fpass1'/'fpass2'/'fstop2'     Filter parameters
%
% Outputs
% - 1 output: pulse timing
% - 2 outputs: [pulse timing, buffer]
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'diff1peak'); method = 'diff1peak';
        elseif strcmp(varargin{arg}, 'diff2peak'); method = 'diff2peak';
        elseif strcmp(varargin{arg}, 'simplediff2peak'); method = 'simplediff2peak';
        elseif strcmp(varargin{arg}, 'simplediff1peak'); method = 'simplediff1peak';
        elseif strcmp(varargin{arg}, 'globalMax'); method = 'globalMax';
        elseif strcmp(varargin{arg}, 'globalMin'); method = 'globalMin';
        elseif strcmp(varargin{arg}, 'threshold'); method = 'threshold'; Threshold = varargin{arg + 1};
        elseif strcmp(varargin{arg}, '%increase'); method = '%increase'; percentIncrease = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'simpleIncrease'); method = 'simpleIncrease'; percentIncrease = varargin{arg + 1};
        elseif strcmp(varargin{arg}, '%decrease'); method = '%decrease'; percentDecrease = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'tangents'); method = 'tangents';
        elseif strcmp(varargin{arg}, 'memory'); memory = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'plot'); Plot = true;
        elseif strcmp(varargin{arg}, 'min'); minIdx = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'max'); maxIdx = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'Fs'); Fs = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'fstop1'); fstop1 = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'fpass1'); fpass1 = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'fpass2'); fpass2 = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'fstop2'); fstop2 = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'smooth'); len = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'forback'); forback = true;
        elseif strcmp(varargin{arg}, 'annotate'); annotate = true;
        elseif strcmp(varargin{arg}, 'init_buffer'); init_buffer = varargin{arg + 1};
        end
    end
end

% Set defaults for optional input arguments
if ~exist('method', 'var'); method = 'diff2peak'; end
if ~exist('memory', 'var'); memory = 0; end
if ~exist('Plot', 'var'); Plot = false; end
if ~exist('minIdx', 'var'); minIdx = 1; end
if ~exist('maxIdx', 'var'); maxIdx = size(signals, 1); end
if ~exist('Fs', 'var'); Fs = 2000; end
if ~exist('fstop1', 'var'); fstop1 = 1; end
if ~exist('fpass1', 'var'); fpass1 = fstop1 + 1; end
if ~exist('fstop2', 'var'); fstop2 = 5; end
if ~exist('fpass2', 'var'); fpass2 = fstop2 - 1; end
if ~exist('len', 'var'); len = 0; end
if ~exist('forback', 'var'); forback = false; end
if ~exist('annotate', 'var'); annotate = false; end

% Set placeholder for return value
numSignals = size(signals, 2); pulse = zeros(numSignals, 1); pulse_raw = pulse;

% Create a filter to smooth the data after each differentiation
% Standard Kaiser window BP filter |fstop1/fpass1-fpass2\fstop2|
Hd = cardio.general.createBPF(fstop1, fstop2, Fs, 'fpass1', fpass1, 'fpass2', fpass2, 'kaiser', 'order', 20);

% -------------------------------------------------------------------------
% Extract Initial Points
% -------------------------------------------------------------------------

% If an init_buffer has been provided, assign it to pulse and continue
if exist('init_buffer', 'var') && memory > 0
    pulse(1:memory) = init_buffer;
    
    % Else, if 'annotate' was selected, annotate the buffer
elseif annotate && memory > 0
    
    % Initialize figure
    f = figure;
    
    % For each signal in the buffer
    for i = 1:memory
        
        % Extract and plot the signal
        sig = signals(:, i); plot(sig); hold on; grid on;
        xlabel("Sample"); ylabel("Amplitude");
        title("Signal " + string(i) + " of " + string(memory))
        
        % Record the user's input in the pulse vector
        [pulse(i), ~] = ginput(1); clf
        
    end
    
    % Close the figure
    close(f);
    
    % Fill the init_buffer
    init_buffer = pulse(1:memory);
    
    % Else, fill the memory buffer automatically
elseif memory > 0
    
    % For each signal in the memory buffer
    for i = 1:memory
    
        % Extract signal
        sig = signals(:, i);

        % Select indicated feature
        switch method
            case 'diff1peak'
                [pulse(i), ~, ~] = diff1Peak(sig);
            case 'diff2peak'
                [pulse(i), ~, ~] = diff2Peak(sig);
            case 'simplediff2peak'
                [pulse(i), ~, ~] = simpleDiff2Peak(sig);
            case 'simplediff1peak'
                [pulse(i), ~, ~] = simpleDiff1Peak(sig);
            case 'globalMax'
                [pulse(i), ~, ~] = globalMax(sig);
            case 'globalMin'
                [pulse(i), ~, ~] = globalMin(sig);
            case 'threshold'
                [pulse(i), ~] = threshold(sig);
            case '%increase'
                [pulse(i), ~] = increase(sig);
            case 'simpleIncrease'
                [pulse(i), ~] = simpleIncrease(sig);
            case '%decrease'
                [pulse(i), ~] = decrease(sig);
            case 'tangents'
                [pulse(i), ~] = tangents(sig);
            otherwise
                [pulse(i), ~] = diff2Peak(sig);
        end

    end
    
    % Fill the init_buffer
    init_buffer = pulse(1:memory);
    
else
    
    % Leave the init_buffer empty if no memory is specified
    init_buffer = [];
    
end

% -------------------------------------------------------------------------
% Extract Subsequent Points (with Re-Sampling)
% -------------------------------------------------------------------------

for i = (memory+1):numSignals
    
    % Extract signal
    sig = signals(:, i);
    
    % Get candidates for indicated feature
    switch method
        case 'diff1peak'
            [bestGuess, ~, candidates] = diff1Peak(sig);
        case 'diff2peak'
            [bestGuess, ~, candidates] = diff2Peak(sig);
        case 'simplediff2peak'
            [bestGuess, ~, candidates] = simpleDiff2Peak(sig);
        case 'simplediff1peak'
            [bestGuess, ~, candidates] = simpleDiff1Peak(sig);
        case 'globalMax'
            [bestGuess, ~, candidates] = globalMax(sig);
        case 'globalMin'
            [bestGuess, ~, candidates] = globalMin(sig);
        case 'threshold'
            [bestGuess, candidates] = threshold(sig);
        case '%increase'
            [bestGuess, candidates] = increase(sig);
        case 'simpleIncrease'
            [bestGuess, candidates] = simpleIncrease(sig);
        case '%decrease'
            [bestGuess, candidates] = decrease(sig);
        case 'tangents'
            [bestGuess, candidates] = tangents(sig);
        otherwise
            [bestGuess, candidates] = diff2Peak(sig);
    end
    
    % If a memory buffer is specified
    if memory > 0 && ~isempty(candidates)
        % Select the candidate that is closest to the prior pulses in memory
        buffer = pulse(~isnan(pulse(1:i-1))); buffer = buffer(max(1, end-memory-1):end);
        dist = zeros(length(candidates), 1);
        % buffer = pulse(i-memory:i-1); dist = zeros(length(candidates), 1);
        for k = 1:length(dist); dist(k) = mean(buffer - candidates(k)); end
        [~, idx] = min(abs(dist)); pulse(i) = candidates(idx); pulse_raw(i) = bestGuess;
    else
        % Select the best guess as the return value
        if ~isempty(bestGuess)
            pulse(i) = bestGuess; pulse_raw(i) = bestGuess;
            % pulse(i) = nan; pulse_raw(i) = nan;
        else
            pulse(i) = nan; pulse_raw(i) = nan;
        end
    end
    
end

% -------------------------------------------------------------------------
% Run Analysis Backwards (to Correct Early Errors)
% -------------------------------------------------------------------------

if forback && memory > 0
    
    % Flip the signals and seed the pulse buffer
    backSignals = fliplr(signals); backPulse = flip(pulse);

    for i = (memory + 1):numSignals
        
        % Extract signal
        sig = backSignals(:, i);

        % Get candidates for indicated feature
        switch method
            case 'diff1peak'
                [bestGuess, ~, candidates] = diff1Peak(sig);
            case 'diff2peak'
                [bestGuess, ~, candidates] = diff2Peak(sig);
            case 'simplediff2peak'
                [bestGuess, ~, candidates] = simpleDiff2Peak(sig);
            case 'simplediff1peak'
                [bestGuess, ~, candidates] = simpleDiff1Peak(sig);
            case 'globalMax'
                [bestGuess, ~, candidates] = globalMax(sig);
            case 'globalMin'
                [bestGuess, ~, candidates] = globalMin(sig);
            case 'threshold'
                [bestGuess, candidates] = threshold(sig);
            case '%increase'
                [bestGuess, candidates] = increase(sig);
            case 'simpleIncrease'
                [bestGuess, candidates] = simpleIncrease(sig);
            case '%decrease'
                [bestGuess, candidates] = decrease(sig);
            case 'tangents'
                [bestGuess, candidates] = tangents(sig);
            otherwise
                [bestGuess, candidates] = diff2Peak(sig);
        end

        % If there are valida candidates...
        if ~isempty(candidates)
            % Select the candidate that is closest to the prior pulses in memory
            buffer = backPulse(~isnan(backPulse(1:i-1))); buffer = buffer(max(1, end-memory-1):end);
            dist = zeros(length(candidates), 1);
            for k = 1:length(dist); dist(k) = mean(buffer - candidates(k)); end
            [~, idx] = min(abs(dist)); backPulse(i) = candidates(idx);
        else
            % Select the best guess as the return value
            if ~isempty(bestGuess)
                backPulse(i) = bestGuess;
            else
                backPulse(i) = nan;
            end
        end

    end
    
    % Flip the back pulse to return the true pulse
    pulse = flip(backPulse);

end

% -------------------------------------------------------------------------
% Visualize Results
% -------------------------------------------------------------------------

% Visualize the results, if indicated
if Plot
    
    % Plot the resulting feature points on a static graph
    figure; hold on; grid on; title("Extracted Pulse Timing"); plot(pulse_raw); 
    plot(pulse); xlabel("Sample"); ylabel("Feature Value"); legend('Raw', 'Corrected');
    
    % Plot an animation of the resulting feature points
    cardio.general.heartbeatMovie(signals, 'features', pulse, 'speed', 'fast');
    
end

% -------------------------------------------------------------------------
% Format output
% -------------------------------------------------------------------------

% Check the number of output arguments
varargout{1} = pulse; if nargout > 1; varargout{2} = init_buffer; end

% -------------------------------------------------------------------------
% Sub-Functions for Feature Extraction
% -------------------------------------------------------------------------

    % First derivative peak (maximum and candidates)
    function [maximum, maxCandidate, candidates] = diff1Peak(signal)
        
        % Set placeholder for return value
        candidates = [];
        
        if len == 0
            % Differentiate signal and filter
            diff1Sig = diff(signal); diff1Sig = filtfilt(Hd.Numerator, 1, diff1Sig);
            diff2Sig = diff(diff1Sig); diff2Sig = filtfilt(Hd.Numerator, 1, diff2Sig);
        else
            % Differentiate signal and smooth
            diff1Sig = diff(signal); diff1Sig = movmean(diff1Sig, len);
            diff2Sig = diff(diff1Sig); diff2Sig = movmean(diff2Sig, len);
        end
        
        % Limit signal to region before global maximum
        [~, maximum] = max(signal);     % Get index of global maximum
        % Limit signal by global maximum
        if maximum <= length(diff2Sig); diff2Sig = diff2Sig(1:maximum); end
        if maximum <= length(diff1Sig); diff1Sig = diff1Sig(1:maximum); end
        
        % Find zero-crossings of second derivative
        for j = 1:length(diff2Sig)-1
            if diff2Sig(j) > 0 && diff2Sig(j+1) < 0; candidates = [candidates j]; end
        end
        
        % Limit valid candidates
        candidates(candidates > maxIdx) = [];
        candidates(candidates < minIdx) = [];
        
        % Determine the optimal candidate, if possible
        if ~isempty(candidates)
            [~, maxCandidate] = max(diff1Sig(candidates));
            maxCandidate = candidates(maxCandidate);
        else; maxCandidate = [];
        end
        
        % Return the maximum index
        [~, maximum] = max(diff1Sig);
        
    end

    % First derivative peak (without limiting by global maximum)
    function [maximum, maxCandidate, candidates] = simpleDiff1Peak(signal)
        
        % Set placeholder for return value
        candidates = [];
        
        if len == 0
            % Differentiate signal and filter
            diff1Sig = diff(signal); diff1Sig = filtfilt(Hd.Numerator, 1, diff1Sig);
            diff2Sig = diff(diff1Sig); diff2Sig = filtfilt(Hd.Numerator, 1, diff2Sig);
        else
            % Differentiate signal and smooth
            diff1Sig = diff(signal); diff1Sig = movmean(diff1Sig, len);
            diff2Sig = diff(diff1Sig); diff2Sig = movmean(diff2Sig, len);
        end
        
        % Find zero-crossings of second derivative
        for j = 1:length(diff2Sig)-1
            if diff2Sig(j) > 0 && diff2Sig(j+1) < 0; candidates = [candidates j]; end
        end
        
        % Limit valid candidates
        candidates(candidates > maxIdx) = [];
        candidates(candidates < minIdx) = [];
        
        % Determine the optimal candidate, if possible
        if ~isempty(candidates)
            [~, maxCandidate] = max(diff1Sig(candidates));
            maxCandidate = candidates(maxCandidate);
        else; maxCandidate = [];
        end
        
        % Return the maximum index
        [~, maximum] = max(diff1Sig);
        
    end

    % Second derivative peak (maximum and candidates)
    function [maximum, maxCandidate, candidates] = diff2Peak(signal)
        
        % Set placeholder for return value
        candidates = [];
        
        if len == 0
            % Differentiate signal and filter
            diff1Sig = diff(signal); diff1Sig = filtfilt(Hd.Numerator, 1, diff1Sig);
            diff2Sig = diff(diff1Sig); diff2Sig = filtfilt(Hd.Numerator, 1, diff2Sig);
            diff3Sig = diff(diff2Sig); diff3Sig = filtfilt(Hd.Numerator, 1, diff3Sig);
        else
            % Differentiate signal and smooth
            diff1Sig = diff(signal); diff1Sig = movmean(diff1Sig, len);
            diff2Sig = diff(diff1Sig); diff2Sig = movmean(diff2Sig, len);
            diff3Sig = diff(diff2Sig); diff3Sig = movmean(diff3Sig, len);
        end
        
        % Limit signal to region before global maximum
        [~, maximum] = max(signal);     % Get index of global maximum
        % Limit signal by global maximum
        if maximum <= length(diff3Sig); diff3Sig = diff3Sig(1:maximum); end
        if maximum <= length(diff2Sig); diff2Sig = diff2Sig(1:maximum); end
        
        % Find zero-crossings of third derivative
        for j = 1:length(diff3Sig)-1
            if diff3Sig(j) > 0 && diff3Sig(j+1) < 0; candidates = [candidates j]; end
        end
        
        % Limit valid candidates
        candidates(candidates > maxIdx) = [];
        candidates(candidates < minIdx) = [];
        
        % Determine the optimal candidate, if possible
        if ~isempty(candidates)
            [~, maxCandidate] = max(diff2Sig(candidates));
            maxCandidate = candidates(maxCandidate);
        else; maxCandidate = [];
        end
        
        % Return the maximum index
        [~, maximum] = max(diff2Sig);
        
    end

    % Second derivative peak (without limiting by global maximum)
    function [maximum, maxCandidate, candidates] = simpleDiff2Peak(signal)
        
        % Set placeholder for return value
        candidates = [];
        
        if len == 0
            % Differentiate signal and filter
            diff1Sig = diff(signal); diff1Sig = filtfilt(Hd.Numerator, 1, diff1Sig);
            diff2Sig = diff(diff1Sig); diff2Sig = filtfilt(Hd.Numerator, 1, diff2Sig);
            diff3Sig = diff(diff2Sig); diff3Sig = filtfilt(Hd.Numerator, 1, diff3Sig);
        else
            % Differentiate signal and smooth
            diff1Sig = diff(signal); diff1Sig = movmean(diff1Sig, len);
            diff2Sig = diff(diff1Sig); diff2Sig = movmean(diff2Sig, len);
            diff3Sig = diff(diff2Sig); diff3Sig = movmean(diff3Sig, len);
        end
        
        % Find zero-crossings of third derivative
        for j = 1:length(diff3Sig)-1
            if diff3Sig(j) > 0 && diff3Sig(j+1) < 0; candidates = [candidates j]; end
        end
        
        % Limit valid candidates
        candidates(candidates > maxIdx) = [];
        candidates(candidates < minIdx) = [];
        
        % Determine the optimal candidate, if possible
        if ~isempty(candidates)
            [~, maxCandidate] = max(diff2Sig(candidates));
            maxCandidate = candidates(maxCandidate);
        else; maxCandidate = [];
        end
        
        % Return the maximum index
        [~, maximum] = max(diff2Sig);
        
    end

    % Maxima (return global maximum and candidates)
    function [maximum, maxCandidate, candidates] = globalMax(signal)
        
        % Set placeholder for return value
        candidates = [];
        
        if len == 0
            % Differentiate signal and filter
            diff1Sig = diff(signal); diff1Sig = filtfilt(Hd.Numerator, 1, diff1Sig);
        else
            % Differentiate signal and smooth
            diff1Sig = diff(signal); diff1Sig = movmean(diff1Sig, len);
        end
        
        % Find zero-crossings of the first derivative
        for j = 1:length(diff1Sig)-1
            if diff1Sig(j) > 0 && diff1Sig(j+1) < 0; candidates = [candidates j]; end
        end
        
        % Limit valid candidates
        candidates(candidates > maxIdx) = [];
        candidates(candidates < minIdx) = [];
        
        % Determine the optimal candidate, if possible
        if ~isempty(candidates)
            [~, maxCandidate] = max(signal(candidates));
            maxCandidate = candidates(maxCandidate);
        else; maxCandidate = [];
        end
        
        % Find global maximum
        [~, maximum] = max(signal);
        
    end

    % Minima (return global minimum and candidates)
    function [minimum, minCandidate, candidates] = globalMin(signal)
        
        % Set placeholder for return value
        candidates = [];
        
        if len == 0
            % Differentiate signal and filter
            diff1Sig = diff(signal); diff1Sig = filtfilt(Hd.Numerator, 1, diff1Sig);
        else
            % Differentiate signal and smooth
            diff1Sig = diff(signal); diff1Sig = movmean(diff1Sig, len);
        end
        
        % Limit signal to region before global maximum
        [~, maximum] = max(signal);     % Get index of global maximum
        % Limit signal by global maximum
        if maximum <= length(diff1Sig); diff1Sig = diff1Sig(1:maximum); end
        if maximum <= length(signal); signal = signal(1:maximum); end
        
        % Find zero-crossings of the first derivative
        for j = 1:length(diff1Sig)-1
            if diff1Sig(j) < 0 && diff1Sig(j+1) > 0; candidates = [candidates j]; end
        end
        
        % Limit valid candidates
        candidates(candidates > maxIdx) = [];
        candidates(candidates < minIdx) = [];
        
        % Determine the optimal candidate, if possible
        if ~isempty(candidates)
            [~, minCandidate] = min(signal(candidates));
            minCandidate = candidates(minCandidate);
        else; minCandidate = [];
        end
        
        % Find global minimum
        [~, minimum] = min(signal);
        
    end

    % Threshold (index and candidates)
    function [index, candidates] = threshold(signal)
        
        % Set placeholder for return value
        candidates = [];
        
        % Limit signal to region before global maximum
        [~, maximum] = max(signal);     % Get index of global maximum
        signal = signal(1:maximum);     % Limit signal by global maximum
        
        % Find threshold crossings of the signal
        for j = 1:length(signal)-1
            if signal(j) < Threshold && signal(j+1) > Threshold; candidates = [candidates j]; end
        end
        
        % Limit valid candidates
        candidates(candidates > maxIdx) = [];
        candidates(candidates < minIdx) = [];
        
        % Find threshold crosing closest to the global maximum
        if ~isempty(candidates); index = candidates(end); else; index = []; end
        
    end

    % Percent increase (index and candidates)
    function [index, candidates] = increase(signal)
        
        % Set placeholder for return value
        candidates = [];
        
        % Get the signal maximum (of actual peaks)
        if len == 0
            % Differentiate signal and filter
            diff1Sig = diff(signal); diff1Sig = filtfilt(Hd.Numerator, 1, diff1Sig);
        else
            % Differentiate signal and smooth
            diff1Sig = diff(signal); diff1Sig = movmean(diff1Sig, len);
        end
        % Find zero-crossings of first derivative
        candidatePeaks = [];
        for j = 1:length(diff1Sig) - 1; if diff1Sig(j) > 0 && diff1Sig(j+1) < 0; ...
                    candidatePeaks = [candidatePeaks j]; end; end
        % Find highest point amongst zero-crossings
        [~, maxVal] = max(signal(candidatePeaks(candidatePeaks ~= 2)));
        maximum = candidatePeaks(maxVal);
        
        % Limit signal to region before global maximum
        % [~, maximum] = max(signal);     % Get index of global maximum
        signal = signal(1:maximum);     % Limit signal by global maximum
        
        % Set threshold value as a percentage of the signal range
        threshold = min(signal) + percentIncrease*range(signal);
        
        % Find candidate indices
        for j = 1:length(signal)-1
            if signal(j) < threshold && signal(j+1) > threshold; candidates = [candidates j]; end
        end
        
        % Limit valid candidates
        candidates(candidates > maxIdx) = [];
        candidates(candidates < minIdx) = [];
        
        % Find candidate closest to global maximum
        if ~isempty(candidates); index = candidates(end); else; index = nan; end
        
    end

    % Simple increase (do not limit by global maximum)
    function [index, candidates] = simpleIncrease(signal)
        
        % Set placeholder for return value
        candidates = [];
        
        % Set threshold value as a percentage of the signal range
        threshold = min(signal) + percentIncrease*range(signal);
        
        % Find candidate indices
        for j = 1:length(signal)-1
            if signal(j) < threshold && signal(j+1) > threshold; candidates = [candidates j]; end
        end
        
        % Limit valid candidates
        candidates(candidates > maxIdx) = [];
        candidates(candidates < minIdx) = [];
        
        % Find candidate closest to global maximum
        if ~isempty(candidates); index = candidates(end); else; index = nan; end
        
    end

    % Percent decrease (index and candidates)
    function [index, candidates] = decrease(signal)
        
        % Set placeholder for return value
        candidates = [];
        
        % Limit signal to region after global maximum
        [~, maximum] = max(signal);     % Get index of global maximum
        signal = signal(maximum:end);   % Limit signal by global maximum
        
        % Set threshold value as a percentage of the signal range
        threshold = max(signal) - percentDecrease*range(signal);
        
        % Find candidate indices
        for j = 1:length(signal)-1
            if signal(j) > threshold && signal(j+1) < threshold; candidates = [candidates j]; end
        end
        
        % Add max index back to candidates
        candidates = candidates + maximum;
        
        % Limit valid candidates
        candidates(candidates > maxIdx) = [];
        candidates(candidates < minIdx) = [];
        
        % Find candidate closest to global maximum
        if ~isempty(candidates); index = candidates(1); else; index = nan; end
        
    end

    % Intersecting tangents method
    function [index, candidates] = tangents(signal)
        
        % Set placeholder for return value
        candidates = [];
        
        % Find the signal minimum
        sigMin = min(signal);
        
        % Differentiate the signal
        diff1Sig = diff(signal);
        
        % Limit the signal to the region before the global maximum, if the
        % global maximum is after the minimum index
        [~, maximum] = max(signal);
        if maximum > minIdx && maximum > 60; signal = signal(1:maximum); end
        
        % Find the first derivative peak and candidates
        [fdPeak, ~, fdPeakCandidates] = simpleDiff1Peak(signal);
        
        % Limit valid candidates
        fdPeakCandidates(candidates > maxIdx) = [];
        fdPeakCandidates(candidates < minIdx) = [];
        
        % Find the slope of the line at each candidate point and return the
        % intersection with the minimum
        candidateSlopes = zeros(length(fdPeakCandidates));
        for j = 1:length(candidateSlopes)
            
            % Find the slope of the line at the candidate point
            candidateSlopes(j) = diff1Sig(fdPeakCandidates(j));
            
            % Find the intersection point
            intersection = (sigMin - (signal(fdPeakCandidates(j)) - ...
                candidateSlopes(j)*fdPeakCandidates(j)))/candidateSlopes(j);
            candidates = [candidates intersection];
            
        end
        
        % Limit valid candidates
        candidates(candidates > maxIdx) = [];
        candidates(candidates < minIdx) = [];
        
        % Find the best guess
        bestSlope = diff1Sig(fdPeak);
        index = (sigMin - (signal(fdPeak) - bestSlope*fdPeak))/bestSlope;
        
    end

end

