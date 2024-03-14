function features = simpleFeatures(scg,  varargin)

% -------------------------------------------------------------------------
% Description:
% This function extracts consistent peaks from SCG signals which are
% relatively free from noise. A specified number of peaks are identified by
% the user, after which the same peak is selected in subsequent intervals.
%
% Arguments (req'd)
% - scg         [LxN]   N SCG signal segments of length L
% Arguments (opt'l)
% - memory (M)          Number of signal segments to store in memory (default 30)
%                       (And to label initially)
% - annotate    FLAG    Annotate the first #memory SCG segments
%                       (Default is automated annotation)
% - verbose     FLAG    Print progress to command window?
% - max         FLAG    Indicate whether the feature of interest is the maximum (default)
% - min         FLAG    Indicate whether the feature of interest is the minimum
% - seed        [Mx1]   Memory buffer seeded with first M features
% - from                Start index for feature selection
% - to                  End index for feature selection 1:30
% - filter      FLAG    Apply low-pass filter to remove extraneous peaks
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% Parse Inputs
% -------------------------------------------------------------------------

% Parse optional arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        switch varargin{arg}
            case 'memory'; memory = varargin{arg + 1};
            case 'annotate'; annotate = true;
            case 'verbose'; verbose = true;
            case 'max'; maximum = true;
            case 'min'; minimum = true;
            case 'seed'; seed = varargin{arg + 1};
            case 'from'; from = varargin{arg + 1};
            case 'to'; to = varargin{arg + 1};
            case 'filter'; filterFLAG = true;
        end
    end
end

% Set defaults for optional arguments if not provided
if ~exist('memory', 'var'); memory = 30; end
if ~exist('annotate', 'var'); annotate = false; end
if ~exist('verbose', 'var'); verbose = false; end
if ~exist('maximum', 'var'); maximum = false; end
if ~exist('minimum', 'var'); minimum = false; end
if ~maximum && ~minimum; maximum = true; end
if maximum && minimum; maximum = true; minimum = false; end
if ~exist('seed', 'var'); seed = []; end
if ~exist('from', 'var'); from = 1; end
if ~exist('to', 'var'); to = size(scg, 1); end
if ~exist('filterFLAG', 'var'); filterFLAG = false ;end

% Initialize memory buffer
buffer = zeros(memory, 1);
if ~isempty(seed); buffer(1:memory) = seed; end

% Set placeholder for return value
features = zeros(size(scg, 2), 1);

% Set filter
Hd = cardio.general.createBPF(1, 5, 2000, 'fpass1', 2, 'fpass2', 4, 'kaiser', 'order', 20);

%% ------------------------------------------------------------------------
% Signal Annotation
% -------------------------------------------------------------------------

% Annotate the first #memory signals if indicated.
% Else, annotate the signals automatically.
if annotate
    
    f = figure;
    
    % For each signal in the buffer...
    for seg = 1:memory
        
        % Extract the signal
        sig = scg(:, seg);
        
        % Differentiate and fileter the signal
        diffSig = diff(sig); if filterFLAG; diffSig = filtfilt(Hd.Numerator, 1, diffSig); end
        
        % Find zero-crossings of the first derivative
        crossings = [];     % Set placeholder for zero-crossings
        for sample = 1:length(diffSig)-1
            if diffSig(sample) > 0 && diffSig(sample + 1) < 0 || ...
                    diffSig(sample) < 0 && diffSig(sample + 1) > 0
                crossings = [crossings sample];
            end
        end
        
        % Plot the signal and zero-crossings
        plot(sig); hold on; title("Select Buffer: " + string(seg) + " of " + string(memory));
        for i = 1:length(crossings); plot(crossings(i), sig(crossings(i)), 'ok'); end
        
        % Accept user input
        [userInput, ~] = ginput(1);
        
        % Find the zero-crossing closest to the coordinate
        [~, crossIdx] = min(abs(crossings - userInput));
        
        % Write the value to the memory buffer
        buffer(seg) = crossings(crossIdx);
        
        % Remove hold on figure
        hold off
        
    end
    
    close(f)
    
end


%% ------------------------------------------------------------------------
% Feature Extraction
% -------------------------------------------------------------------------

% If there was no seed provided, seed the buffer
if ~annotate && isempty(seed)
    
    % For each buffer interval
    for seg = 1:memory
        
        % Extract signal
        sig = scg(from:to, seg);
        
        % Find minimum or maximum
        if minimum; [~, buffer(seg)] = min(sig); ...
        else; [~, buffer(seg)] = max(sig); end
        
    end
    
end

% Write buffer to placeholder
if ~annotate; buffer = buffer + from - 1; features(1:memory) = buffer; end

% For each remaining segment...
for seg = (memory + 1):size(scg, 2)
    
    % Extract the signal segment
    sig = scg(:, seg);
    
    % Initialize placeholder for candidate features
    candidates = [];
    
    % Differentiate and filter the signal
    diff1Sig = diff(sig); diff1Sig = filtfilt(Hd.Numerator, 1, diff1Sig);
    
    % Find proper zero-crossings
    for i = 1:length(diff1Sig) - 1
        if minimum && diff1Sig(i) < 0 && diff1Sig(i + 1) > 0
            candidates = [candidates i];
        elseif maximum && diff1Sig(i) > 0 && diff1Sig(i + 1) < 0
            candidates = [candidates i];
        end
    end
    
    % Remove candidates that are invalid
%     new_candidates = candidates;
%     new_candidates(new_candidates < from) = [];
%     new_candidates(new_candidates > to) = [];
%     if ~isempty(new_candidates); candidates = new_candidates; end
    
    % Compute the expected next beat from the buffer
    % trajectory = (buffer(end) - buffer(1))/memory;
    % expectedPoint = buffer(end) + trajectory;
    
    % Find the candidate feature closest to the buffer points
    [~, idx] = min(abs(candidates - mean(buffer))); newPoint = candidates(idx);
    % [~, idx] = min(abs(candidates - expectedPoint)); newPoint = candidates(idx);
    
    % Write the feature to the return vector and update buffer
    if isempty(newPoint); newPoint = nan; end
    features(seg) = newPoint; buffer = features((seg - memory + 1):seg);
    
end

end

