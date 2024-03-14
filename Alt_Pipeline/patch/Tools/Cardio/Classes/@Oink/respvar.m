function obj = respvar(obj, varargin)

% -------------------------------------------------------------------------
% This function returns a beatwise estimate of respiratory variability in
% the specified signal.
%
% Arguments (optional)
% - levels      [Level]     Levels for which to compute variability
% - subjects                List of subjects to process
% - modalities  [Mode]      Modalities for which to compute variability
% - features    [Feature]   Features for which to compute variability
% - respRate                Respiration rate in cycles/min (default: 10)
% - varType                 Type of max/min comparison strategy
%   - 'SPV'                 SPV style: max - min
%   - 'PEPv'                PEPv style: (max - min)/(0.5*(max + min))
%   - 'PPV'                 PPV style (default)
% - no_interp   FLAG        Do not interpolate data (default: false)
% - verbose     FLAG        Print updates to console? (default: false)
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'levels'); levels = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'subjects'); subjects = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'modalities'); modalities = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'features'); features = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'respRate'); respRate = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'varType'); varType = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'no_interp'); no_interp = true;
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        end
    end
end

% Set defaults for optional input arguments
if ~exist('levels', 'var'); levels = obj.levels; end
if ~exist('subjects', 'var'); subjects = obj.subjects; end
if ~exist('modalities', 'var'); modalities = Mode.femoralPressure; end
if ~exist('features', 'var'); features = []; end
if ~exist('respRate', 'var'); respRate = 10; end
if ~exist('varType', 'var'); varType = 'PPV'; end
if ~exist('no_iterp', 'var'); no_interp = false; end
if ~exist('verbose', 'var'); verbose = false; end

% Combine modalities and features
modalities = [modalities features];

% Ensure that a valid varType is chosen
if ~strcmp(varType, 'SPV') && ~strcmp(varType, 'PEPv') && ~strcmp(varType, 'PPV')
    disp("Error in Oink.respvar(): Please enter a valid varType (SPV, PEPv, or PPV)")
    return
end

% For each subject...
for subject = 1:length(subjects)
    
    % Print progress, if applicable
    if verbose; disp("-> Processing Subject " + string(subject) + " of " + ...
            string(length(subjects))); end
    
    % For each level...
    for level = levels
        
        % If the level does not exist for the subject, continue
        if ~isfield(obj.dataset{subject}, string(level)); continue; end
        
        % Print progress, if applicable
        if verbose; disp("--> Processing Level: " + string(level)); end
        
        % For each modality...
        for mode = modalities
            
            % If the modality does not exist for the subject/level, continue
            if ~isfield(obj.dataset{subject}.(string(level)), string(mode)); continue; end

            % Extract the data for the current subject/level/modality
            signal = obj.dataset{subject}.(string(level)).(string(mode));
            
            % Get the heart rate
            HR = obj.dataset{subject}.(string(level)).HR_biopac;
            
            % -------------------------------------------------------------
            % Process Data <begin>
            % -------------------------------------------------------------
            
            % If SPV was selected, find the maxima of each signal segment
            if strcmp(varType, 'SPV'); signal = max(signal); end
            
            % Smooth the heart rate to remove ectopic beats and find the min
            smoothHR = movmean(HR, 10); minHR = min(smoothHR);
            
            % Find the minimum number of beats required (+ 25% error bound)
            minBeats = ceil(0.75*floor(minHR/respRate));
            
            % If PPV was selected, find the signal minima and maxima
            % for respiration cycle extraction
            if strcmp(varType, 'PPV'); minSig = movmean(min(signal), 3); ...
                    signal = movmean(max(signal), 5); end
            
            % Find peaks in the signal
            [~, maxLocs] = findpeaks(signal);
            [~, minLocs] = findpeaks(-signal);
            
            % Append by the length of the final peak/valley
            if length(maxLocs) > length(minLocs); minLen = length(minLocs); else; ...
                    minLen = length(maxLocs); end
            minLocs = minLocs(1:minLen); maxLocs = maxLocs(1:minLen);
            
            % Set placeholder for final maxima locations
            maxLocsKeep = [maxLocs(1)];
            
            % Determine which maxima to keep
            k = 2; last_k = length(maxLocs);
            beatFLAG = zeros(1, last_k); beatFLAG(1) = 1;
            while k < last_k
                while (maxLocs(k) - maxLocsKeep(end) > minBeats) && ((k + 1) < last_k)
                    maxLocsKeep = [maxLocsKeep maxLocs(k)];
                    beatFLAG(k) = 1; k = k + 1;
                end; k = k + 1;
            end
            
            % Determine minima to keep
            minLocsKeep = minLocs(logical(beatFLAG));
            
            % Evaluate the maxima and minima of the signal during the
            % maxima and minima of the respiratory cycles at the specified
            % peak locations
            if strcmp(varType, 'PPV')
                maxMax = signal(maxLocsKeep);	% Maximum of maxima
                maxMin = minSig(maxLocsKeep);   % Maximum of minima
                minMax = signal(minLocsKeep);	% Minimum of maxima
                minMin = minSig(minLocsKeep);   % Minimum of minima
                % Determine spread during respiratory cycle
                maxVal = maxMax - maxMin;
                minVal = minMax - minMin;
            else
                % Determine spread during respiratory cycle
                maxVal = signal(maxLocsKeep);
                minVal = signal(minLocsKeep);
            end
            
            % Compute the preload correlate based on the selected type
            if strcmp(varType, 'SPV'); preload = maxVal - minVal;
            else; preload = (maxVal - minVal)./(0.5*(maxVal + minVal));
            end
            
            % Interpolate the signal, if indicated
            if no_interp; result = preload;
            else
                interpBeats = [maxLocsKeep(1):maxLocsKeep(end)];
                result = [preload(1)];
                for i = 2:length(maxLocsKeep)
                    beats = [maxLocsKeep(i - 1):maxLocsKeep(i)];
                    range = preload(i) - preload(i - 1);
                    interval = range/(length(beats) - 1);
                    for j = 1:length(beats) - 1
                        result = [result preload(i-1) + j*interval];
                    end
                end
            end
            
            % -------------------------------------------------------------
            % Process Data <end>
            % -------------------------------------------------------------
            
            % Save the result
            obj.dataset{subject}.(string(level)).respvar.(string(mode)) = result;
            
        end
        
    end
    
end

end

