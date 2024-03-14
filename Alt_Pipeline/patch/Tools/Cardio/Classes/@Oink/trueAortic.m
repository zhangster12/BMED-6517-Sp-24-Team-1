function obj = trueAortic(obj, varargin)

% -------------------------------------------------------------------------
% Description:
% This function extracts the rAO/C interval from pressure waveforms. The
% feature extraction method is optimized for waveforms that are assumed
% clean, as data was collected in this study in a controlled environment.
%
% Arguments (opt'l):
% - levels  [Level]     Levels for which to estimate rAO/C interval
% - subjects            List of subjects for which to estimate rAO/C interval
% - AO      FLAG        Extract AO only?
% - AC      FLAG        Extract AC only?
% - verbose FLAG        Print progress?
%
% Returns:
% - obj.dataset{subject}.(level).trueRAO and/or
% - obj.dataset{subject}.(level).trueRAO
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'levels'); levels = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'subjects'); subjects = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'AO'); AO = true;
        elseif strcmp(varargin{arg}, 'AC'); AC = true;
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        end
    end
end

% Set defaults for optional arguments
if ~exist('levels', 'var'); levels = obj.levels; end
if ~exist('subjects', 'var'); subjects = obj.subjects; end
if ~exist('AO', 'var'); AO = false; end
if ~exist('AC', 'var'); AC = false; end
if ~exist('verbose', 'var'); verbose = false; end
if ~AO && ~AC; both = true; else; both = false; end

% If the list contains relative hypovolemia sub-levels, move
% Levels.allRelative to the front.
levels = levels(:)';
if ismember(Level.relBaseline1, levels) || ismember(Level.relBaseline2, levels) || ...
        ismember(Level.relative5, levels) || ismember(Level.relative10, levels) || ...
        ismember(Level.relative20, levels)
    levels(levels == Level.allRelative) = []; levels = [Level.allRelative levels];
end
% If the list contains absolute hypovolemia sub-levels, move
% Levels.allAbsolute to the front.
if ismember(Level.absBaseline1, levels) || ismember(Level.absBaseline2, levels) || ...
        ismember(Level.absDecrease7, levels) || ismember(Level.absDecrease14, levels) || ...
        ismember(Level.absDecrease21, levels) || ismember(Level.absDecrease28, levels) || ...
        ismember(Level.absIncrease7, levels) || ismember(Level.absIncrease14, levels) || ...
        ismember(Level.absIncrease21, levels) || ismember(Level.absIncrease28, levels)
    levels(levels == Level.allAbsolute) = []; levels = [Level.allAbsolute levels];
end

% For each subject...
for subject = 1:length(subjects)
    
    % Print progress (if applicable)
    if verbose; disp("-> Processing Subject " + string(subject) + " of " + string(length(subjects))); end
    
    % For each level...
    for level = levels
        
        % If the level does not exist for the subject, continue
        if ~isfield(obj.dataset{subject}, string(level)); continue; end
        
        % If the level is a relative hypovolemia sub-level, select the
        % appropriate beats from the pre-computed interval
        if level == Level.relBaseline1 || level == Level.relBaseline2 || ...
                level == Level.relative5 || level == Level.relative10 || ...
                level == Level.relative20
            beats = obj.dataset{subject}.(string(level)).indices_biopac;
            beats = beats - obj.dataset{subject}.allRelative.indices_biopac(1) + 1;
            rao = obj.dataset{subject}.allRelative.trueRAO(beats(1):beats(end));
            rac = obj.dataset{subject}.allRelative.trueRAC(beats(1):beats(end));
            obj.dataset{subject}.(string(level)).trueRAO = rao;
            obj.dataset{subject}.(string(level)).trueRAC = rac;
            continue
        end
        
        % If the level is an absolute hypovolemia sub-level, select the
        % appropriate beats from the pre-computed interval
        if level == Level.absBaseline1 || level == Level.absBaseline2 || ...
                level == Level.absDecrease7 || level == Level.absDecrease14 || ...
                level == Level.absDecrease21 || level == Level.absDecrease28 || ...
                level == Level.absIncrease7 || level == Level.absIncrease14 || ...
                level == Level.absIncrease21 || level == Level.absIncrease28
            beats = obj.dataset{subject}.(string(level)).indices_biopac;
            beats = beats - obj.dataset{subject}.allAbsolute.indices_biopac(1) + 1;
            rao = obj.dataset{subject}.allAbsolute.trueRAO(beats(1):beats(end));
            rac = obj.dataset{subject}.allAbsolute.trueRAC(beats(1):beats(end));
            obj.dataset{subject}.(string(level)).trueRAO = rao;
            obj.dataset{subject}.(string(level)).trueRAC = rac;
            continue
        end
        
        % Get the aortic pressure data (front and back split)
        aortic_front = obj.dataset{subject}.(string(level)).aorticPressure;
        aortic_back = obj.dataset{subject}.(string(level)).aorticPressure_back;
        aortic_front = normalize(aortic_front); aortic_back = normalize(aortic_back);
        
        % -----------------------------------------------------------------
        % ABSOLUTE HYPOVOLEMIA
        % -----------------------------------------------------------------
        
        % Extract the rAO interval, if indicated
        if (AO || both) && level == Level.allAbsolute
            
            % Print progress if necessary
            if verbose; disp("--> Absolute Hypovolemia (AO)"); end
        
            % Compute rAO interval for each subject
            switch obj.subjects(subject)
                
                case 1
                    
                    % Subject 1
                    rao = cardio.general.pulseTime(aortic_front, '%increase', 0.2, 'min', 150, 'max', 300, 'memory', 30);
                    
                case 2
                    
                    % Subject 2
                    rao = cardio.general.pulseTime(aortic_front, '%increase', 0.2, 'min', 125, 'max', 250, 'memory', 30);
                    
                case 3
                    
                    % Subject 3
                    rao = cardio.general.pulseTime(aortic_front, '%increase', 0.2, 'min', 175, 'max', 300, 'memory', 30);
                    
                case 4
                    
                    % Subject 4
                    rao = cardio.general.pulseTime(aortic_front, '%increase', 0.2, 'min', 150, 'max', 250, 'memory', 30);
                    
                case 5
                    
                    % Subject 5
                    rao = cardio.general.pulseTime(aortic_front, '%increase', 0.2, 'min', 125, 'max', 250, 'memory', 30);
                    
                case 6
                
                    % Subject 6
                    rao = cardio.general.pulseTime(aortic_front, '%increase', 0.2, 'min', 200, 'max', 300, 'memory', 30);
                    
            end
        
        end
        
        % Extract the rAC interval, if indicated
        if (AC || both) && level == Level.allAbsolute
        
            % Print progress, if necessary
            if verbose; disp("--> Absolute Hypovolemia (AC)"); end
            
            % Compute rAC for each subject
            switch obj.subjects(subject)
                case 1
                    
                    % Subject 1
                    rac = cardio.general.pulseTime(flipud(aortic_front), 'diff2peak', 'min', 150, 'memory', 30);
                    rac = size(aortic_front, 1) - rac;
                    
                case 2
                    
                    % Subject 2
                    rac = cardio.general.pulseTime(flipud(aortic_front), 'diff2peak', 'memory', 30);
                    rac = size(aortic_front, 1) - rac;
                    
                case 3
                    
                    % Subject 3
                    rac = cardio.general.pulseTime(flipud(aortic_back), 'diff2peak', 'memory', 30);
                    
                    % Compute beat lengths
                    bb = diff(obj.dataset{subject}.(string(level)).beats_t3);
                    beat_lengths = [bb(1) bb];
                    % Correct the rAC interval since back-beats were used
                    for i = 1:length(rac); rac(i) = beat_lengths(i) - rac(i); end
                    
                case 4
                    
                    % Subject 4
                    rac = cardio.general.pulseTime(aortic_front, '%decrease', 0.2, 'memory', 30);
                    
                case 5
                    
                    % Subject 5
                    rac = cardio.scg.simpleFeatures(aortic_front, 'min', 'from', 850, 'to', 950);
                    
                case 6
                    
                    % Subject 6
                    rac = cardio.general.pulseTime(flipud(aortic_front), 'min', 225, 'max', 450, 'memory', 30);
                    rac = size(aortic_front, 1) - rac;
                    
            end
        
        end
        
        % -----------------------------------------------------------------
        % RELATIVE HYPOVOLEMIA
        % -----------------------------------------------------------------
        
        % Extract the rAO interval, if indicated
        if (AO || both) && level == Level.allRelative
            
            % Print progress, if necessary
            if verbose; disp("--> Relative Hypovolemia (AO)"); end
            
            % Compute the rAO interval for each subject
            switch obj.subjects(subject)
                
                case 1
                    
                    % Subject 1
                    rao = cardio.general.pulseTime(aortic_front, '%increase', 0.2, 'min', 150, 'max', 300, 'memory', 30);
                    
                case 2
                    
                    % Subject 2
                    rao = cardio.general.pulseTime(aortic_front, '%increase', 0.2, 'min', 150, 'max', 300, 'memory', 30);
                    
                case 3
                    
                    % Subject 3
                    rao = cardio.general.pulseTime(aortic_front, '%increase', 0.2, 'min', 175, 'max', 300, 'memory', 30);
                    
                case 4
                    
                    % Subject 4
                    rao = cardio.general.pulseTime(aortic_front, '%increase', 0.2, 'min', 150, 'max', 300, 'memory', 30);
                    
                case 5
                    
                    % Subject 5
                    rao = cardio.general.pulseTime(aortic_front, '%increase', 0.2, 'min', 100, 'max', 200, 'memory', 30);
                    
                case 6
                    
                    % Subject 6
                    rao = cardio.general.pulseTime(aortic_front, '%increase', 0.2, 'min', 200, 'max', 300, 'memory', 30);
                    
            end
            
        end
        
        % Extract the rAC interval, if indicated
        if (AC || both) && level == Level.allRelative
            
            % Print progress, if necessary
            if verbose; disp("--> Relative Hypovolemia (AC)"); end
            
            % Compute ther rAC interval for each subject
            switch obj.subjects(subject)
                case 1
                    
                    % Subject 1
                    rac = cardio.general.pulseTime(flipud(aortic_front), 'diff2peak', 'memory', 30);
                    rac = size(aortic_front, 1) - rac;
                    
                case 2
                    
                    % Subject 2
                    rac = cardio.general.pulseTime(flipud(aortic_front), 'diff2peak', 'memory', 30);
                    rac = size(aortic_front, 1) - rac;
                    
                case 3
                    
                    % Subject 3
                    rac = cardio.general.pulseTime(flipud(aortic_back), 'diff2peak', 'memory', 30);
                    
                    % Compute beat lengths
                    bb = diff(obj.dataset{subject}.(string(level)).beats_t3);
                    beat_lengths = [bb(1) bb];
                    % Correct the rAC interval since back-beats were used
                    for i = 1:length(rac); rac(i) = beat_lengths(i) - rac(i); end
                    
                case 4
                    
                    % Subject 4
                    rac = cardio.general.pulseTime(flipud(aortic_back), 'diff2peak', 'memory', 30);
                    
                    % Compute beat lengths
                    bb = diff(obj.dataset{subject}.(string(level)).beats_t3);
                    beat_lengths = [bb(1) bb];
                    % Correct the rAC interval since back-beats were used
                    for i = 1:length(rac); rac(i) = beat_lengths(i) - rac(i); end
                    
                case 5
                    
                    % Subject 5
                    rac = cardio.general.pulseTime(flipud(aortic_front), 'diff2peak', 'min', 100, 'memory', 30);
                    rac = size(aortic_front, 1) - rac;
                    
                case 6
                    
                    % Subject 6
                    rac = cardio.general.pulseTime(flipud(aortic_front), 'diff2peak', 'min', 200, 'max', 350, 'memory', 30);
                    rac = size(aortic_front, 1) - rac;
                    
            end
            
        end
        
        % Interpolate NaNs (linear)
        if AO || both; rao = fillmissing(rao, 'linear'); end
        if AC || both; rac = fillmissing(rac, 'linear'); end
        
        % Convert to ms
        if AO || both; rao = (rao/obj.Fs)*1000; end
        if AC || both; rac = (rac/obj.Fs)*1000; end
        
        % Save output
        if AO || both; obj.dataset{subject}.(string(level)).trueRAO = rao; end
        if AC || both; obj.dataset{subject}.(string(level)).trueRAC = rac; end
        
    end
end

end