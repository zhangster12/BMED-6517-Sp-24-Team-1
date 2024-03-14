function obj = truePTT(obj, varargin)

% -------------------------------------------------------------------------
% Description:
% This function extracts the PTT interval from pressure waveforms. The
% feature extraction method is optimized for waveforms that are assumed
% clean, as data was collected in this study in a controlled environment.
%
% Arguments (opt'l):
% - levels  [Level]     Levels for which to estimate PTT
% - subjects            List of subjects for which to estimate PTT
% - verbose FLAG        Print progress?
%
% Returns:
% - obj.dataset{subject}.(level).truePTT
% - obj.dataset{subject}.(level).truePAT
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'levels'); levels = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'subjects'); subjects = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        end
    end
end

% Set defaults for optional arguments
if ~exist('levels', 'var'); levels = obj.levels; end
if ~exist('subjects', 'var'); subjects = obj.subjects; end
if ~exist('verbose', 'var'); verbose = false; end

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
            ptt = obj.dataset{subject}.allRelative.truePTT(beats(1):beats(end));
            pat = obj.dataset{subject}.allRelative.truePAT(beats(1):beats(end));
            obj.dataset{subject}.(string(level)).truePTT = ptt;
            obj.dataset{subject}.(string(level)).truePAT = pat;
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
            ptt = obj.dataset{subject}.allAbsolute.truePTT(beats(1):beats(end));
            pat = obj.dataset{subject}.allAbsolute.truePAT(beats(1):beats(end));
            obj.dataset{subject}.(string(level)).truePTT = ptt;
            obj.dataset{subject}.(string(level)).truePAT = pat;
            continue
        end
        
        % Get the femoral pressure data
        femoral = obj.dataset{subject}.(string(level)).femoralPressure;
        femoral = normalize(femoral);
        
        % Compute AO if it hasn't been computed already
        if ~isfield(obj.dataset{subject}.(string(level)), 'trueRAO')
            obj = obj.trueAortic('subjects', subject, 'levels', level, 'AO', 'verbose');
        end
        
        % -----------------------------------------------------------------
        % ABSOLUTE HYPOVOLEMIA
        % -----------------------------------------------------------------
        
        % Extract the PTT during absolute hypovolemia
        if level == Level.allAbsolute
            
            % Print progress if necessary
            if verbose; disp("--> Absolute Hypovolemia (PTT)"); end
        
            % Compute PTT for each subject
            switch obj.subjects(subject)
                
                case 1
                    
                    % Subject 1
                    pat = cardio.general.pulseTime(femoral, '%increase', 0.1, 'min', 400, 'max', 620);
                    
                case 2
                    
                    % Subject 2
                    pat = cardio.general.pulseTime(femoral, '%increase', 0.1, 'min', 275, 'max', 425);
                    
                case 3
                    
                    % Subject 3
                    pat = cardio.general.pulseTime(femoral, 'diff2peak', 'min', 300, 'memory', 30);
                    
                case 4
                    
                    % Subject 4
                    pat = cardio.general.pulseTime(femoral, '%increase', 0.1, 'min', 350, 'max', 500);
                    
                case 5
                    
                    % Subject 5
                    pat = cardio.general.pulseTime(femoral, '%increase', 0.1, 'min', 310, 'max', 450);
                    
                case 6
                
                    % Subject 6
                    pat = cardio.general.pulseTime(femoral, '%increase', 0.1, 'min', 400, 'max', 550);
                    
            end
        
        end
        
        % -----------------------------------------------------------------
        % RELATIVE HYPOVOLEMIA
        % -----------------------------------------------------------------
        
        % Extract the PTT during absolute hypovolemia
        if level == Level.allRelative
            
            % Print progress if necessary
            if verbose; disp("--> Relative Hypovolemia (PTT)"); end
        
            % Compute PTT for each subject
            switch obj.subjects(subject)
                
                case 1
                    
                    % Subject 1
                    pat = cardio.general.pulseTime(femoral, '%increase', 0.1, 'min', 450, 'max', 650);
                    
                case 2
                    
                    % Subject 2
                    pat = cardio.general.pulseTime(femoral, '%increase', 0.1, 'min', 300, 'max', 400);
                    
                case 3
                    
                    % Subject 3
                    pat = cardio.general.pulseTime(femoral, '%increase', 0.1, 'min', 300, 'max', 500);
                    
                case 4
                    
                    % Subject 4
                    pat = cardio.general.pulseTime(femoral, '%increase', 0.1, 'min', 350, 'max', 475);
                    
                case 5
                    
                    % Subject 5
                    pat = cardio.general.pulseTime(femoral, '%increase', 0.1, 'min', 350, 'max', 400);
                    
                case 6
                
                    % Subject 6
                    pat = cardio.general.pulseTime(femoral, '%increase', 0.1, 'min', 400, 'max', 550);
                    
            end
        
        end
        
        % Interpolate NaNs (linear)
        pat = fillmissing(pat, 'linear');
        
        % Convert to ms
        pat = (pat/obj.Fs)*1000;
        
        % Compute PTT from AO and PAT
        ptt = pat - obj.dataset{subject}.(string(level)).trueRAO;
        
        % Save output
        obj.dataset{subject}.(string(level)).truePTT = ptt;
        obj.dataset{subject}.(string(level)).truePAT = pat;
        
    end
end

end