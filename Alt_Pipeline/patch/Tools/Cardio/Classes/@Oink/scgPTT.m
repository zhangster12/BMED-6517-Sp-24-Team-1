function obj = scgPTT(obj, varargin)

% -------------------------------------------------------------------------
% Description:
% This function extracts the PTT interval from SCG and PPG waveforms. The
% feature extraction method is optimized for waveforms that are assumed
% clean, as data was collected in this study in a controlled environment.
%
% Arguments (opt'l):
% - levels  [Level]     Levels for which to estimate PTT
% - subjects            List of subjects for which to estimate PTT
% - verbose FLAG        Print progress?
%
% Returns:
% - obj.dataset{subject}.(level).scgPTT
% - obj.dataset{subject}.(level).scgPAT
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
            ptt = obj.dataset{subject}.allRelative.scgPTT(beats(1):beats(end));
            pat = obj.dataset{subject}.allRelative.scgPAT(beats(1):beats(end));
            obj.dataset{subject}.(string(level)).scgPTT = ptt;
            obj.dataset{subject}.(string(level)).scgPAT = pat;
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
            ptt = obj.dataset{subject}.allAbsolute.scgPTT(beats(1):beats(end));
            pat = obj.dataset{subject}.allAbsolute.scgPAT(beats(1):beats(end));
            obj.dataset{subject}.(string(level)).scgPTT = ptt;
            obj.dataset{subject}.(string(level)).scgPAT = pat;
            continue
        end
        
        % Get the PPG data
        ppg = obj.dataset{subject}.(string(level)).femoralPPG;
        ppg = normalize(ppg);
        
        % Compute AO if it hasn't been computed already
        if ~isfield(obj.dataset{subject}.(string(level)), 'scgRAO')
            obj = obj.scgAortic('subjects', subject, 'levels', level, 'AO', 'verbose');
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
                    pat = cardio.general.pulseTime(ppg, 'simplediff2peak', 'min', 650, 'memory', 30);
                    
                case 2
                    
                    % Subject 2
                    pat = cardio.general.pulseTime(ppg, 'simplediff2peak', 'min', 500, 'max', 800, 'memory', 30);
                    
                case 3
                    
                    % Subject 3
                    pat = cardio.general.pulseTime(ppg, 'diff2peak', 'smooth', 200, 'min', 500, 'memory', 30);
                    
                case 4
                    
                    % Subject 4
                    pat = cardio.general.pulseTime(ppg, 'tangents', 'smooth', 10, 'min', 500, 'memory', 30);
                    
                case 5
                    
                    % Subject 5
                    pat = cardio.general.pulseTime(ppg, 'simplediff2peak', 'smooth', 20, 'min', 700, 'memory', 60);
                    
                case 6
                
                    % Subject 6
                    pat = cardio.general.pulseTime(ppg, 'simplediff2peak', 'min', 600, 'memory', 30);
                    
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
                    pat = cardio.general.pulseTime(ppg, 'simplediff2peak', 'min', 700, 'memory', 30);
                    
                case 2
                    
                    % Subject 2
                    pat = cardio.general.pulseTime(ppg, 'simplediff2peak', 'min', 500, 'max', 800, 'memory', 30);
                    
                case 3
                    
                    % Subject 3
                    pat = cardio.general.pulseTime(ppg, 'simplediff2peak', 'min', 500, 'memory', 30);
                    
                case 4
                    
                    % Subject 4
                    pat = cardio.general.pulseTime(ppg, 'simplediff2peak', 'smooth', 50, 'min', 500, 'max', 800, 'memory', 30);
                    
                case 5
                    
                    % Subject 5
                    pat = cardio.general.pulseTime(ppg, 'simplediff2peak', 'min', 500, 'memory', 30);
                    
                case 6
                
                    % Subject 6
                    pat = cardio.general.pulseTime(ppg, 'simplediff2peak', 'smooth', 25, 'min', 600, 'max', 900, 'memory', 50);
                    
            end
        
        end
        
        % Interpolate NaNs (linear)
        pat = fillmissing(pat, 'linear');
        
        % Convert to ms
        pat = (pat/obj.Fs)*1000;
        
        % Compute PTT from AO and PAT
        rao = obj.dataset{subject}.(string(level)).scgRAO;
        if length(pat) > length(rao); pat = pat(1:length(rao)); end
        ptt = pat - obj.dataset{subject}.(string(level)).scgRAO;
        
        % Save output
        obj.dataset{subject}.(string(level)).scgPTT = ptt;
        obj.dataset{subject}.(string(level)).scgPAT = pat;
        
    end
end

end