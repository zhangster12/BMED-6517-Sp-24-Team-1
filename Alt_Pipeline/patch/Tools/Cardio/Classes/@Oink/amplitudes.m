function obj = amplitudes(obj, varargin)

% -------------------------------------------------------------------------
% This function computes waveform amplitudes
%
% Arguments (optional)
% - levels      [Level] Levels for which to compute amplitudes
% - subjects            List of subjects for which to compute amplitudes
% - modalities  [Mode]  Field names in Oink object for which to compute amplitudes
% - verbose     FLAG    Print progress?
%
% Returns
% - obj.dataset{subject}.(level).amplitude.(Mode)
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'levels'); levels = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'subjects'); subjects = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'modalities'); modalities = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        end
    end
end

% Set defaults for optional arguments
if ~exist('levels', 'var'); levels = obj.levels; end
if ~exist('subjects', 'var'); subjects = obj.subjects; end
if ~exist('modalities', 'var'); modalities = [Mode.sternumSCG_z, Mode.femoralPPG]; end
if ~exist('verbose', 'var'); verbose = false; end

% For each subject
for subject = 1:length(subjects)
    
    % Print progress (if applicable)
    if verbose; disp("-> Processing Subject " + string(subject) + " of " + string(length(subjects))); end
    
    % For each level...
    for level = levels
        
        % If the level does not exist for the subject, continue
        if ~isfield(obj.dataset{subject}, string(level)); continue; end
        
        % For each mode...
        for mode = modalities
            
            % If the modality does not exist for the level, continue
            if ~isfield(obj.dataset{subject}.(string(level)), string(mode)); continue; end
            
            % Extract data for the current modality
            data = obj.dataset{subject}.(string(level)).(string(mode));
            
            % Obtain the amplitude of the each signal segment
            amplitude = max(data) - min(data);
            
            % Save the result
            obj.dataset{subject}.(string(level)).amplitude.(string(mode)) = amplitude;
            
        end
        
    end
    
end

end

