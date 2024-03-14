function obj = trueMAP(obj, varargin)

% -------------------------------------------------------------------------
% This function computes the mean arterial pressure in each of the four
% pressure catheters during each cardiac cycle.
%
% Arguments (opt'l)
% - levels      [Level] Levels for which to compute pressures
% - subjects            List of subjects for which to compute pressures
% - aortic      FLAG    Compute aortic MAP only?
% - femoral     FLAG    Compute femoral MAP only?
% - wedge       FLAG    Compute wedge pressure only?
% - ra          FLAG    Compute right atrial pressure only?
% - verbose     FLAG    Print progress?
%
% Returns
% - obj.dataset{subject}.(level).aorticMAP
% - obj.dataset{subject}.(level).femoralMAP
% - obj.dataset{subject}.(level).wedgeMAP
% - obj.dataset{subject}.(level).rightAtrialMAP
% - obj.dataset{subject}.(level).systolic   (if aorticMAP)
% - obj.dataset{subject}.(level).diastolic  (if aorticMAP)
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'levels'); levels = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'subjects'); subjects = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'aortic'); aortic = true;
        elseif strcmp(varargin{arg}, 'femoral'); femoral = true;
        elseif strcmp(varargin{arg}, 'wedge'); wedge = true;
        elseif strcmp(varargin{arg}, 'ra'); ra = true;
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        end
    end
end

% Set defaults for optional arguments
if ~exist('levels', 'var'); levels = obj.levels; end
if ~exist('subjects', 'var'); subjects = obj.subjects; end
if ~exist('aortic', 'var'); aortic = false; end
if ~exist('femoral', 'var'); femoral = false; end
if ~exist('wedge', 'var'); wedge = false; end
if ~exist('ra', 'var'); ra = false; end
if ~aortic && ~femoral && ~wedge &&~ra; ...
        allPressures = true; else; allPressures = false; end

% For each subject...
for subject = 1:length(subjects)
    
    % Print progress (if applicable)
    if verbose; disp("-> Processing Subject " + string(subject) + " of " + string(length(subjects))); end
    
    % For each level...
    for level = levels
        
        % If the level does not exist for the subject, continue
        if ~isfield(obj.dataset{subject}, string(level)); continue; end
        
        % Process the aortic pressure waveforms, if indicated
        if aortic || allPressures
            
            % Import waveforms
            aorticSignals = obj.dataset{subject}.(string(level)).aorticPressure;
            % Compute aortic MAP
            aorticMAP = mean(aorticSignals);
            
            % Compute systolic and diastolic pressures
            systolic = max(aorticSignals);
            diastolic = min(aorticSignals);
            
        end
        
        % Process the femoral pressure waveform, if indicated
        if femoral || allPressures
            
            % Import waveforms
            femoralSignals = obj.dataset{subject}.(string(level)).femoralPressure;
            % Compute femoral MAP
            femoralMAP = mean(femoralSignals);
            
        end
        
        % Process the wedge pressure waveform, if indicated
        if wedge || allPressures
            
            % Import waveforms
            wedgeSignals = obj.dataset{subject}.(string(level)).wedgePressure;
            % Compute mean wedge pressure
            wedgeMAP = mean(wedgeSignals);
            
        end
        
        % Process the right atrial pressure waveform, if indicated
        if ra || allPressures
            
            % Import waveforms
            raSignals = obj.dataset{subject}.(string(level)).rightAtrialPressure;
            % Compute mean right atrial pressure
            raMAP = mean(raSignals);
            
        end
        
        % Save output values
        if aortic || allPressures
            obj.dataset{subject}.(string(level)).aorticMAP = aorticMAP; 
            obj.dataset{subject}.(string(level)).systolic = systolic;
            obj.dataset{subject}.(string(level)).diastolic = diastolic;
        end
        if femoral || allPressures; obj.dataset{subject}.(string(level)).femoralMAP = femoralMAP; end
        if wedge || allPressures; obj.dataset{subject}.(string(level)).wedgeMAP = wedgeMAP; end
        if ra || allPressures; obj.dataset{subject}.(string(level)).rightAtrialMAP = raMAP; end
        
    end
    
end

end

