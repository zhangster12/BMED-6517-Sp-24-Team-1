function obj = ecgHRV(obj, varargin)

% -------------------------------------------------------------------------
% This function computes heart rate variability using ECG via the
% difference method, spectral method, or Poincare method.
%
% Arguments (optional)
% - levels      [Level] Levels for which to compute HRV
% - subjects            List of subjects for which to compute HRV
% - difference	FLAG    Compute HRV using difference method (default)
% - spectral    FLAG    Compute HRV using spectral method
% - poincare    FLAG    Compute HRV using Poincare method
% - biopac      FLAG    Compute HRV using Biopac ECG (default)
% - t3          FLAG    Compute HRV using T3 ECG
% - winLength           Window length to compute HRV with sliding window
% - overlap             Percent overlap [0, 1] between windows
% - verbose     FLAG    Print progress?
% - tol                 Tolerance for outliers (in standard deviations)
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'levels'); levels = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'subjects'); subjects = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'difference'); difference = true;
        elseif strcmp(varargin{arg}, 'spectral'); spectral = true;
        elseif strcmp(varargin{arg}, 'poincare'); poincare = true;
        elseif strcmp(varargin{arg}, 'biopac'); biopac = true;
        elseif strcmp(varargin{arg}, 't3'); t3 = true;
        elseif strcmp(varargin{arg}, 'winLength'); winLength = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'overlap'); overlap = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        elseif strcmp(varargin{arg}, 'tol'); tol = varargin{arg + 1};
        end
    end
end

% Set defaults for optional input arguments
if ~exist('levels', 'var'); levels = obj.levels; end
if ~exist('subjects', 'var'); subjects = obj.subjects; end
if ~exist('difference', 'var'); difference = false; end
if ~exist('spectral', 'var'); spectral = false; end
if ~exist('poincare', 'var'); poincare = false; end
if ~exist('biopac', 'var'); biopac = false; end
if ~exist('t3', 'var'); t3 = false; end
if ~exist('verbose', 'var'); verbose = false; end
if ~difference && ~spectral && ~poincare; difference = true; end
if ~biopac && ~t3; biopac = true; end
if ~exist('winLength', 'var'); winLength = 100; end
if ~exist('overlap', 'var'); overlap = 0.5; end
if ~exist('tol', 'var'); tol = inf; end

% For each subject...
for subject = 1:length(subjects)
    
    % Print progress (if applicable)
    if verbose; disp("-> Processing Subject " + string(subject) + " of " + string(length(subjects))); end
    
    % For each level...
    for level = levels
        
        % If the level does not exist for the subject, continue
        if ~isfield(obj.dataset{subject}, string(level)); continue; end
        
        % Get the ECG beat indices
        if biopac; indices = obj.dataset{subject}.(string(level)).beats_biopac; ...
        elseif t3; indices = obj.dataset{subject}.(string(level)).beats_t3; end
        
        % Compute the HRV using the difference method, if indicated
        if difference
            if verbose
                hrv = cardio.ecg.slidingWindowHRV(indices, obj.Fs, 'difference', ...
                    'winLength', winLength, 'overlap', overlap, 'verbose', 'ends', 'tol', tol);
            else
                hrv = cardio.ecg.slidingWindowHRV(indices, obj.Fs, 'difference', ...
                    'winLength', winLength, 'overlap', overlap, 'ends', 'tol', tol);
            end
        end
        
        % Compute the HRV using the spectral method, if necessary
        if spectral
            % Compute HRV
            if verbose
                [~, LF, HF] = cardio.ecg.slidingWindowHRV(indices, obj.Fs, 'spectral', ...
                    'winLength', winLength, 'overlap', overlap, 'verbose', 'ends', 'tol', tol);
            else
                [~, LF, HF] = cardio.ecg.slidingWindowHRV(indices, obj.Fs, 'spectral', ...
                    'winLength', winLength, 'overlap', overlap, 'ends', 'tol', tol);
            end
            
            % Compute HRV
            hrv = LF./HF;
            
        end
        
        % Compute the HRV using the Poincare method, if necessary
        if poincare
            if verbose
                [~, Cov1, Cov2] = cardio.ecg.slidingWindowHRV(indices, obj.Fs, 'poincare', ...
                    'winLength', winLength, 'overlap', overlap, 'verbose', 'ends', 'tol', tol);
            else
                [~, Cov1, Cov2] = cardio.ecg.slidingWindowHRV(indices, obj.Fs, 'poincare', ...
                    'winLength', winLength, 'overlap', overlap, 'ends', 'tol', tol);
            end
            
            % Compute HRV
            hrv = Cov1./Cov2;
            
        end
        
        % Remove outliers that may occur at the beginning and end
        hrv(hrv > 2*median(hrv)) = nan; hrv = fillmissing(hrv, 'nearest');
        
        % Save the result
        if difference
            obj.dataset{subject}.(string(level)).hrv.difference = hrv;
        elseif spectral
            obj.dataset{subject}.(string(level)).hrv.spectral = hrv;
            obj.dataset{subject}.(string(level)).hrv.spectralLF = LF;
            obj.dataset{subject}.(string(level)).hrv.spectralHF = HF;
        else
            obj.dataset{subject}.(string(level)).hrv.poincare = hrv;
            obj.dataset{subject}.(string(level)).hrv.poincareEig1 = Cov1;
            obj.dataset{subject}.(string(level)).hrv.poincareEig2 = Cov2;
        end
        
    end
    
end

end

