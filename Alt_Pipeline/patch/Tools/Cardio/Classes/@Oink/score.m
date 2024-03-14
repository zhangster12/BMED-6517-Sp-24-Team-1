function obj = score(obj, varargin)

% -------------------------------------------------------------------------
% Description:
% This function generates SQI scores for each segment in the Oink object.
% If there is no TemplateSet object already, one is created.
%
% Arguments (opt'l):
% - verbose     FLAG    Print progress to console?
% - levels      [Level] Levels for which to generate SQI scores
% - subjects            List of subjects for which to generate SQI scores
% - modalities  [Mode]  Modalities for which to generate SQI scores
%
% Returns:
% - obj.dataset{subject}.(level).sqi.(modality)
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'verbose'); verbose = true;
        elseif strcmp(varargin{arg}, 'levels'); levels = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'subjects'); subjects = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'modalities'); modalities = varargin{arg + 1};
        end
    end
end

% Set defaults for optional arguments
if ~exist('verbose', 'var'); verbose = false; end
if ~exist('levels', 'var'); levels = obj.levels; end
if ~exist('subjects', 'var'); subjects = obj.subjects; end
if ~exist('modalities', 'var'); modalities = [Mode.sternumSCG_z]; end

% For each modality...
for mode = modalities

    % Print progress, if indicated
    if verbose; disp("-> Processing Modality: " + string(mode)); end
    
    % If TemplateSet object does not exist, create it
    if ~isfield(obj.templateSet, string(mode))
        % Print to console, if indicated, and create object
        if verbose
            disp("Creating TemplateSet object...")
            obj = obj.createTemplateSet('verbose', 'subjects', subjects, 'levels', levels, 'modalities', mode);
        else
            obj = obj.createTemplateSet('subjects', subjects, 'levels', levels, 'modalities', mode);
        end
    end

    % Generate SQI scores for each subject
    for subject = 1:length(subjects)
        % Generate SQI scores for each level
        for level = levels

            % If the level does not exist for the subject, continue
            if ~isfield(obj.dataset{subject}, string(level)); continue; end

            % Print update to console (if indicated)
            if verbose; disp("--> Generating SQIs for Level: " + string(level)); end

            % Extract data
            data = normalize(obj.dataset{subject}.(string(level)).(string(mode)));
            
            % Generate scores
            if verbose
                scores = obj.templateSet.(string(mode)).score('data', {data}, 'verbose');
            else
                scores = obj.templateSet.(string(mode)).score('data', {data});
            end

            % Save scores
            obj.dataset{subject}.(string(level)).sqi.(string(mode)) = ...
                obj.templateSet.(string(mode)).meanScore(scores);

        end
    end

end

end

