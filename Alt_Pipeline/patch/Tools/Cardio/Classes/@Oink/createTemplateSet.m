function obj = createTemplateSet(obj, varargin)

% -------------------------------------------------------------------------
% Description:
% This function creates TemplateSet objects for the desired modalities.
%
% Arguments (opt'l):
% - level       Level   Level from which SCG data is selected
% - subjects    [Nx1]	Subjects with which to generate TemplateSet
% - verbose     FLAG    Print progress to console?
% - lambda              Lambda parameter for DTFM
% - modalities  [Mode]  Field names in Oink object for which to generate TemplateSet
% - metric      Index   Distance metric for signal quality index
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'level'); level = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'subjects'); subjects = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        elseif strcmp(varargin{arg}, 'lambda'); lambda = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'modalities'); modalities = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'metric'); metric = varargin{arg + 1};
        end
    end
end

% Set default values for optional arguments
if ~exist('level', 'var'); level = Level.relBaseline1; end
if ~exist('subjects', 'var'); subjects = obj.subjects; end
if ~exist('verbose', 'var'); verbose = false; end
if ~exist('lambda', 'var'); lambda = 25; end
if ~exist('modalities', 'var'); modalities = [Mode.sternumSCG_z]; end
if ~exist('metric', 'var'); metric = Index.FMDTW; end

% Convert baseline interval to string
level = string(level);

%% ------------------------------------------------------------------------
% Create TemplateSet Object for each Modality
% -------------------------------------------------------------------------

for mode = modalities

    % Display progress, if indicated
    if verbose; disp("-> Processing Modality: " + string(mode)); end
    
    % Collect baseline SCG data from each subject
    data = cell(length(subjects), 1);
    for subject = 1:length(subjects)
        data{subject} = obj.dataset{subject}.(level).(string(mode));
    end

    % Initialize TemplateSet object
    templateSet = TemplateSet({data}, metric, lambda, obj.Fs);

    % Create template set
    if verbose; obj.templateSet.(string(mode)) = templateSet.create('verbose');
    else; obj.templateSet.(string(mode)) = templateSet.create();
    end

    % Clear the "data" field in the TemplateSet object (redundant)
    obj.templateSet.(string(mode)).data = [];

end

end