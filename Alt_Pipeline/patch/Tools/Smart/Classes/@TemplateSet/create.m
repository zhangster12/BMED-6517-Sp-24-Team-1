function obj = create(obj, varargin)

% -------------------------------------------------------------------------
% Description:
% This function constructs and characterizes a template set.
%
% Arguments (opt'l)
% - compute FLAG    Compute SQI?
% - verbose FLAG	Print progress to console?
% - char    FLAG	Re-characterize model?
% -------------------------------------------------------------------------

% Parse optional arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'verbose'); verbose = true;
        elseif strcmp(varargin{arg}, 'compute'); compute = true;
        elseif strcmp(varargin{arg}, 'char'); char = true;
        end
    end
end

% Set defaults
if ~exist('verbose', 'var'); verbose = false; end
if ~exist('char', 'var'); char = false; end
if ~exist('compute', 'var'); compute = false; end

% Process input argument
C = length(obj.data); S = length(obj.data{1});
obj.numClasses = C; obj.numSets = S;

% -------------------------------------------------------------------------
% Construct Template Sets
% -------------------------------------------------------------------------

% Display updates if indicated
if verbose; disp("Creating Templates"); disp("------------------"); end

% Set placeholder for return value
templates = cell(C, 1); templates(:) = {cell(S, 1)}; counter = 0;

% For each class...
for c = 1:C
    % For each set...
    for s = 1:S
        if verbose; disp("Progress: " + string(100*counter/(C*S)) + "%"); end
        % Construct a template and assign to proper place
        [~, temp] = cardio.general.template(obj.data{c}{s}, obj.Fs);
        templates{c}{s} = temp(:, end); counter = counter + 1;
    end
end

% -------------------------------------------------------------------------
% Evaluate Template Sets
% -------------------------------------------------------------------------

if compute

    % Display updates if indicated
    if verbose; disp("Evaluating Sets"); disp("------------------"); end

    % Set placeholders
    SQI = cell(S, 1); SQI(:) = {cell(S, C, C)};     % Signal quality index

    % Perform the analysis for each set
    for train = 1:S

        if verbose; disp("Progress: " + string(100*(train - 1)/S) + "%"); end

        % For each test set...
        for t = 1:C
            % For each template in the set...
            for ts = 1:S
                % Compute the SQI at each true class
                for c = 1:C
                    SQI{train}{ts, t, c} = cardio.sqi.filter(obj.data{c}{ts}, ...
                            obj.Fs, inf, 'template', templates{t}{train}, 'lambda', obj.lambda, 'type', obj.type);
                end
            end
        end

    end

    if verbose; disp("Progress: 100%"); end
    obj.SQI = SQI;

end

% -------------------------------------------------------------------------
% Organize Return Values and Characterize Set
% -------------------------------------------------------------------------

% Return templates
obj.set = templates;

% Characterize set, if indicated
if char; if verbose; obj = obj.characterize('verbose'); else; obj = obj.characterize(); end; end

end