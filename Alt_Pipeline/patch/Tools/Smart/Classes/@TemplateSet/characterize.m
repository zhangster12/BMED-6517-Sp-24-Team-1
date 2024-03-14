function obj = characterize(obj, varargin)

% This function characterizes the template set.

if isempty(obj.SQI)
    disp("Error in TemplateSet.characterize():")
    disp(" -> Run TemplateSet.create() before characterizing set")
    return
end

% Parse optional arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'without'); ho = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        end
    end
end

% Set defaults
if ~exist('verbose', 'var'); verbose = false; end

if verbose; disp(" "); ...
        disp("Characterizing Ensemble"); disp("-----------------------"); end

% Set placeholders and constants
C = obj.numClasses; S = obj.numSets;
Sigma = cell(S, 1); Sigma(:) = {cell(C, 1)};    % Covariance of scores
mu = cell(S, 1); mu(:) = {zeros(C, C)};         % Mean values for each set

% Specify held-out sets if indicated
hi = 1:S; if exist('ho', 'var'); hi(ho) = []; end

% Construct covariance matrix and mean value matrix
% For each template set that is held in...
for train = hi
    
    if verbose; disp("-> Progress: " + string(100*(train-1)/length(hi)) + "%"); end
    
    % Set placeholder for SQIs
    SQI = cell(C, C); SQI(:) = {[]};
    
    % Combine the SQIs for all valid test sets
    for test = hi
        % For each template in the set...
        for t = 1:C
            % Compute the SQI at each true class
            for c = 1:C
                SQI{t, c} = [SQI{t, c}; obj.SQI{train}{test, t, c}];
            end
        end
    end
    
    % For each true class...
    for c = 1:C
        % Get the covariance
        Sigma{train}{c} = nancov([SQI{:, c}]);
        % For each template...
        for t = 1:C
            mu{train}(t, c) = nanmean(SQI{t, c});
        end
    end
    
end

% Write outputs
obj.mu = mu; obj.sigma = Sigma;

end