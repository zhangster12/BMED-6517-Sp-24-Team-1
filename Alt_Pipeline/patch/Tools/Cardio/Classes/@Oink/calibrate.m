function obj = calibrate(obj, varargin)

% -------------------------------------------------------------------------
% This function calibrates the subject's signals by obtaining a mapping
% between each metric and the aortic MAP. The mapping is:
% BP = C_1/PTT + C2
%
% Arguments (optional)
% - subjects                Subjects for which to calibrate signals
% - levels      [Level]     Levels for which to calibrate signals
% - features    [Feature]   Features for which to calibrate signals
%                           (Default: all PAT/PTT features)
%                           (Note: may be provided as string)
% - baseline    Level       Baseline for which to calibrate signals
%                           (Default: Level.allAbsolute)
% - smoothSignal            Smoothing factor for signal (default: 1)
% - smoothTarget            Smoothing factor for target (default: 1)
% - quad        FLAG        Quadratic mapping?
% - log         FLAG        Logarithmic mapping?
% - under       String      If the desired feature is a sub-feature (e.g.,
%                           hrvDifference is a sub-feature of hrv), specify
%                           the parent feature
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'subjects'); subjects = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'levels'); levels = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'features'); features = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'baseline'); baseline = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'smoothSignal'); sigM = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'smoothTarget'); tarM = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'quad'); quad = true;
        elseif strcmp(varargin{arg}, 'log'); logarithmic = true;
        elseif strcmp(varargin{arg}, 'under'); under = varargin{arg + 1};
        end
    end
end

% Set defaults for optional arguments
if ~exist('subjects', 'var'); subjects = obj.subjects; end
if ~exist('levels', 'var'); levels = obj.levels; end
if ~exist('features', 'var'); features = [Feature.scgPAT, Feature.scgPTT, ...
        Feature.truePAT, Feature.truePTT]; end
if ~exist('baseline', 'var'); baseline = Level.allAbsolute; end
if ~exist('sigM', 'var'); sigM = 1; end
if ~exist('tarM', 'var'); tarM = 1; end
if ~exist('quad', 'var'); quad = false; end
if ~exist('logarithmic', 'var'); logarithmic = false; end
if ~exist('under', 'var'); under = []; end

% Convert features to strings
if isa(features, 'Feature')
    temp = features; features = strings(size(temp));
    for i = 1:length(temp); features(i) = string(temp(i)); end
end

% For each subject...
for subject = 1:length(subjects)
    
    % Create a mapping placeholder for each feature
    mapping = cell(length(features), 1); counter = 1;
    
    % For each feature...
    for feature = features
        
        % If the baseline level does not exist for the subject, continue
        if ~isfield(obj.dataset{subject}, string(baseline)); continue; end
        
        % Get the data for the baseline level(s) for the current feature
        data = [];
        if length(baseline) < 2
            if isempty(under)
                data = obj.getSubinterval(feature, string(baseline), 'subjects', subject);
            else
                data = obj.getSubinterval(feature, string(baseline), 'subjects', subject, 'under', under);
            end
            if sigM > 1; data = movmedian(data, sigM); end
        else
            for i = 1:length(baseline)
                if isempty(under)
                    temp = obj.getSubinterval(feature, string(baseline(i)), 'subjects', subject);
                else
                    temp = obj.getSubinterval(feature, string(baseline(i)), 'subjects', subject, 'under', under);
                end
                if sigM > 1; temp = movmedian(temp, sigM); end; data = [data; temp(:)];
            end
        end
        
        % Remove outliers
        upperBound = mean(data) + 2*std(data); lowerBound = mean(data) - 2*std(data);
        data(data < lowerBound | data > upperBound) = nan; data = fillmissing(data, 'linear');
        
        % Convert data to a column and append column of 1's (offsets)
        data = data(:); data = 1./data;
        if ~quad && ~logarithmic
            data = [data ones(length(data), 1)];
        elseif ~logarithmic
            data = [data.^2 data ones(length(data), 1)];
        else
            data = [log(1./data) ones(length(data), 1)];
        end
        data(isinf(data)) = 0;  % Correct errors
        
        % Get the target data for the baseline level
        target = [];
        if length(baseline) < 2
            target = obj.getSubinterval('aorticMAP', string(baseline), 'subjects', subject);
            if tarM > 1; target = movmedian(target, tarM); end
        else
            for i = 1:length(baseline)
                temp = obj.getSubinterval('aorticMAP', string(baseline(i)), 'subjects', subject);
                if tarM > 1; temp = movmedian(temp, tarM); end; target = [target; temp(:)];
            end
        end
        target = target(:);
        
        % Compute the mapping
        mapping{counter} = data\target; counter = counter + 1;
        
    end
    
    % For each level...
    for level = levels
       
        % Initialize feature counter
        counter = 1;
        
        % For each feature...
        for feature = features
        
            % If the level does not exist for the subject, continue
            if ~isfield(obj.dataset{subject}, string(level)); continue; end

            % If the field does not exist for the level, continue
            if ~isfield(obj.dataset{subject}.(string(level)), feature); continue; end

            % Extract the data for the current level
            if isempty(under)
                data = obj.getSubinterval(feature, string(level), 'subjects', subject);
            else
                data = obj.getSubinterval(feature, string(level), 'subjects', subject, 'under', under);
            end
            if sigM > 1; data = movmedian(data, sigM); end
            
            % Append column of 1's (offsets)
            data = data(:); data = 1./data;
            if ~quad && ~logarithmic
                data = [data ones(length(data), 1)];
            elseif ~logarithmic
                data = [data.^2 data ones(length(data), 1)];
            else
                data = [log(1./data) ones(length(data), 1)];
            end
            
            % Compute the calibrated data
            calibrated = data*mapping{counter}; counter = counter + 1;
            
            % Save the calibrated data for the feature
            % Create feature name
            ftName = feature + "Calibrated";
            if isempty(under)
                obj.dataset{subject}.(string(level)).(ftName) = calibrated;
            else
                obj.dataset{subject}.(string(level)).(under).(ftName) = calibrated;
            end
        
        end
        
    end
    
end

end

