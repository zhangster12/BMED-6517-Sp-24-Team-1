function [data, labels] = trainingData(obj, varargin)

% -------------------------------------------------------------------------
% This function produces a matrix of training features and corresponding
% labels indicating blood volume status.
%
% Arguments (optional)
% - subjects        [Sx1]       List of subjects to include in training set
% - levels          [Level]     List of levels to include in training set
% - features        [Feature]   List of features to include in training set
% - labelMethod     Label       Method for generating labels (default .raw)
% - sampleMethod    Sample      Method for generating samples
% - normalize                   Normalize feature values against baseline?
%                               Specify level (e.g. Level.absBaseline1)
% - average                     Number of samples to average per datapoint
% - interval        [Int, Int]  Starting sample and space between each
%                               sample (if Sample.fixedInterval selected)
% - numSamples                  Number of samples for Samples.randomize
% - verbose         FLAG        Display progress?
% <<< Rolling Window Outlier Removal >>
% (To disable, set tol = inf)
% - winRad                      Window radius (default 500)
% - overlap                     Percent overlap in [0, 1] (default 0.5)
% - tol                         Tolerance in stddevs (default 2)
%
% Returns
% - data            [MxF]       Array of M values for F features
% - labels          [Mx1]       Vector of labels for each of M samples
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        try
            switch varargin{arg}
                case 'subjects'; subjects = varargin{arg + 1};
                case 'levels'; levels = varargin{arg + 1};
                case 'features'; features = varargin{arg + 1};
                case 'labelMethod'; labelMethod = varargin{arg + 1};
                case 'sampleMethod'; sampleMethod = varargin{arg + 1};
                case 'normalize'; normalize = true; baselineLevel = varargin{arg + 1};
                case 'average'; average = varargin{arg + 1};
                case 'interval'; interval = varargin{arg + 1};
                case 'numSamples'; numSamples = varargin{arg + 1};
                case 'verbose'; verbose = true;
                case 'winRad'; winRad = varargin{arg + 1};
                case 'overlap'; overlap = varargin{arg + 1};
                case 'tol'; tol = varargin{arg + 1};
            end
        catch
            continue
        end
    end
end

% Set default values for optional arguments
if ~exist('subjects', 'var'); subjects = obj.subjects; end
if ~exist('levels', 'var'); levels = obj.levels; end
if ~exist('features', 'var'); features = Feature.all; end
if ~exist('labelMethod', 'var'); labelMethod = Label.raw; end
if ~exist('sampleMethod', 'var'); sampleMethod = Sample.allSamples; end
if ~exist('normalize', 'var'); normalize = false; end
if ~exist('average', 'var'); average = 0; end
if ~exist('verbose', 'var'); verbose = false; end
if sampleMethod == Sample.fixedInterval && ~exist('interval', 'var'); interval = [1, 1]; end
if sampleMethod == Sample.randomized && ~exist('numSamples', 'var'); numSamples = 100; end
if ~exist('winRad', 'var'); winRad = 500; end
if ~exist('overlap', 'var'); overlap = 0.5; end
if ~exist('tol', 'var'); tol = 2; end

% If Feature.all was selected, construct the list of all features
if features == Feature.all
    features = enumeration('Feature'); features(features == Feature.all) = [];
end

% Ensure that the baseline level is the first level processed
levels = levels(:)';
if normalize; levels = [baselineLevel levels(levels ~= baselineLevel)]; end

% Display progress, if indicated
if verbose; f = waitbar(0, "Extracting Parent Data..."); end

% Create global variables for Level.allAbsolute and Level.allRelative data
absoluteData = cell(length(subjects), 1); relativeData = cell(length(subjects), 1);

% If the levels contain a Level.allAbsolute/Relative sub-level, extract the parent
for parent = [Level.allAbsolute, Level.allRelative]
    if ~isempty(find(obj.getParent(levels) == parent, 1))

        % For each subject...
        for subject = 1:length(subjects)
            
            % If the parent for this level does not exist, continue
            if ~isfield(obj.dataset{subject}, string(parent)); continue; end

            % Initialize a placeholder for the data
            temp = [];

            % For each feature...
            for feature = features

                % Extract the feature vector
                ftVector = obj.extractFeature(feature, parent, 'subjects', subject);

                % Ensure the feature vector is not too long; if it is too
                % short, append a NaN on the end
                if ~isempty(temp) && (length(ftVector) > size(temp, 1))
                    ftVector = ftVector(1:size(temp, 1));
                elseif ~isempty(temp) && (length(ftVector) < size(temp, 1))
                    diffLength = size(temp, 1) - length(ftVector);
                    ftVector = [ftVector(:); nan*ones(diffLength, 1)];
                end
                
                % Apply rolling filter to the data
                ftVector = cardio.general.rollingFilter(ftVector, ...
                    'winRad', winRad, 'overlap', overlap, 'tol', tol);
                
                % Enforce the averaging protocol
                if average > 0
                    ftVector = movmean(ftVector, average);
                end
                
                % Append the feature vector to the placeholder
                temp = [temp ftVector(:)];
                
            end

            % Assign the placeholder to the return value
            if parent == Level.allAbsolute; absoluteData{subject} = temp; end
            if parent == Level.allRelative; relativeData{subject} = temp; end

        end

    end
end

% Define placeholders for return values
data = []; labels = []; baselineValues = zeros(1, length(features));

% Initialize counter (for progress bar)
counter = 0;

% For each specified level, get the proper sub-interval of the parent
% For each subject...
for subject = 1:length(subjects)
    
    % Initialize placeholder for subject data
    subject_data = [];
    
    % For each level...
    for level = levels
        
        % Determine the parent for this level
        parent = obj.getParent(level);
        
        % Extract the parent's data
        if isfield(obj.dataset{subject}, string(parent))
            if parent == Level.allAbsolute; parentData = absoluteData{subject}; end
            if parent == Level.allRelative; parentData = relativeData{subject}; end
        end
        
        % If the current level does not exist for this subject, continue
        if ~isfield(obj.dataset{subject}, string(level)); continue; end
        
        % Initialize placeholder for subject + level data
        subject_level_data = [];
        
        % For each feature...
        for ft = 1:length(features)
            
            % Update waitbar, if necessary
            if verbose
                % Compute Progress
                progress = counter/(length(subjects)*length(labels)*length(features));
                % Update waitbar and message
                waitbar(progress, f, "Processing Subject " + string(subject) + ...
                    " of " + string(length(subjects))); counter = counter + 1;
            end
            
            % Extract the feature
            if isfield(obj.dataset{subject}, string(parent))
                ftVector = obj.getSubinterval([], level, 'of', parentData(:, ft));
            else
                ftVector = obj.extractFeature(feature, level, 'subjects', subject);
            end
            
            % Apply rolling filter to the data
            ftVector = cardio.general.rollingFilter(ftVector, ...
                'winRad', winRad, 'overlap', overlap, 'tol', tol);
            
            % Enforce the averaging protocol
            if average > 0
                ftVector = movmean(ftVector, average);
            end
            
            % Ensure the feature vector is not too long; if it is too
            % short, append a NaN on the end
            if ~isempty(subject_level_data) && (length(ftVector) > size(subject_level_data, 1))
                ftVector = ftVector(1:size(subject_level_data, 1));
            elseif ~isempty(subject_level_data) && (length(ftVector) < size(subject_level_data, 1))
                diffLength = size(subject_level_data, 1) - length(ftVector);
                ftVector = [ftVector(:); nan*ones(diffLength, 1)];
            end
            
            % Write the feature to the subject + level matrix
            subject_level_data = [subject_level_data ftVector(:)];
            
        end
        
        % If the current level is the baseline level, find the mean value
        % for each feature
        if normalize && level == baselineLevel
            % Record the mean values for each feature
            for ft = 1:length(features)
                baselineValues(ft) = mean(subject_level_data(:, ft));
            end
        end
        
        % Append subject + level to subject matrix
        subject_data = [subject_data; subject_level_data];
        
        % Append labels to label vector
        labels = [labels; repmat(level, size(subject_level_data, 1), 1)];
        
    end
    
    % Adjust the subject data by the baseline data
    if normalize
        subject_data = (subject_data - baselineValues)./baselineValues;
    end
    
    % Enforce the sampling protocol for the subject
    switch sampleMethod
        case Sample.fixedInterval
            
            % Set temporary placeholder
            temp = [];
            
            % Re-sample the data for the subject
            idx = interval(1);
            while idx <= size(subject_data, 1)
                temp = [temp; subject_data(idx, :)];
                idx = idx + interval(2);
            end
            
            % Update subject_data
            subject_data = temp;
            
        case Sample.randomized
            
            % Select a random permutation of samples
            perm = randperm(size(subject_data, 1), numSamples);
            
            % Update subject_data
            subject_data = subject_data(perm, :);
            
    end
    
    % Append subject matrix to data matrix
    data = [data; subject_data];
    
end

% Enforce the labeling protocol for the dataset
switch labelMethod
    case Label.bloodVolume
        
        % Set a placeholder of integers the same size as "labels"
        temp = zeros(size(labels));     % Baseline levels are 0
        
        % Set 7% levels to 1
        temp(labels == Level.absDecrease7 | labels == Level.absIncrease7) = 1;
        
        % Set 14% levels to 2
        temp(labels == Level.absDecrease14 | labels == Level.absIncrease14) = 2;
        
        % Set 21% levels to 3
        temp(labels == Level.absDecrease21 | labels == Level.absIncrease21) = 3;
        
        % Set 28% levels to 4
        temp(labels == Level.absDecrease28 | labels == Level.absIncrease28) = 4;
        
        % Set relative 5% BP drop to 5
        temp(labels == Level.relative5) = 5;
        
        % Set relative 10% BP drop to 6
        temp(labels == Level.relative10) = 6;
        
        % Set relative 20% BP drop to 7
        temp(labels == Level.relative20) = 7;
        
        % Write temporary placeholder to "labels"
        labels = temp;
        
    case Label.relativeVolume % Absolute hypovolemia only
        
        % Set a placeholder of integers the same size as "labels"
        temp = zeros(size(labels));     % Baseline levels are 0
        
        switch obj.subjects(subject)
            case 1
                
                % 1: N/A, 2: 7%
                temp(labels == Level.absDecrease7 | labels == Level.absIncrease7) = 2;
                
                % 3: N/A, 4: 14%
                temp(labels == Level.absDecrease14 | labels == Level.absIncrease14) = 4;
                
                % 5: N/A 6: 21%
                temp(labels == Level.absDecrease21 | labels == Level.absIncrease21) = 6;
                
            case 2
                
                % 1: 7%, 2: N/A
                temp(labels == Level.absDecrease7 | labels == Level.absIncrease7) = 1;
                
                % 3: 14%, 4: N/A
                temp(labels == Level.absDecrease14 | labels == Level.absIncrease14) = 3;
                
                % 5: 21%
                temp(labels == Level.absDecrease21 | labels == Level.absIncrease21) = 5;
                
                % 6: 28%
                temp(labels == Level.absDecrease28 | labels == Level.absIncrease28) = 6;
                
            case 3
                
                % 1: N/A, 2: 7%
                temp(labels == Level.absDecrease7 | labels == Level.absIncrease7) = 2;
                
                % 3: N/A, 4: 14%
                temp(labels == Level.absDecrease14 | labels == Level.absIncrease14) = 4;
                
                % 5: N/A 6: 21%
                temp(labels == Level.absDecrease21 | labels == Level.absIncrease21) = 6;
                
            case 4
                
                % 1: N/A, 2: 7%
                temp(labels == Level.absDecrease7 | labels == Level.absIncrease7) = 2;
                
                % 3: N/A, 4: 14%
                temp(labels == Level.absDecrease14 | labels == Level.absIncrease14) = 4;
                
                % 5: N/A 6: 21%
                temp(labels == Level.absDecrease21 | labels == Level.absIncrease21) = 6;
                
            case 5
                
                % 1: N/A, 2: N/A, 3: 7%
                temp(labels == Level.absDecrease7 | labels == Level.absIncrease7) = 3;
                
                % 4: N/A, 5: N/A, 6: 14%
                temp(labels == Level.absDecrease14 | labels == Level.absIncrease14) = 6;
                
            case 6
                
                % 1: 7%, 2: N/A
                temp(labels == Level.absDecrease7 | labels == Level.absIncrease7) = 1;
                
                % 3: 14%, 4: N/A
                temp(labels == Level.absDecrease14 | labels == Level.absIncrease14) = 3;
                
                % 5: 21%
                temp(labels == Level.absDecrease21 | labels == Level.absIncrease21) = 5;
                
                % 6: 28%
                temp(labels == Level.absDecrease28 | labels == Level.absIncrease28) = 6;
                
        end
        
        % Write temporary placeholder to "labels"
        labels = temp;
        
end

% Close the waitbar
close(f)

end

