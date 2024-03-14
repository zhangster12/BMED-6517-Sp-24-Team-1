function data = getSubinterval(obj, fields, levels, varargin)

% -------------------------------------------------------------------------
% This function pulls data from a sub-interval of either absolute or
% relative hypovolemia test intervals.
%
% Arguments (required)
% - fields      [String]    Fields for which to segment data
% - levels      [Level]     List of levels for which to extract data
%
% Arguments (optional)
% - subjects                List of subjects for which to extract data
% - under       String      Struct name under parent level
% - of          [NxM]       Specify M data vectors of length N to process 
%                           (one subject/level only)
% - t3          FLAG        Use T3 system for segmentation (for
%                           non-beat-separated signals)
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        try
            switch varargin{arg}
                case 'subjects'; subjects = varargin{arg + 1};
                case 'under'; under = varargin{arg + 1};
                case 'of'; givenData = varargin{arg + 1};
                case 't3'; t3 = true;
            end
        catch
            continue
        end
    end
end

% Set default values for optional arguments
if ~exist('subjects', 'var'); subjects = obj.subjects; end
if ~exist('t3', 'var'); t3 = false ;end

% -------------------------------------------------------------------------
% Provided Data
% -------------------------------------------------------------------------

% If the data vector is provided, segment the vector
if exist('givenData', 'var')
    
    % Correct the subject number
    if length(subjects) > 1; subjects = subjects(1); end
    subject = find(obj.subjects == subjects);
    
    % Ensure that only one level provided
    if length(levels) > 1; levels = levels(1); end; level = levels;
    
    % Determine the parent for the level
    parent = obj.getParent(level);
    
    % Extract beats from the Biopac data
    beats = obj.dataset{subject}.(string(level)).indices;
    beats = beats - obj.dataset{subject}.(string(parent)).indices(1) + 1;
    
    % Extract the feature vector for the selected beats
    if size(givenData, 2) > 1
        data = givenData(beats(1):beats(end), :);
    else
        data = givenData(beats(1):beats(end));
    end
    
    return
    
end

% Ensure "fields" is a string
if isa(fields, 'char'); fields = convertCharsToStrings(fields); end

% -------------------------------------------------------------------------
% Non-Segmented Data
% -------------------------------------------------------------------------

if ~obj.beatSeparated
    
    % For each subject...
    for subject = 1:length(subjects)
    
        % For each level...
        for level = levels
            
            % Determine the parent for the level
            parent = obj.getParent(level);
            
            % Get the beginning/ending indices
            if ~t3
                % Use the Biopac reference for segmentation
                idx = obj.dataset{subject}.(string(level)).samples_b;
                idx = idx - obj.dataset{subject}.(string(parent)).samples_b(1) + 1;
            else
                % Use the T3 reference for segmentation
                idx = obj.dataset{subject}.(string(level)).samples_t;
                idx = idx - obj.dataset{subject}.(string(parent)).samples_t(1) + 1;
            end

            % For each field...
            for field = fields

                % Extract parent feature vector for the selected beats
                if ~exist('under', 'var')
                    allData = obj.dataset{subject}.(string(parent)).(field);
                else
                    allData = obj.dataset{subject}.(string(parent)).(under).(field);
                end
                
                % If there is more than one field to be specified, return
                % the data as a struct
                if length(fields) > 1; data.(field) = allData(idx(1):idx(2));
                else; data = allData(idx(1):idx(2));
                end

            end

        end
    
    end
    
    % Return when complete
    return
    
end

% -------------------------------------------------------------------------
% Segmented Data
% -------------------------------------------------------------------------

% For each subject...
for subject = 1:length(subjects)
    
    % For each level...
    for level = levels
        
        % Determine the parent for the level
        parent = obj.getParent(level);
        
        % For each field...
        for field = fields
            
            % Extract beats from the Biopac data
            beats = obj.dataset{subject}.(string(level)).indices;
            beats = beats - obj.dataset{subject}.(string(parent)).indices(1) + 1;

            % Extract feature vector for the selected beats
            if ~exist('under', 'var')
                allData = obj.dataset{subject}.(string(parent)).(field);
            else
                allData = obj.dataset{subject}.(string(parent)).(under).(field);
            end
            if length(fields) > 1
                if size(allData, 2) > 1
                    data.(field) = allData(:, beats(1):beats(end));
                else
                    data.(field) = allData(beats(1):beats(end));
                end
            else
                if size(allData, 2) > 1
                    if beats(end) <= length(allData)
                        data = allData(:, beats(1):beats(end));
                    else
                        data = allData(:, beats(1):end);
                    end
                else
                    if beats(end) <= length(allData)
                        data = allData(beats(1):beats(end));
                    else
                        data = allData(beats(1):end);
                    end
                end
            end
        
        end
        
    end
    
end

end

