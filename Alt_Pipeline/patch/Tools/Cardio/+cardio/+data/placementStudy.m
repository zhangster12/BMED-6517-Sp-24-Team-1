function result = placementStudy(subjects, placements, varargin)

% -------------------------------------------------------------------------
% SUMMARY
% This function extracts data from the dataset associated with the Ashouri
% 2017 paper entitled "Automatic Detection of Seismocardiogram Sensor
% Misplacement for Robust Pre-Ejection Period Estimation in Unsupervised 
% Settings".  Data may be extracted for the specified subjects for specified
% placements, and optionally filtered with a 40Hz lowpass filter.
% 
% ARGUMENTS (REQ'D)
% - subjects        [Sx1]       List of subjects to extract (1-10)
% - placements      [Positions] List of valid positions
% 
% ARGUMENTS (OPT'L)
% - 'path'          String      Path to dataset directory
% - 'filter'        FLAG        Filter results?
% - 'int'           Interval    Interval to extract
% - 'separate'      FLAG        Beat-separation of results?
% - 'mindist'                   Minimum distance between R-peaks
% - 'ensembleLen'               Fixed length of ensemble
% - 'M'                         Smoothing factor for EMA
% - 'offset'                    Offset (in +/- samples)
% - 'verbose'       FLAG        Print updates to console?
% 
% OUTPUT
% The output is a cell array of length SUBJECTS, each containing fields
% corresponding to each of the placements as well as their respective ECG.
% e.g. result{subject}.position.ecg, result{subject}.position.scg
% -------------------------------------------------------------------------

% Parse optional arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'path'); addpath(genpath(varargin{arg + 1}));
        elseif strcmp(varargin{arg}, 'filter'); Filter = true;
        elseif strcmp(varargin{arg}, 'int'); interval = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'separate'); Separate = true;
        elseif strcmp(varargin{arg}, 'mindist'); mindist = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'ensembleLen'); ensembleLen = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'M'); M = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        elseif strcmp(varargin{arg}, 'offset'); offset = varargin{arg + 1};
        end
    end
end

% Set defaults
if ~exist('Filter', 'var'); Filter = false; end
if ~exist('interval', 'var'); interval = Intervals.All; end
if ~exist('Separate', 'var'); Separate = false; end
if ~exist('verbose', 'var'); verbose = false; end
if ~exist('mindist', 'var'); mindist = 800; end
if ~exist('ensembleLen', 'var'); ensembleLen = 800; end
if ~exist('M', 'var'); M = 0; end
if ~exist('offset', 'var'); offset = 0; end
Fs = 2000;  % Sampling frequency (property of dataset)

% Set the start and stop indices based on the selected interval
load('contents.mat', 'contents');   % Load the file
% Switch over interval for each subject to obtain start/stop indices
switch interval
    case Intervals.All
        start = contents.restBegin;
        stop = contents.recoveryEnd;
    case Intervals.Rest
        start = contents.restBegin;
        stop = contents.restEnd;
    case Intervals.Exercise
        start = contents.exerciseBegin;
        stop = contents.exerciseEnd;
    case Intervals.Recovery
        start = contents.recoveryBegin;
        stop = contents.recoveryEnd;
end

% Determine which intervals to extract
Central = ~isempty(find(placements == Positions.Center, 1));
Above = ~isempty(find(placements == Positions.Above, 1));
Below = ~isempty(find(placements == Positions.Below, 1));
Left = ~isempty(find(placements == Positions.Left, 1));
Right = ~isempty(find(placements == Positions.Right, 1));
All = ~isempty(find(placements == Positions.All, 1));

% Set placeholder for return value
result = cell(length(subjects), 1);

% If the data should be filtered, create a filter
if Filter
    Hd = cardio.general.createBPF(1, 46, Fs, 'kaiser', 'Fpass1', 2, 'Fpass2', 45);
    Hd_icg = cardio.general.createBPF(0.5, 35, Fs, 'kaiser', 'Fpass1', 1, 'Fpass2', 34);
end

c = 1;  % Initialize counter

% Print an update if indicated
if verbose; disp("Extracting Dataset"); disp("------------------"); end

% For each subject...
for s = subjects
    
    % Print an update if indicated
    if verbose; disp("Processing Subject " + string(s)); end
    
    % Set the range for data collection
    rng = start(2*s-1):stop(2*s-1);
    
    % Get the filenames for subject
    filename1 = "Subject " + string(s) + " position 1.mat";
    filename2 = "Subject " + string(s) + " position 2.mat";
    
    % Extract the left/right positions if indicated
    if Left || Right || All
        
        % Load the data
        load(filename1, 'data')
        
        % Extract the left/right/central data if indicated (and ECG)
        if Left || All; [result{c}.left.ecg, result{c}.left.scg] = extract(8, 7); end
        if Right || All; [result{c}.right.ecg, result{c}.right.scg] = extract(8, 12); end
        if Central || All; [result{c}.central_1.ecg, result{c}.central_1.scg] = extract(8, 4); end
        
        % Extract ICG
        [~, result{c}.icg_1] = extract(8, 11);
        
    end
    
    % Extract the upper/lower positions if indicated
    if Above || Below || All
        
        % Set the range for data collection
        rng = start(2*s):stop(2*s);
        
        % Load the data
        load(filename2, 'data');
        
        % Correct for the quirk with Subject 9
        if s ~= 9
        
            % Extract the upper/lower/central data if indicated (and ECG)
            if Above || All; [result{c}.above.ecg, result{c}.above.scg] = extract(8, 7); end
            if Below || All; [result{c}.below.ecg, result{c}.below.scg] = extract(8, 12); end
            if Central || All; [result{c}.central_2.ecg, result{c}.central_2.scg] = extract(8, 4); end
            
        else
            
            % Extract the upper/lower/central data if indicated (and ECG)
            if Above || All; [result{c}.above.ecg, result{c}.above.scg] = extract(8, 7); end
            if Below || All; [result{c}.below.ecg, result{c}.below.scg] = extract(8, 4); end
            if Central || All; [result{c}.central_2.ecg, result{c}.central_2.scg] = extract(8, 12); end
            
        end
        
    end
    
    % Extract the central position if it hasn't been already
    if Central && ~Left && ~Right && ~All
        load(filename1, 'data'); rng = start(2*s-1):stop(2*s-1);
        [result{c}.central_1.ecg, result{c}.central_1.scg] = extract(8, 4);
    end
    if Central && ~Above && ~Below && ~All
        load(filename2, 'data'); rng = start(2*s):stop(2*s);
        [result{c}.central_2.ecg, result{c}.central_2.scg] = extract(8, 4);
    end
     
    % Extract ICG
    [~, result{c}.icg_2] = extract(8, 11);
    
    c = c + 1;  % Increment counter
    
end

% Print an update if indicated
if verbose; disp("Extraction Complete"); end

% Nested function for separating beats
function [ecg, scg] = separateBeats(ecg, scg)

    % Threshold the heart rate
    [~, ~, threshold] = cardio.ecg.automatedHR(ecg, 100, Fs, mindist, false);
    % Beat-separate the data
    [ecg, scg] = cardio.general.separateBeats(scg, 'ecg', ecg, ...
        'threshold', threshold, 'minDist', mindist, 'samples', ensembleLen, 'offset', offset);

end

% Nested function for extracting data
function [ecg, scg] = extract(ecgCol, dataCol)
    
    % Extract data
    ecg = data(rng, ecgCol);    % Set ECG
    if Filter && dataCol == 11; scg = filtfilt(Hd_icg.Numerator, 1, data(rng, dataCol));
    elseif Filter; scg = filtfilt(Hd.Numerator, 1, data(rng, dataCol)); 
    else; scg = data(rng, dataCol);
    end
    % Separate beats if indicated
    if Separate
        [ecg, scg] = separateBeats(ecg, scg);
        % Smooth beats if indicated
        if M > 0; scg = cardio.general.ema(scg, M, false); end
    end
    
end

end
