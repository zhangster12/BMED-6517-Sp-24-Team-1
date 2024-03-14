function [ecg, acc, gyr_x, icg, varargout] = exerciseStudy(subjects, interval, varargin)

% -------------------------------------------------------------------------
% SUMMARY
% This function extracts data from the exercise study dataset. Data is
% extracted from the interval specified and optionally beat-segmented.
% 
% ARGUMENTS (REQ'D)
% - subjects  [1xN]   List of subjects to process (102-117)
% - interval          Interval to extract ('rest', 'exercise', or 'recovery)
%     - 'squat' or 'walk' for exercise sub-intervals
% 
% ARGUMENTS (OPT'L)
% - 'path'                Path to dataset
% - 'segment' FLAG        Return beat-segmented data?
% - 'minDist'             Minimum distance (in samples) between beats
% - 'ensembleLen'         Length of beats for segmentation
% - 'M'                   Coefficient for exponential moving average
% - 'verbose' FLAG        Print updates?
% - 'train'   FLAG        Training dataset (subjects 102-129)
% - 'test'    FLAG        Testing dataset (subjects 201-210)
% - 'offset'              Offset (in +/- samples) of acc, gyr_x, and icg
% 
% OUTPUTS
% - indices   [(N+1)x1]	Indices for each subject
% -------------------------------------------------------------------------

% Parse optional arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'path'); path = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'segment'); segment = true;
        elseif strcmp(varargin{arg}, 'minDist'); minDist = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'ensembleLen'); ensembleLen = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'M'); M = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        elseif strcmp(varargin{arg}, 'train'); Train = true;
        elseif strcmp(varargin{arg}, 'test'); Train = false;
        elseif strcmp(varargin{arg}, 'offset'); offset = varargin{arg + 1};
        end
    end
end

% Set defaults
if ~exist('segment', 'var'); segment = false; end
if exist('path', 'var'); addpath(path); end
if ~exist('minDist', 'var'); minDist = 800; end
if ~exist('ensembleLen', 'var'); ensembleLen = 800; end
if ~exist('verbose', 'var'); verbose = false; end
if ~exist('Train', 'var'); Train = true; end
if ~exist('offset', 'var'); offset = 0; end

% Set sampling frequency
Fs = 2000;

% Load project data
if Train; info = readtable('subject_info_training.csv');
else; info = readtable('subject_info_testing.csv');
end

% Set placeholders
data = cell(length(subjects), 1);           % Subject data
ecg = cell(size(data));                     % ECG data
acc = cell(size(data));                     % Combined acceleration data
acc_x = cell(size(data));                   % X-acceleration data
acc_y = cell(size(data));                   % Y-acceleration data
acc_z = cell(size(data));                   % Z-acceleration data
gyr_x = cell(size(data));                   % X-axis rotation data
icg = cell(size(data));                     % ICG data
indices = zeros(length(subjects) + 1, 1);   % Indices for each subject

% Print update (if indicated)
if verbose; disp("Extracting Data"); end

% Load subject data from files
for s = 1:length(subjects)
    
    % Print update (if indicated)
    if verbose; disp(" -> Processing " + string(s) + " of " + string(length(subjects))); end
    
    % Set filename and import data
    filename = 'Filtered_Subject_' + string(subjects(s)) + ...
        '_Mid_Sternum_Rest_Exer_Rec.mat';
    data = importdata(filename);
    
    %{
    Table of Contents:
    data.t                  time vector                 {1}
    data.ax_f               x acceleration (sternum)    {2}
    data.ay_f               y acceleration (sternum)    {3}
    data.az_f               z acceleration (sternum)    {4}
    data.bcg_filtered       BCG (sternum)               {5}
    data.bp_filtered        blood pressure              {6}
    data.ecg_filtered       ECG 1                       {7}
    data.ecg_soft_filtered  ECG 2                       {8}
    data.gyro_x_filtered    pitch (sternum)             {9}
    data.gyro_y_filtered    roll (sternum)              {10}
    data.gyro_z_filtered    yaw (sternum)               {11}
    data.icg_filtered       ICG                         {12}
    %}
    
    % Get valid indices for analysis
    % Create subtable for subject
    rows = (info.SUBJECTID == subjects(s));
    
    if Train
        
        % Extract data from table
        vars = {'RESTSTART','RESTEND','RECOVERYSTART','RECOVERYEND', 'EXERCISESTART', 'EXERCISEEND', 'SQUATSTART', 'SQUATEND'};
        subtable = info(rows,vars);
        
        % Split dataset based on selected interval
        subtable{1, :} = subtable{1, :}.*Fs;    % Convert seconds to samples
        if interval == Intervals.Rest; idx = subtable.RESTSTART+1:subtable.RESTEND;
        elseif interval == Intervals.Recovery; idx = subtable.RECOVERYSTART:subtable.RECOVERYEND;
        elseif interval == Intervals.Exercise; idx = subtable.EXERCISESTART:subtable.SQUATEND;
        elseif interval == Intervals.Walk; idx = subtable.EXERCISESTART:subtable.EXERCISEEND;
        elseif interval == Intervals.Squat; idx = subtable.SQUATSTART:subtable.SQUATEND;
        end
        
        % Assign output data
        ecg{s} = data.ecg_filtered(idx); acc_z{s} = data.az_f(idx);
        acc_x{s} = data.ax_f(idx); acc_y{s} = data.ay_f(idx);
        gyr_x{s} = data.gyro_x_filtered(idx); icg{s} = data.icg_filtered(idx);
        
    else
        
        % Extract data from table
        vars = {'REST1START','REST1END','VALSALVASTART','VALSALVAEND', 'REST2START', 'REST2END', 'EXERCISESTART', 'EXERCISEEND'...
            'SQUATSSTART', 'SQUATSEND', 'REST3START', 'REST3END', 'COLDPRESSSTART', 'COLDPRESSEND', 'REST4START', 'REST4END'};
        subtable = info(rows,vars);
        
        % Split dataset based on selected interval
        subtable{1, :} = subtable{1, :}.*Fs;    % Convert seconds to samples
        if interval == Intervals.Rest; idx = subtable.REST1START+1:subtable.REST1END;
        elseif interval == Intervals.Valsalva; idx = subtable.VALSALVASTART:subtable.VALSALVAEND;
        elseif interval == Intervals.ValsalvaRecovery; idx = subtable.REST2START:subtable.REST2END;
        elseif interval == Intervals.Exercise; idx = subtable.EXERCISESTART:subtable.EXERCISEEND;
        elseif interval == Intervals.Squat; idx = subtable.SQUATSSTART:subtable.SQUATSEND;
        elseif interval == Intervals.ExerciseRecovery || interval == Intervals.Recovery
            idx = subtable.REST3START:subtable.REST3END;
        elseif interval == Intervals.ColdPress; idx = subtable.COLDPRESSSTART:subtable.COLDPRESSEND;
        elseif interval == Intervals.ColdPressRecovery; idx = subtable.REST4START:subtable.REST4END;
        end
        
        % Assign output data
        ecg{s} = data.ecg_filtered(idx); acc_z{s} = data.az_f(idx);
        acc_x{s} = data.ax_f(idx); acc_y{s} = data.ay_f(idx);
        gyr_x{s} = data.gyro_x_filtered(idx); icg{s} = data.icg_filtered(idx);
        
    end
    
    % Segment the dataset if indicated
    if segment

        % Compute heart rate
        % Arguments: ecg, numThresholds, Fs, mindist, Plot
        [~, ~, threshold] = cardio.ecg.automatedHR(ecg{s}, 100, Fs, minDist, false);
        
        % Separate SCG by beats
        % Arguments: FilteredECG, Threshold, MinDist, Signal, Samples, SignalName, plotECG, Plot
        [ecg_slice, ax_slice] = cardio.general.separateBeats(acc_x{s}, 'ecg', ...
            ecg{s}, 'threshold', threshold, 'minDist', minDist, 'samples', ensembleLen, 'offset', offset);
        [~, ay_slice] = cardio.general.separateBeats(acc_y{s}, 'ecg', ...
            ecg{s}, 'threshold', threshold, 'minDist', minDist, 'samples', ensembleLen, 'offset', offset);
        [~, az_slice] = cardio.general.separateBeats(acc_z{s}, 'ecg', ...
            ecg{s}, 'threshold', threshold, 'minDist', minDist, 'samples', ensembleLen, 'offset', offset);
        [~, gx_slice] = cardio.general.separateBeats(gyr_x{s}, 'ecg', ...
            ecg{s}, 'threshold', threshold, 'minDist', minDist, 'samples', ensembleLen, 'offset', offset);
        [~, icg_slice] = cardio.general.separateBeats(icg{s}, 'ecg', ...
            ecg{s}, 'threshold', threshold, 'minDist', minDist, 'samples', ensembleLen, 'offset', offset);
        
        % Assign slices to variables
        ecg{s} = ecg_slice; acc_z{s} = az_slice; gyr_x{s} = gx_slice; 
        icg{s} = icg_slice; acc_x{s} = ax_slice; acc_y{s} = ay_slice;
        
        % Obtain EMA if indicated
        if exist('M', 'var')
            acc_x{s} = cardio.general.ema(acc_x{s}, M, false);  % X-Acceleration
            acc_y{s} = cardio.general.ema(acc_y{s}, M, false);  % Y-Acceleration
            acc_z{s} = cardio.general.ema(acc_z{s}, M, false);	% Z-Acceleration
            gyr_x{s} = cardio.general.ema(gyr_x{s}, M, false);	% X-Gyro
            icg{s} = cardio.general.ema(icg{s}, M, false);      % ICG
        end

    end
    
    % Record indices for subject
    indices(s + 1) = indices(s) + size(acc_z{s}, 2);
    
end

% Combine accelerometer dimensions
for s = 1:length(subjects); acc{s} = cat(3, acc_x{s}, acc_y{s}, acc_z{s}); end

% Provide optional outputs
indices(1) = 1; varargout{1} = indices;

end

