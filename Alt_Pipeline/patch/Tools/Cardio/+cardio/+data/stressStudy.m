function [indices, scg, icg] = stressStudy(subjects, varargin)

% -------------------------------------------------------------------------
% This function extracts data collected from the study Fusing Near-Infrared
% Spectroscopy with Wearable Hemodynamic Measurements Improves
% Classification of Mental Stress by Gurel et al (2018). 
%
% ARGS (REQ'D)
% - subjects [Sx1]  Subjects to include (102-119)
%
% ARGS (OPT'L)
% - 'filter'    FLAG    Filter the data?
% - 'separate'  FLAG    Beat-separate the data?
% - 'path'      String  Path/to/dataset (if not on server)
% - 'verbose'   FLAG    Print updates?
% - 'M'         Int     Smoothing factor for EMA
% - 'offset'    Int     Offset (in +/- samples) of scg and icg
% -------------------------------------------------------------------------

% Data file organization:
% 1 - RSP
% 2 - ICG_Z
% 3 - PPG
% 4 - EDA
% 5 - ECG
% 6 - AX
% 7 - ICG_DZDT
% 8 - AY
% 9 - AZ
% 10 - Custom PPG
        
% Parse input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'filter'); Filter = true;
        elseif strcmp(varargin{arg}, 'separate'); Separate = true;
        elseif strcmp(varargin{arg}, 'path'); Path = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'verbose'); Verbose = true;
        elseif strcmp(varargin{arg}, 'M'); M = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'offset'); offset = varargin{arg + 1};
        end
    end
end

% Set default arguments
if ~exist('Filter', 'var'); Filter = false; end
if ~exist('Separate', 'var'); Separate = false; end
if ~exist('Path', 'var'); Path = "/media/Data/Stress_Study/"; end
if ~exist('Verbose', 'var'); Verbose = false; end
if ~exist('M', 'var'); M = 1; end
if ~exist('offset', 'var'); offset = 0; end
mindist = 800; ensembleLen = 800; Fs = 2000;

% Create placeholders for the raw SCG data and indices / set counter
indices = zeros(length(subjects)+1, 1);
raw_scg_x = cell(length(subjects), 1);
raw_scg_y = cell(length(subjects), 1);
raw_scg_z = cell(length(subjects), 1);
raw_icg = cell(length(subjects), 1); c = 1;

% Generate filter if needed
if Filter
    Hd_scg = cardio.general.createBPF(1, 46, Fs, 'kaiser', 'Fpass1', 2, 'Fpass2', 45);
    Hd_icg = cardio.general.createBPF(1, 31, Fs, 'kaiser', 'Fpass1', 2, 'Fpass2', 30);
end

% Print updates (if indicated)
if Verbose; disp("Extracting Data"); end

% For each subject
for s = subjects
    
    % Print updates (if indicated)
    if Verbose
        disp(" -> Subject " + string(c) + " of " + string(length(subjects)));
    end
    
    data = importdata(Path + string(s) + ".mat");
    data = data.data; raw_scg_x{c} = data(:, 6);
    raw_scg_y{c} = data(:, 8); raw_scg_z{c} = data(:, 9);
    raw_icg{c} = data(:, 7); raw_ecgs = data(:, 5);
    
    % Filter raw data
    if Filter
        raw_scg_x{c} = filtfilt(Hd_scg.Numerator, 1, raw_scg_x{c});
        raw_scg_y{c} = filtfilt(Hd_scg.Numerator, 1, raw_scg_y{c});
        raw_scg_z{c} = filtfilt(Hd_scg.Numerator, 1, raw_scg_z{c});
        raw_icg{c} = filtfilt(Hd_icg.Numerator, 1, raw_icg{c});
    end
    
    % Beat-separate raw data
    % Threshold the heart rate
    [~, ~, threshold] = cardio.ecg.automatedHR(raw_ecgs, 100, Fs, mindist, false);
    % Beat-separate the data
    [~, raw_scg_x{c}] = cardio.general.separateBeats(raw_scg_x{c}, 'ecg', ...
        raw_ecgs, 'threshold', threshold, 'minDist', mindist, 'samples', ensembleLen, 'offset', offset);
    [~, raw_scg_y{c}] = cardio.general.separateBeats(raw_scg_y{c}, 'ecg', ...
        raw_ecgs, 'threshold', threshold, 'minDist', mindist, 'samples', ensembleLen, 'offset', offset);
    [~, raw_scg_z{c}] = cardio.general.separateBeats(raw_scg_z{c}, 'ecg', ...
        raw_ecgs, 'threshold', threshold, 'minDist', mindist, 'samples', ensembleLen, 'offset', offset);
    [~, raw_icg{c}] = cardio.general.separateBeats(raw_icg{c}, 'ecg', ...
        raw_ecgs, 'threshold', threshold, 'minDist', mindist, 'samples', ensembleLen, 'offset', offset);
    
    % Smooth the data if indicated
    if M > 1; raw_scg_x{c} = cardio.general.ema(raw_scg_x{c}, M, false); end
    if M > 1; raw_scg_y{c} = cardio.general.ema(raw_scg_y{c}, M, false); end
    if M > 1; raw_scg_z{c} = cardio.general.ema(raw_scg_z{c}, M, false); end
    if M > 1; raw_icg{c} = cardio.general.ema(raw_icg{c}, M, false); end
    
    % Record indices
    c = c + 1; indices(c) = indices(c - 1) + size(raw_scg_z{c-1}, 2);
    
end

% Combine the SCG data
scg = cell(length(subjects), 1);
for i = 1:length(subjects)
    scg{i} = cat(3, raw_scg_x{i}, raw_scg_y{i}, raw_scg_z{i});
end

% Return data
icg = raw_icg; indices(1) = 1;

end
