function [HR, RA_Pressure, Delta_P] = starling(obj, level, varargin)

% -------------------------------------------------------------------------
% This function returns mean values for parameters relevant to
% Frank-Starling curve generation: heart rate, RA pressure, and the
% pressure gradient between the aorta and left ventricle (PCWP).
%
% Arguments (required)
% - level       [Level]
%
% Arguments (optional)
% - subject             Subject for which to extract data
% - path        String  Path to co_recording_indices.mat
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'path'); path = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'subject'); subject = varargin{arg + 1};
        end
    end
end

% Set defaults for optional input arguments
if ~exist('path', 'var')
    path = '/media/Data/Hypovolemia_Study/'; addpath(path);
else; addpath(path);
end
if ~exist('subject', 'var'); subject = 1; end

% Load table with recording indices
idx = load('co_recording_indices.mat'); idx = idx.tab;

% Extract the beginning and end times
anchorTime = idx.anchor(idx.pig == obj.subjects(subject) & idx.level == level);
startTime = idx.start(idx.pig == obj.subjects(subject) & idx.level == level);
endTime = idx.finish(idx.pig == obj.subjects(subject) & idx.level == level);
set = idx.set(idx.pig == obj.subjects(subject) & idx.level == level);
if set == 1; interval = "allRelative"; else; interval = "allAbsolute"; end

% Convert times to samples
startIdx = obj.timeToSample(anchorTime, obj.Fs, startTime);
endIdx = obj.timeToSample(anchorTime, obj.Fs, endTime);

% Convert samples to beats
start = obj.sampleToBeat(startIdx, obj.dataset{subject}.(interval).beats_biopac);
finish = obj.sampleToBeat(endIdx, obj.dataset{subject}.(interval).beats_biopac);

% Todo: Adjust for subject 1
if obj.subjects(subject) == 1; start = startTime; finish = endTime; end

% Compute HR, RA pressure, and pressure gradient
HR = mean(obj.dataset{subject}.(interval).HR_biopac(start:finish));
RA_Pressure = mean(obj.dataset{subject}.(interval).rightAtrialMAP(start:finish));
aortic = obj.dataset{subject}.(interval).aorticMAP(start:finish);
pcwp = obj.dataset{subject}.(interval).wedgeMAP(start:finish);
Delta_P = mean(aortic - pcwp);

end