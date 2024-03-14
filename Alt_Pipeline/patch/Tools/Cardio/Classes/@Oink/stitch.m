function stitched = stitch(obj, field, level, varargin)

% -------------------------------------------------------------------------
% This function stitches together beat-separated data to form the original
% continuous waveform.
%
% Arguments (required)
% - field       String      Field name for signal to stitch
% - level       Level       Level for which to stitch signal
%
% Arguments (optional)
% - subject                 Subject for which to stitch signal
% - under       String      Parent field name if the field is a subfield
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'under'); under = varargin{arg + 1}; 
        elseif strcmp(varargin{arg}, 'subject'); subject = varargin{arg + 1};
        end
    end
end

% Set defaults for optional input arguments
if ~exist('under', 'var'); under = []; end
if ~exist('subject', 'var'); subject = 1; end

% Extract beat-separated data
if isempty(under)
    orig = obj.getSubinterval(field, level, 'subjects', subject);
else
    orig = obj.getSubinterval(field, level, 'subjects', subject, 'under', under);
end

% Set placeholder for stitched data
stitched = [];

% For each interval in the original data...
for i = 1:size(orig, 2)
    % Append the vector to the stitched data
    stitched = [stitched; orig(:, i)];
end

end

