function equalized = equalize(data, varargin)

% -------------------------------------------------------------------------
% This function equalizes the length of variable-length vectors via
% truncation or re-sampling. The truncated vectors of non-zero length are
% returned in a column matrix.
%
% Arguments (required)
% - data        {Nx1}	Cell vector of signal vectors to equalize
%
% Arguments (optional)
% - 'resample'  FLAG    Use re-sampling method rather than truncation
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'resample'); truncate = false; end
    end
end

% Set defaults for optional arguments
if ~exist('truncate', 'var'); truncate = true; end

% Initialize minimum length and initialize counter
minLen = inf; validIdx = [];

% For each vector in the cell vector
for i = 1:length(data)
    
    % Extract the vector
    vector = data{i};
    
    % If the length of the vector is greater than 0, continue
    if length(vector) < 1; continue; end
    
    % Is the vector less than the minimum length? Update minLen.
    if length(vector) < minLen; minLen = length(vector); end
    
    % Increment the counter for valid vectors
    validIdx = [validIdx i];
    
end

% Initialize placeholder for return value
equalized = zeros(minLen, length(validIdx));

% For each valid vector...
for i = 1:length(validIdx)
    
    % Extract the vector
    vector = data{validIdx(i)};
    
    % Truncate the vector, if indicated. Else, resample.
    if truncate
        vector = vector(1:minLen);
    else
        vector = resample(vector, minLen, length(vector));
    end
    
    % Append the vector to the return matrix
    equalized(:, i) = vector(:); 
    
end

end

