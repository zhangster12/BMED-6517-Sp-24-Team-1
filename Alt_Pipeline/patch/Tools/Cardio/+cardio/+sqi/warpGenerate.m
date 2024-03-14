function [new_signal, varargout] = warpGenerate(signal, varargin)

% -------------------------------------------------------------------------
% This function generates a new signal from a template by generating a
% random DTW warp path and deriving the new signal from the warp path. In
% this manner, the old and new signals are DTW-equivalents.
%
% Arguments (req'd)
% - signal  [Nx1]   Template signal
%
% Arguments (opt'l)
% - prob    (0, 1]  Probability of changing warp path at each step
%                   (1 = random path, 0 = fixed path in random direction)
%                   (default 0.5)
%
% Outputs (opt'l)
% - varargout{1}    Warp grid for generating new signal
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        switch varargin{arg}
            case 'prob'; prob = varargin{arg + 1};
        end
    end
end

% Set defaults for optional arguments
if ~exist('prob', 'var'); prob = 0.5; end

% Return an error if the probability is 0
if prob == 0; disp("Error in warpGenerate.m: prob must be > 0"); return; end

% Get the length of the signal
len = length(signal);

% Initialize path
path = [1; 1];

% 1 = row; 2 = col; 3 = diag;
% Set starting momentum
momentum = int64(randi(3, 1));

% Set loop FLAG
FLAG = true;

% Create path
while FLAG
    
    % If keeping the same momentum, just update the path and continue
    if rand(1) > prob; path = updatePath(); continue; end
    
    % Else, update the momentum before updating the path
    momentum = int64(randi(3, 1)); path = updatePath();
    
    % Stopping criteria
    % 1. If the bottom corner of the warp grid has been reached, stop
    % 2. If the maximum length of the path has been reached, stop
    if isequal(path(:, end), [len; len]); FLAG = false; end
    if size(path, 2) > 2*len; FLAG = false; end

end

% Generate warp grid
grid = zeros(len);
for i = 1:size(path, 2); grid(path(1, i), path(2, i)) = 1; end

% Get new path length
newLen = size(path, 2);

% Initialize and generate t_prime (warped template)
t_prime = zeros(newLen, 1);
for i = 1:newLen; t_prime(i) = signal(path(1, i)); end

% Reverse the generation to get the new signal
new_signal = zeros(len, 1);     % Initialize the new signal
for i = 1:len
    % For horizontal paths on the warp grid, select the grid point farthest
    % to the right (such that the lower-right grid point will always be
    % reached). Otherwise just map the point on the new signal directly to
    % the corresponding point on the old signal.
    idx = find(path(1, :) == i, 1, 'last');
    new_signal(i) = signal(path(2, idx));
end

% Parse output arguments
if nargout > 1; varargout{1} = grid; end

    % Helper function for updating the warp path
    function newPath = updatePath()
        
        % Is the momentum valid?
        % 1. If you're in the bottom row, you must proceed row-wise
        % 2. If you're in the rightmost column, you must proceed columnwise
        if path(1, end) == len && momentum ~= 1; newPath = path; return; end
        if path(2, end) == len && momentum ~= 2; newPath = path; return; end
        
        % Update the path based on the momentum
        switch momentum
            case 1  % Row-wise momentum
                newPath = [path path(:, end) + [0; 1]];
            case 2  % Columnwise momentum
                newPath = [path path(:, end) + [1; 0]];
            case 3  % Diagonal momentum
                newPath = [path path(:, end) + [1; 1]];
        end
    end
    
end