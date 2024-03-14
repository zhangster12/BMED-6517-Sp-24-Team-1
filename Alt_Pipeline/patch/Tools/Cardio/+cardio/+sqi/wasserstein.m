function [transformed, varargout] = wasserstein(signal, template, varargin)

% -------------------------------------------------------------------------
% This function computes the Wasserstein distance between a signal and
% template. The function optionally returns the transformed template and
% the transformation for converting the template to signal vector.
% (This method is also called "earth mover distance")
%
% This method is adapted from
% https://vincentherrmann.github.io/blog/wasserstein/
%
% Arguments (required)
% - signal      [Nx1]   Signal vector
% - template    [Nx1]   Template vector
%
% Arguments (optional)
% - 'gamma'     [NxN]   Transformation matrix (provided)
% - 'plot'      FLAG    Plot results?
% - 'verbose    FLAG    Print progress to command window?
%
% Outputs:
% - transformed     [Nx1]   Transformed template vector
% - varargout{1}            Distance
% - varargout{2}    [NxN]   Transformation matrix
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'gamma'); gamma = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'plot'); Plot = true;
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        end
    end
end

% Set defaults for optional arguments
if ~exist('Plot', 'var'); Plot = false; end
if ~exist('verbose', 'var'); verbose = false; end

% Ensure that both signals have same integral (1)
% signal = normalize(signal); template = normalize(template);
signal = signal + abs(min(signal)); template = template + abs(min(template));
signal = signal./sum(signal); template = template./sum(template);

% Get length of signal
sigLen = length(signal);

% Create matrix A and initialize counter
A = zeros(sigLen^2, 2*sigLen); i = 1;
for counter = 1:sigLen
    A1 = zeros(sigLen); A1(:, counter) = ones(sigLen, 1);
    A2 = eye(sigLen);
    A(i:(i + sigLen - 1), :) = [A1 A2]; i = i + sigLen;
end

% Create matrix D and initialize counter/loop flag
D = zeros(sigLen); flag = false; counter = 1;
while ~flag
    try
        D = D + counter*diag(ones(sigLen - counter, 1), counter);
        D = D + counter*diag(ones(sigLen - counter, 1), -counter);
        counter = counter + 1;
    catch
        flag = true;
    end
end

% Create vector b
b = [template(:); signal(:)];

% Create vector c
c = D(:);

% Find gamma, if necessary
if ~verbose && ~exist('gamma', 'var')
    options = optimoptions('linprog', 'Display', 'none');
    gamma = linprog(c, [], [], A', b, zeros(sigLen^2, 1), [], options);
elseif ~exist('gamma', 'var')
    gamma = linprog(c, [], [], A', b, zeros(sigLen^2, 1), []);
end

% Contstruct transformation
transformation = zeros(sigLen);
for i = 1:sigLen
    transformation(i, :) = gamma(((i-1)*sigLen + 1):((i-1)*sigLen + sigLen));
end

% Apply transformation to template
transformed = template;
for i = 1:sigLen
    for j = 1:sigLen
        transformed(i) = transformed(i) - transformation(i, j);
        transformed(j) = transformed(j) + transformation(i, j);
    end
end

% Plot results, if indicated
if Plot
    figure; subplot(1, 2, 1); imagesc(transformation > 0); title("Transformation")
    xlabel("To"); ylabel("From"); subplot(1, 2, 2); imagesc(D); title("Distance")
    xlabel("To"); ylabel("From");
end

% Return output arguments
if nargout > 1; varargout{1} = sum(transformation.*D, 'all'); end   % EMD
if nargout > 2; varargout{2} = transformation; end

end

