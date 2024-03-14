function R2 = rsquared(signal1, signal2, varargin)

% -------------------------------------------------------------------------
% This function returns the coefficient of determination between two signals,
% handling NaN values if necessary. The function also aligns the signals if
% indices have been removed. Input arguments include:
% - signal1     [Nx1]   Signal Vector 1
% - signal2     [Nx1]   Signal Vector 2
% - varargin
%   - ri1       [Mx1]   Vector of indices from signal1 that have been removed
%   - ri2       [Lx1]   Vector of indices from signal2 that have been removed
% -------------------------------------------------------------------------

% Parse inputs
if nargin == 4; ri1 = varargin{1}; ri2 = varargin{2}; end

% Align vectors to adjust for removed indices
% Arguments: vector1, vector2, remIdx1, remIdx2
if exist('ri1', 'var') && exist('ri2', 'var')
    [aligned1, aligned2] = cardio.general.alignVectors(signal1, signal2, ri1, ri2);
else
    aligned1 = signal1; aligned2 = signal2;
end

% Compute covariance
% cov_a = nancov([aligned1 aligned2]);

% Compute correlation coefficent
% c = (cov_a(1,2)/(nanstd(signal1)*nanstd(signal2)))^2;

% Generate a linear regression fitting signal 1 to signal2
mdl = fitlm(aligned1, aligned2);

% Return the R^2 value (coefficient of determination)
R2 = mdl.Rsquared.Adjusted;

end