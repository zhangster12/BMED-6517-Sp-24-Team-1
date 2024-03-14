function [coeffs, score] = regPCA(data, varargin)

% -------------------------------------------------------------------------
% This function implements regularized principal component analysis (PCA)
% to find uniform subspace trajectories. This is accomplished by adding a
% distance penalty to PCA for trajectory distance.
%
% ARGUMENTS (REQ'D)
% - data    [MxN]   M row-wise observations of N features
%
% ARGUMENTS (OPT'L)
% - step            Step size for gradient descent
% - lambda          Penalty for distance between subspace trajectories
% - indices [Lx1]   Indices at which the data matrix is separated by trial
% - penalties       Number of PCs on which to impose penalty
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Part 1: Setup
% -------------------------------------------------------------------------

% Parse variable inputs
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'step'); step = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'lambda'); lambda = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'indices'); indices = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'penalties'); penalties = varargin{arg + 1};
        end
    end
end

% Set default variable values
if ~exist('step', 'var'); step = 0.001; end
if ~exist('lambda', 'var'); lambda = 0; end
if ~exist('indices', 'var'); indices = 1; end
if ~exist('penalties', 'var'); penalties = 0; end

% Get number of features
numFeatures = size(data, 2);

% Mean-center the original dataset
% Set a placeholder for the matrix X (which will be modified)
orig_data = data; data = data - mean(data); X = data;

% Initialize a placeholder for the coefficients
coeffs = rand(numFeatures, numFeatures);    % Initialize random coeffs
coeffs = coeffs./vecnorm(coeffs);           % Normalize columns

% -------------------------------------------------------------------------
% Part 2: Compute Coefficients
% -------------------------------------------------------------------------

% A vector of coefficients is needed for each score vector.
for i = 1:numFeatures
    
    FLAG = true; % Flag for continuing loop
    
    % Extract the current coefficient vector
    current_coeff = coeffs(:, i);
    
    % Update the data matrix based on previous coefficient vectors
    if i > 1; X = X - X*coeffs(:, i-1)*coeffs(:, i-1)'; end
    
    % Perform gradient descent to find coefficient vector
    while FLAG
        
        % Do not repeat loop unless otherwise indicated
        FLAG = false;

        % Update each element in coefficient vector
        for element = 1:numFeatures

            % Record the initial loss
            init_loss = computeLoss(current_coeff);

            % Create temporary coeff vectors with the selected elements
            % modified in different ways.
            coeff_max = current_coeff; coeff_max(element) = coeff_max(element) + step;
            coeff_min = current_coeff; coeff_min(element) = coeff_min(element) - step;
            coeff_max = coeff_max/norm(coeff_max); coeff_min = coeff_min/norm(coeff_min);

            % Get loss from modified coefficient vectors
            loss_max = computeLoss(coeff_max);
            loss_min = computeLoss(coeff_min);
            
            % Select updated coefficients based on loss maximization
            if loss_max > init_loss || loss_min > init_loss
                if loss_max > loss_min; current_coeff = coeff_max; 
                else; current_coeff = coeff_min;
                end; FLAG = true;   % If an update was made, continue the loop
            end

        end
        
    end
    
    % Update coefficients
    coeffs(:, i) = current_coeff;
    
end

% -------------------------------------------------------------------------
% Part 3: Return Scores
% -------------------------------------------------------------------------

score = (orig_data - mean(orig_data))*coeffs;

% -------------------------------------------------------------------------
% APPENDIX
% -------------------------------------------------------------------------

    function loss = computeLoss(coefficients)
        
        % -----------------------------------------------------------------
        % First Loss Term
        % -----------------------------------------------------------------
        
        % The first term of the loss function is the norm-squared error
        term1 = norm(X*coefficients)^2;
        
        % -----------------------------------------------------------------
        % Second Loss Term
        % -----------------------------------------------------------------
        
        % Determine whether the second loss term is indicated
        if i <= penalties
            % The second term is the norm-squared error of specified points
            % First, compute the scores for each specified point
            penaltyScores = X(indices, :)*coefficients;
            % Then, compute the second term
            term2 = norm(penaltyScores - mean(penaltyScores))^2;
        else; term2 = 0; % If it is not indicated, return 0
        end
        
        % Return the loss
        loss = term1 - lambda*term2;
        
    end

end
