function obj = confidence(obj, B)

% -------------------------------------------------------------------------
% Description:
% This function returns the confidence of predicting a class over B
% observations given the confusion matrix A.
% -------------------------------------------------------------------------

if isempty(obj.cumulativePerf)
    disp("Error in TemplateSet.confidence():")
    disp(" -> Run TemplateSet.create() before computing error")
    return
end

% Set placeholders
A = obj.cumulativePerf; n = obj.numClasses;

% Set placeholders for return values
Pe = zeros(n, B);

% Perform the analysis for each target class
for j = 1:n
    
    % Iterate over each prediction not matching the true class
    idx = 1:n; idx(j) = [];
    for i = idx
        
        % Determine the maximum value among all values in the confusion
        % matrix in the same column, besides the current element under
        % observation, and return its probability.
        probs = A(:, j); probs(i) = []; p_star = max(probs);
        
        % A unique Pe is obtained for each observation
        for b = 1:B
        for k = 1:b
            % Obtain Pe for each beat
            Pe(i, b) = Pe(i, b) + binopdf(k, b, A(i, j))*sum(binopdf(0:k-1, b, p_star));
        end
        end
        
    end
    
end

% Return final value
obj.error = Pe;

end