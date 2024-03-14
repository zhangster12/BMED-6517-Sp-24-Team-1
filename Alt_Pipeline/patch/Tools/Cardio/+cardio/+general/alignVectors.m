function [alignedV1, alignedV2] = alignVectors(vector1, vector2, ...
    remIdx1, remIdx2)

% -------------------------------------------------------------------------
% This function aligns two vectors that both had a fixed original length,
% but have been processed such that the specified indices have been
% removed. The resulting vectors have the same length, with NaN used to
% replace values corresponding to removed samples. Input arguments include:
% - vector1     [Mx1] First vector of samples
% - vector2     [Nx1] Second vector of samples
% - remIdx1     [Rx1] Removed indices from original first vector
% - remIdx2     [Sx1] Removed indices from original second vector
%
% This function returns:
% - alignedV1   [(M+R)x1]
% - alignedV2   [(N+S)x1]
% where (M+R) = (N+S).
% -------------------------------------------------------------------------

% Initialize placeholders for return values
alignedV1 = zeros(length(vector1)+length(remIdx1), 1);
alignedV2 = zeros(length(vector2)+length(remIdx2), 1);

% -------------------------------------------------------------------------
% Reconstruct Vector 1
% -------------------------------------------------------------------------

% Initialize counter
counter = 1;

for i = 1:length(alignedV1)
    
    % If the element was not removed, write the value from vector1
    if isempty(find(remIdx1 == i, 1))
        alignedV1(i) = vector1(counter); counter = counter + 1;
    else
        % Else, write NaN
        alignedV1(i) = NaN;
    end
    
end

% -------------------------------------------------------------------------
% Reconstruct Vector 2
% -------------------------------------------------------------------------

% Initialize counter
counter = 1;

for i = 1:length(alignedV2)
    
    % If the element was not removed, write the value from vector2
    if isempty(find(remIdx2 == i, 1))
        alignedV2(i) = vector2(counter); counter = counter + 1;
    else
        % Else, write NaN
        alignedV2(i) = NaN;
    end
    
end

end

