function [outputArray1, outputArray2] = recombine(inputArray1, inputArray2, C)

% -------------------------------------------------------------------------
% This function creates a new output array from two input array by splicing
% them and crossing over a maximum of R times. The function accepts the
% following input arguments: 
% - inputArray1     {NxM}   Cell array containing M vectors of N alleles
% - inputArray2     {NxM}   Cell array containing M vectors of N alleles
% - C                       Maximum number of recombinations
% -------------------------------------------------------------------------

% Get parameters
len = size(inputArray1, 1); % Length of input arrays

% Set return values placeholders
outputArray1 = inputArray1; outputArray2 = inputArray2;

% Determine number of recombinations (0 <= recombinations <= C)
recombinations = randi([0, C]);

if C > 0 && recombinations > 0
    
    % Generate indices for recombination
    indices = randi(len-1, [recombinations, 1]);

    % For each index...
    for i = 1:recombinations

        % Splice arrays at index
        a1 = outputArray1(1:indices(i), :); b1 = outputArray1(indices(i)+1:end, :);
        a2 = outputArray2(1:indices(i), :); b2 = outputArray2(indices(i)+1:end, :);

        % Recombine arrays
        outputArray1(:, :) = [a2; b1];
        outputArray2(:, :) = [a1; b2];

    end

end

end

