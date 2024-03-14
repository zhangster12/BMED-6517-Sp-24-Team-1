function outputArray = mutate(inputArray, M)

% -------------------------------------------------------------------------
% This function creates a new output array by mutating an input array a
% maximum of M times. A mutation is defined as swapping an element of two
% columns of the input cell array. Input arguments include:
% - inputArray  {NxM}   Cell array containing M vectors of N alleles
% - M           Maximum number of mutations
% -------------------------------------------------------------------------

% Get parameters
len = size(inputArray, 1); % Length of input arrays
num = size(inputArray, 2); % Number of input arrays

% Set return value placeholder
outputArray = inputArray;

% Determine number of mutations (0 <= mutations <= M)
mutations = randi([0, M]);

% For each vector...
if M > 0 && mutations > 0
    
    % Generate indices for mutations
    indices = randi(len-1, [mutations, 1]);
    
    % For each index...
    for j = 1:mutations

        % Mutate array at index
        fromIndex = randi(num); toIndex = randi(num);   % Vectors between which to swap values
        fromVal = outputArray(indices(j), fromIndex);   % Value to swap from first vector
        toVal = outputArray(indices(j), toIndex);       % Value to swap from second vector
        outputArray(indices(j), fromIndex) = toVal;     % Overwrite first vector with second value
        outputArray(indices(j), toIndex) = fromVal;     % Overwrite second vector with first value

    end

end

end

