function outputArray = generateAlleles(inputArray, B)

% -------------------------------------------------------------------------
% This function divides a columnwise numerical array into a columnwise cell
% array, with each cell containing a specified number of array elements.
% The last cell may be abbreviated if the number of elements and elements
% per cell do no divide evenly. Input arguments include:
% - inputArray  [NxM]   A columnwise array of vectors
% - B                   Elements per cell, or vector of breakpoints
% -------------------------------------------------------------------------

% Set parameters
inputLen = size(inputArray, 1); % Length of input array
numInput = size(inputArray, 2); % Number of input arrays

% Set placeholder for return value
if isscalar(B)
    outputArray = cell(ceil(inputLen/B), numInput); % Evenly dividing bases
else
    outputArray = cell(length(B), numInput);    % Dividing by index
end

% Partition input vectors into output array
for i = 1:numInput                  % For each vector
    
    % If B is a scalar, partition the bases evenly into alleles
    if isscalar(B)
        
        currentIdx = 1;                 % Set index counter for for-loop

        for j = 1:ceil(inputLen/B)      % For each element

            if j < ceil(inputLen/B)     % For every element but the last one

                endIdx = currentIdx + (B-1);    % Update end index
                % Write values to output array
                outputArray{j, i} = inputArray(currentIdx:endIdx, i);
                currentIdx = endIdx + 1;    % Update current index

            else    % For the last element

                % Write values to output array (abbreviate result)
                outputArray{j, i} = inputArray(currentIdx:end, i);

            end

        end
    
    else
        
        currentIdx = 1; % Set index counter for for-loop
        
        % If B is a vector, generate alleles according to indices provided
        for index = 1:length(B)
           
            % For all indices but the last one...
            if index < length(B)
                % Write result to output array
                outputArray{index, i} = inputArray(currentIdx:B(index), i);
                currentIdx = B(index) + 1;  % Update current index
                % For the last index...
            else
                % Write result to output array
                outputArray{index, i} = inputArray(currentIdx:end, i);
            end
            
        end
        
    end
    
end

end

