function output = confusion(predictions, targets)

% This function creates a confusion matrix given vectors of predicted and
% target classes. Both vectors should have contiguous integer values.

% Set placeholder for return value
output = zeros(length(unique(targets)));

% Fill the confusion matrix
for i = 1:length(targets)
    output(targets(i), predictions(i)) = output(targets(i), predictions(i)) + 1;
end

end

