function [output, remIdx] = removeRegion(signal)

% -------------------------------------------------------------------------
% Remove signal regions based on an imagesc plot of the data.
%
% Arguments (required)
% - signal      [MxN]   N signal segments of length M
% -------------------------------------------------------------------------

% Set placeholder for indices for removal
remIdx = [];

% Normalize the signal and plot on imagesc plot
figure; imagesc(normalize(signal)); hold on; grid on; colormap(flipud(gray));

% Ask user the number of regions to remove
numRegions = input("Number of regions for removal: ");

% For each region...
for i = 1:numRegions
    
    % Prompt the user to select a start and end point
    title("Select Region " + string(i) + " Start/End Points")
    
    % Return start and end points selected by the user
    [x, ~] = ginput(2); x = round(x);
    if x(1) < 1; x(1) = 1; end
    
    % Add region to removal indices placeholder
    remIdx = [remIdx x(1):x(2)];
    
end

% Remove specified indices from the signal
signal(:, remIdx) = []; output = signal;

end

