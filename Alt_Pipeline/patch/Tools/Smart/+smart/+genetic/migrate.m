function new_populations = migrate(old_populations, pop, Mi, Plot)

% -------------------------------------------------------------------------
% This function swaps a maximum proportion of the current population with
% other populations being analyzed. Input parameters are as follows:
% - old_populations     {Px1}   Cell vector containing each population
% - pop                 Index of current population
% - Mi                  Maximum proportion of current population to swap
%                       0 < Mi < 1
% - Plot                Bool    Plot results?
% -------------------------------------------------------------------------

% Set placeholder for return value
new_populations = old_populations;

% Get number of populations
numPop = length(old_populations);

% Determine number of organisms to swap
totOrg = size(old_populations{pop}, 1); % Organisms available to swap
maxOrg = round(Mi*totOrg);              % Maximum organisms to swap
numOrg = randi([0, maxOrg]);            % Organisms to swap

% Record swaps
swaps = zeros(totOrg, numPop);

% Swap each organism
for i = 1:numOrg
    
    % Select organism at random to swap, from source and destination
    sourceIdx = randi(totOrg); destIdx = randi(totOrg);
    
    destPop = pop;  % Set destination to current population
    % Determine the population to which to send current organism
    % Repeat until the destination population isn't the current one
    while destPop == pop; destPop = randi(numPop); end
    
    % Get the organisms
    sourceOrg = old_populations{pop}{sourceIdx};
    destOrg = old_populations{destPop}{destIdx};
    
    % Swap the organisms
    new_populations{pop}{sourceIdx} = destOrg;
    new_populations{destPop}{destIdx} = sourceOrg;
    
    % Record the swap
    swaps(sourceIdx, pop) = swaps(sourceIdx, pop) + 1;
    swaps(destIdx, destPop) = swaps(destIdx, destPop) - 1;
    
end

if Plot
    
    figure(numPop + 1)   % Initialize figure
    % For each population...
    for i = 1:numPop
        % Create subplot
        subplot(1, numPop, i); hold on
        temp = swaps(:, i);     % Create placeholder for population
        imagesc(reshape(temp, [sqrt(totOrg), sqrt(totOrg)])); colorbar
        title("Swaps for Population " + string(i));
        xlim([0,sqrt(totOrg)+1]); ylim([0,sqrt(totOrg)+1]); hold off
    end
    
    
end

end

