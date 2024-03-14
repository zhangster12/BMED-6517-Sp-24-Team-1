% Saving variables from within a parfor loop
function parsave(subject, varargin)

% Set save file name
savefile = string(subject) + ".mat";

% Get variable names and values
names = varargin(1:2:(length(varargin) - 1));
values = varargin(2:2:length(varargin));

% Create a struct with names and values
for i = 1:length(varargin)/2
    results.(names{i}) = values{i};
end

% Save results
save(savefile, 'results')

end