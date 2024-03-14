function obj = confusion(obj, varargin)

% -------------------------------------------------------------------------
% Description:
% This function generates a confusion matrix from the template set.
%
% Arguments (opt'l)
% - without                 List of indices of held-out templates
% - verbose FLAG            Print progress to console?
% - B                       Number of beats to average
% -------------------------------------------------------------------------

% Parse optional arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'without'); ho = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        elseif strcmp(varargin{arg}, 'B'); B = varargin{arg + 1};
        end
    end
end

% Set defaults
if ~exist('verbose', 'var'); verbose = true; end
if ~exist('B', 'var'); B = 100; end
if ~exist('ho', 'var'); hoflag = false; ho = 0; else; hoflag = true; end

if verbose; disp(" "); ...
        disp("Building Confusion Matrix"); disp("-------------------------"); end

% Process input arguments
C = obj.numClasses; S = obj.numSets; hi = 1:S; if hoflag; hi(ho) = []; end

% Subject-specific and cumulative performance
cumulative = zeros(C, C);	% Cumulative performance

% Make a directory for saving data
mkdir metadata

parfor s = hi
    
    % Subject-specific performance
    ss = zeros(C, C);
    
    if verbose; disp("Processing Set: " + string(s) + " of " + string(length(hi))); end
    
    % Get predictions
    if hoflag; [~, pred] = obj.predict('on', s, 'without', ho, 'nochar');
    else; [~, pred] = obj.predict('on', s, 'nochar');
    end
    
    % Populate placeholders
    for i = 1:C; for j = 1:C; ss(i, j) = sum(pred{j} == i); end; end
    
    % Save data
    smart.general.parsave("metadata/ss_" + string(s), 'ss', ss);
    
end

% Get the cumulative data
subspec = cell(S, 1);
for s = hi
    filename = "metadata/ss_" + string(s) + ".mat";
    load(filename, 'results'); subspec{s} = results.ss;
    cumulative = cumulative + subspec{s};
end

% Normalize columns
for i = 1:C; cumulative(:, i) = cumulative(:, i)./sum(cumulative(:, i)); end

% Write outputs
obj.setPerf = subspec; obj.cumulativePerf = cumulative; obj = obj.confidence(B);

end