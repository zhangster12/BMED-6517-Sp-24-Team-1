function scores = score(obj, varargin)

% -------------------------------------------------------------------------
% Description:
% This function constructs a template set and returns SQI scores.
%
% Arguments (opt'l)
% - data    {[MxN] x C}     Dataset for which to generate scores
% - on                      List of elements of "data" for which to generate scores
% - without                 List of indices of held-out templates for score generation
% - with                    List of indices of held-in templates for score generation
% - verbose FLAG            Print progress to console?
% - char    FLAG            Re-characterize model?
% -------------------------------------------------------------------------

if isempty(obj.set)
    disp("Error in TemplateSet.score():")
    disp(" -> Run TemplateSet.create() before computing scores")
    return
end

% Parse optional arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'data'); data = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'on'); on = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'without'); ho = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'with'); hi = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        elseif strcmp(varargin{arg}, 'char'); char = true;
        end
    end
end

% Process input arguments
C = obj.numClasses; S = obj.numSets;

% Set defaults
if ~exist('verbose', 'var'); verbose = false; end
if ~exist('char', 'var'); char = false; end
if ~exist('data', 'var'); data = [obj.data{:}]; end
if ~exist('on', 'var'); on = []; end
if ~exist('ho', 'var'); ho = []; end
if ~exist('hi', 'var'); hi = 1:S; end

% Determine which sets are held in
if ~isempty(ho); hi(ho) = []; end

% Modify dataset and template sets based on settings
if (~isempty(on) || ~isempty(ho)) && char
    if verbose; obj = obj.characterize('without', [on, ho], 'verbose');
    else; obj = obj.characterize('without', [on, ho]);
    end; hi([on, ho]) = [];
elseif char
    if verbose; obj = obj.characterize('verbose');
    else; obj = obj.characterize();
    end
elseif ~isempty(on) || ~isempty(ho); hi([on, ho]) = [];
end; if ~isempty(on); data = data(on, :); end

if verbose; disp(" "); ...
        disp("Generating Scores"); disp("-----------------"); end

% Assess the held-out subject
test_SQI = cell(S, 1); test_SQI(:) = {cell(C, C)};
for ts = hi
    
    if verbose; disp(" -> Progress: " + string(100*(ts-1)/length(hi)) + "%"); end
    
    % For each template in each set...
    for t = 1:C
        for c = 1:C
            test_SQI{ts}{t, c} = cardio.sqi.filter(data{c}, obj.Fs, inf, ...
                'template', obj.set{t}{ts}, 'lambda', obj.lambda, 'type', obj.type);
        end
    end
end

% Set return value
scores = test_SQI;

end