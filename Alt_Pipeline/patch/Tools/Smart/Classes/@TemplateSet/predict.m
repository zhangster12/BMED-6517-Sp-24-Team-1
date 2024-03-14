function [poll, pred] = predict(obj, varargin)

% -------------------------------------------------------------------------
% Description:
% This function returns a vector of predictions from a template ensemble.
%
% Arguments (opt'l)
% - data    {[MxN] x C}     Dataset for which to generate scores
% - on                      List of elements of "data" for which to generate predictions
% - without                 List of indices of held-out templates for prediction
% - with                    List of indices of held-in templates for prediction
% - verbose FLAG            Print progress to console?
% - nochar  FLAG            Do not re-characterize model?
% -------------------------------------------------------------------------

if isempty(obj.set)
    disp("Error in TemplateSet.predict():")
    disp(" -> Run TemplateSet.create() before computing predictions")
    return
end

% Parse optional arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'data'); data = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'on'); on = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'without'); ho = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'verbose'); verbose = true;
        elseif strcmp(varargin{arg}, 'nochar'); nochar = true;
        end
    end
end

% Set defaults
if ~exist('verbose', 'var'); verbose = false; end
if ~exist('nochar', 'var'); nochar = false; end
if ~exist('on', 'var'); on = []; end
if ~exist('ho', 'var'); ho = []; end
if ~exist('data', 'var'); data = [obj.data{:}]; end

% Process input arguments
C = obj.numClasses; S = obj.numSets;

% Determine which sets are held in
hi = 1:S;

if (~isempty(on) || ~isempty(ho)) && ~nochar
    if verbose; obj = obj.characterize('without', [on, ho], 'verbose');
    else; obj = obj.characterize('without', [on, ho]);
    end; hi([on, ho]) = [];
elseif ~nochar
    if verbose; obj = obj.characterize('verbose'); else; obj = obj.characterize(); end
elseif ~isempty(on) || ~isempty(ho); hi([on, ho]) = [];
end; if ~isempty(on); data = data(on, :); end

if verbose; disp(" "); ...
        disp("Generating Predictions"); disp("-----------------------"); end

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

% Initialize polling matrix for each true class
poll = cell(C, 1);

% Perform analysis for each true class
for c = 1:C
    
    disp(" -> Performing analysis for class " + string(c))

    % Get the number of beats and construct a polling matrix
    poll{c} = zeros(size(test_SQI{1}{1, c}, 1), S);

    % For each template set...
    for ts = hi
        % Get the vector of SQIs
        temp = [test_SQI{ts}{:, c}];

        % For each row in the result, get the most likely class
        for row = 1:size(temp, 1)

            % If the value is nan, continue and predict 1
            if ~isempty(find(isnan(temp(row, :)), 1)); poll{c}(row, ts) = 1; continue; end

            % For each possible class...
            max = 0; idx = 1;
            for class = 1:C
                % Get the probability
                prob = mvncdf(temp(row, :)' - 0.05, temp(row, :)' + 0.05, ...
                    obj.mu{ts}(:, class), obj.sigma{ts}{class});
                if prob > max; max = prob; idx = class; end
            end

            % Update the polling matrix
            poll{c}(row, ts) = idx;

        end

    end

end

% Consolidate the polling matrix into a prediction matrix
pred = cell(C, 1);
for c = 1:C; pred{c} = mode(poll{c}, 2); end

end