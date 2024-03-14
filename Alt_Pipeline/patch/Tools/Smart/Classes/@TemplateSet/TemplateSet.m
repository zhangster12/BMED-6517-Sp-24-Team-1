% Template Set Class

classdef TemplateSet
    
    % Class properties (semi-private)
    properties (SetAccess = public, GetAccess = public)
        data            % Raw datasets
        type Index      % SQI type of template set
        lambda          % Lambda parameter for SQI
        Fs              % Sampling frequency of data
        set             % Template sets
        setPerf         % Performance on training data (set-specific)
        cumulativePerf  % Performance on training data (cumulative)
        mu              % Mean SQI for each template
        sigma           % Covariance of each template in set
        error           % Probability of error given performance matrix
        numSets         % Number of template sets
        numClasses      % Number of templates per set / classes
        SQI             % Placeholder for cumulative SQIs for template set
    end
    
    % Class methods (public)
    methods
        
        % Class constructor
        function obj = TemplateSet(dataSet, sqiType, lambda, Fs)
            obj.data = dataSet; obj.type = sqiType; obj.lambda = lambda; obj.Fs = Fs;
        end
        
        % Create template set
        obj = create(obj, varargin)
        
        % Characterize template set (on all templates or a subset)
        obj = characterize(obj, varargin)
        
        % Generate predictions (on new data or a template set)
        [poll, predictions] = predict(obj, varargin)
        
        % Generate scores (on new data or a template set)
        scores = score(obj, varargin)
        
        % Determine prediction error of template set
        obj = confidence(obj, B)
        
        % Create confusion matrix
        obj = confusion(obj, B, varargin)
        
        % Return the mean score for each segment
        function scores = meanScore(~, allScores)
            temp = [];  % Compute average SQI for each segment
            for i = 1:length(allScores); temp = [temp allScores{i}{1}]; end
            scores = mean(temp, 2);
        end
        
    end
    
end