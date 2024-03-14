% Class for extracting, organizing, and processing data from the
% hypovolemia dataset.

classdef Oink
    
    % Class properties (public)
    properties (SetAccess = public, GetAccess = public)
        description     string  % Description of object
    end
    
    % Class properties (semi-private)
    properties (SetAccess = public, GetAccess = public)
        % Data and metadata
        subjects                % Subjects included in object
        dataset                 % Struct containing processed data
        levels          Level   % List of blood volume levels
        path            string  % Path to dataset (if not on IRL server)
        Fs                      % Sampling frequency
        
        % TemplateSet objects for quality indexing
        templateSet             % Struct w/ TemplateSet object fields
        
        % FLAGS
        filtered        logical % Is the data filtered?
        beatSeparated   logical % Is the data beat-separated?
    end
    
    % Class methods (public)
    methods (Access = public)
        
        % Class constructor
        function obj = Oink(varargin)
            
            % It is possible to set the description and subjects using the
            % class constructor. Else, default values will be assigned.
            
            % Initialize flags
            subjectsFLAG = false; descriptionFLAG = false;
            
            % Process optional input arguments
            if ~isempty(varargin)
                for i = 1:length(varargin)
                    if ~isa(varargin{i}, 'string'); obj.subjects = varargin{i}; subjectsFLAG = true;
                    else; obj.description = varargin{i}; descriptionFLAG = true;
                    end
                end
            end
            
            % Assign defaults if no values were given
            if ~subjectsFLAG; obj.subjects = 1; end
            if ~descriptionFLAG; obj.description = "Hypovolemia Study Dataset"; end
            
            % Set sampling frequency
            obj.Fs = 2000; % Hz
            
        end
        
        % Extracting data for the specified subjects
        obj = extractData(obj, varargin)
        
        % Generate TemplateSet object for quality indexing of each signal
        obj = createTemplateSet(obj, varargin)
        
        % Assign a pre-made TemplateSet object for a certain modality
        function obj = assignTemplateSet(obj, set, modality)
            obj.templateSet.(string(modality)) = set;
        end
        
        % Assign a pre-made SQI vector for a certain subject/level/modality
        function obj = assignSQI(obj, sqi, subject, level, modality)
            obj.dataset{subject}.(string(level)).sqi.(string(modality)) = sqi;
        end
        
        % Score all segments in the selected levels and modalities
        % (Creates TemplateSet)
        obj = score(obj, varargin)
        
        % -----------------------------------------------------------------
        % Ground-Truth Feature Extraction Functions
        % -----------------------------------------------------------------
        
        % Extract rAO (PEP) and/or rAC from Aortic Pressure
        obj = trueAortic(obj, varargin)
        
        % Extract Pulse Transit Time (PTT) from Femoral and Aortic Pressure
        obj = truePTT(obj, varargin)
        
        % Calibrate PAT/PTT based on Aortic MAP (Wearable and Catheter)
        obj = calibrate(obj, varargin)
        
        % Extract Arterial Blood Pressure from pressure waveforms
        obj = trueMAP(obj, varargin)
        
        % -----------------------------------------------------------------
        % Wearable Sensor Feature Extraction Functions
        % -----------------------------------------------------------------
        
        % Extract rAO (PEP) and/or rAC from SCG
        obj = scgAortic(obj, varargin)
        
        % Extract Pulse Transit Time (PTT) from SCG and PPG
        obj = scgPTT(obj, varargin)
        
        % Extract Heart Rate Variability (HRV) from ECG
        obj = ecgHRV(obj, varargin)
        
        % Extract SCG and PPG amplitude 
        obj = amplitudes(obj, varargin)
        
        % Obtain respiratory variability for specified signal
        obj = respvar(obj, varargin)
        
        % -----------------------------------------------------------------
        % Data Formatting Functions
        % -----------------------------------------------------------------
        
        % Generate training data with held-out subjects
        [data, labels] = trainingData(obj, varargin)
        
        % Pull data from a specific sub-interval
        obj = getSubinterval(obj, fields, levels, varargin)
        
        % Determine the parent for a specified level
        obj = getParent(obj, level)
        
        % Normalize a feature vector by data from a specified level
        output = normalizeByLevel(obj, feature, varargin)
        
        % Save data in a parfor loop
        parsave(obj, filename, oink)
        
        % Stitch beat-separated data
        stitched = stitch(obj, field, level, varargin)
        
    end
    
    % Class methods (private)
    methods (Access = public)
        
        % Obtain sample number from time-step
        sample = timeToSample(obj, start_time, Fs, time)
        
        % Obtain beat number from sample number
        beat = sampleToBeat(obj, samples, beats)
        
    end
    
end