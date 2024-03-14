function [SQI, varargout] = ...
    filter(signals, Fs, tol, varargin)

% -------------------------------------------------------------------------
% SUMMARY
% This function computes the signal quality index as the percentage of
% template featurs matched by signal features, penalized by signal noise.
% 
% ARGUMENTS (REQ'D)
% - signals   [NxM]   Signal vectors
% - Fs                Sampling frequency (Hz)
% - tol               SQI cutoff
% 
% ARGUMENTS (OPT'L)
% - 'plot'    FLAG    Plot results?
% - 'type'    Index   Distance calculation
% - 'tolType'         Type of tolerance cutoff to use
%     - 'raw'         Cutoff is a raw SQI score
%     - 'percentage'  Cutoff is a percentage of segments
% - 'template' [Nx1]  Template signal
% - 'lambda'          Lambda hyperparameter for SQI
% - 'modelobj'        Model object with which to predict SQI
% - 'theta'           Parameters to classify with model
% - 'maxDist'         Maximum distance for DTFM feature matching
% -------------------------------------------------------------------------
    
% Parse optional arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'tolType'); tolType = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'type'); type = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'plot'); Plot = true;
        elseif strcmp(varargin{arg}, 'template'); template = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'lambda'); lambda = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'modelobj'); modelobj = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'theta'); theta = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'maxDist'); maxDist = varargin{arg + 1};
        end
    end
end

% Set defaults
if ~exist('tolType', 'var'); tolType = 'raw'; end
if ~exist('type', 'var'); type = Index.FMDTW; end
if ~exist('Plot', 'var'); Plot = false; end
if ~exist('lambda', 'var'); lambda = 1; end
if ~exist('maxDist', 'var'); maxDist = 50; end

% -------------------------------------------------------------------------
% Generate Template and Align Signals
% -------------------------------------------------------------------------

if ~exist('template', 'var') && type == Index.Model
    
    % Generate template signals
    [signals, templates] = ...
    cardio.general.template(signals, Fs, 'maxlag', 0.1, 'smooth', true);

    % Fetch overall template
    template = templates(:, end);
    
elseif type == Index.Model
    
    % Generate template signals
    [signals, templates] = ...
    cardio.general.template(signals, Fs, 'maxlag', 0.1, 'smooth', true);
    
    % Bring the template into alignment
    [template, ~, ~, ~] = cardio.general.alignSignals(template, ...
        templates(:, end), Fs, false, 'method', 'pad', 'modify', 'consistent');
    
end

% -------------------------------------------------------------------------
% Compute SQI
% -------------------------------------------------------------------------

SQI = zeros(size(signals, 2), 1);           % Set placeholder

% If a ML model was provided, use the model
if type == Index.Model
    % Predict the SQI
    SQI = 1- predict(modelobj, theta);
else
    
    % Compute signal distance
    if type == Index.FMDTW
        % Compute distance for modified SQI
        % Arguments: signals, template, varargin
        dist = cardio.sqi.peakMatch(signals, template, 'maxDist', maxDist);
    end


    for s = 1:size(signals, 2)                  % For each signal...
        % Compute modified or traditional SQI
        if type == Index.FMDTW; SQI(s) = exp(-lambda*dist(s));	% Modified
        else
            try
                temp = dtw(signals(:, s), template);
                SQI(s) = exp(-temp/length(signals(:, s)));      % Traditional
            catch
                % Preventing termination in execution when signal is invalid
                SQI(s) = nan;
            end
        end
    end
    
end

% Constrain SQI to 0 - 1
SQI(SQI < 0) = 0; SQI(SQI > 1) = 1;

% -------------------------------------------------------------------------
% Filter Segments
% -------------------------------------------------------------------------

if ~isinf(tol)
    % Determine accepted and rejected segments
    if strcmp(tolType, 'raw')
        accepted = find(SQI >= tol); rejected = find(SQI < tol);
        signals_c = signals(:, accepted); signals_r = signals(:, rejected);
    else
        numAccepted = floor(tol*size(signals, 2));	% Number accepted
        [~, rank] = sort(SQI, 'descend');           % Rank segments
        accepted = sort(rank(1:numAccepted), 'ascend');     % Accepted
        rejected = sort(rank(numAccepted+1:end), 'ascend'); % Rejected
        signals_c = signals(:, accepted); signals_r = signals(:, rejected);
    end
else
    accepted = 1:size(signals, 2); rejected = [];
    signals_c = signals; signals_r = [];
end

% -------------------------------------------------------------------------
% Visualize Results
% -------------------------------------------------------------------------

if Plot
    
    % Plot Accepted
    figure; hold on; grid on; title("Accepted Segments")
    xlabel("Segment"); ylabel("Amplitude"); plot(signals(:, accepted))
    
    % Plot Rejected
    figure; hold on; grid on; title("Rejected Segments")
    xlabel("Segment"); ylabel("Amplitude"); plot(signals(:, rejected))
    
    % Plot SQI
    figure; hold on; grid on; title("Signal Quality Index")
    xlabel("Segment"); ylabel("SQI"); plot(movmean(SQI, 25))
    
end

% -------------------------------------------------------------------------
% Format Output
% -------------------------------------------------------------------------

varargout{1} = signals_c; varargout{2} = signals_r;
varargout{3} = accepted; varargout{4} = rejected;

end

