function varargout = compensatoryReserve(signals, templates, varargin)

% -------------------------------------------------------------------------
% This function is an in-house formulation of the Compensatory Reserve
% Index (CRI) developed by Convertino et al. It returns the CRI of a signal
% (or set of signals) as a function of the signals' distance from a set of
% templates, as defined by the SQI. The CRI returned by this function is
% intended to be included as part of a model of blood volume status
% estimation using PPG signals.
%
% Arguments (req'd)
% - signals         [NxM]   M signals of length N
% - templates       [NxL]   L template signals of length N
%
% Arguments (opt'l)
% - lambda                  Lambda parameter for SQI (default: 25)
% - distance        Index   Distance metric for SQI (default: DTFM)
% - Fs                      Sampling frequency (default: 2kHz)
% - verbose         FLAG    Print progress?
%
% Outputs
% - varargout{1}    [Mx1]   Compensatory Reserve Index (CRI)
% - varargout{2}    [MxL]   Signal Quality Index (SQI) per template
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        switch varargin{arg}
            case 'lambda'; lambda = varargin{arg + 1};
            case 'distance'; distance = varargin{arg + 1};
            case 'Fs'; Fs = varargin{arg + 1};
            case 'verbose'; verbose = true;
        end
    end
end

% Set defaults for optional input arguments
if ~exist('lambda', 'var'); lambda = 25; end
if ~exist('distance', 'var'); distance = Index.FMDTW; end
if ~exist('Fs', 'var'); Fs = 2000; end
if ~exist('verbose', 'var'); verbose = false; end

% Initialize constants
numTemplates = size(templates, 2);  % Number of templates
numSignals = size(signals, 2);      % Number of signals

% Initialize placeholder for SQI scores for each template
SQI = zeros(numTemplates, numSignals);

% Create waitbar, if necessary
if verbose; f = waitbar(0, 'Computing CRI...'); end

% For each template...
for t = 1:numTemplates
    % Update waitbar, if necessary
    if verbose; waitbar(t/numTemplates, f, "Processing Template " + ...
            string(t) + " of " + string(numTemplates)); end
    % Compute the SQI of each signal versus the template
    SQI(t, :) = cardio.sqi.filter(signals, Fs, inf, 'type', distance, ...
        'template', templates(:, t), 'lambda', lambda);
end

% Close waitbar, if necessary
if verbose; close(f); end

% Calculate the CRI as the average SQI for each signal across all templates
CRI = mean(SQI, 2);

% Return output arguments
if nargout > 0; varargout{1} = CRI; end
if nargout > 1; varargout{2} = SQI; end

end