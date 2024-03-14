function [response, performance, features, targets] = impulse(template, signals, Fs, varargin)

% -------------------------------------------------------------------------
% SUMMARY
% This function models the impulse response of distance filters for SCG. The
% error parameters in the Fourier series domain are mapped to the final
% distance, with the coefficient assigned to each feature serving as an
% importance weight. The features are randomly divided into training and
% testing sets, with the linear model being trained on the training set and
% the performance being evaluated on the testing set. 
% 
% ARGUMENTS (REQ'D)
% - template          [Nx1]   Template signal vector
% - signals                   Number of synthetic signals to generate
% - Fs                        Sampling frequency (Hz)
% 
% ARGUMENTS (OPT)
% - varargin
%     - 'type'        INDEX   Distance metric
%     - 'order'               Order of Fourier series
%     - 'windows'             Number of signal windows
%     - 'period'              Period of Fourier series
%     - 'training'            Percent of synthetic set to train upon
%     - 'testing'     [NxM]   Testing dataset
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Pre-Processing
% -------------------------------------------------------------------------

% Parse arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'type'); type = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'order'); order = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'windows'); windows = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'period'); period = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'training'); training = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'testing'); testing = varargin{arg + 1};
        end
    end
end

% Set defaults
if ~exist('type', 'var'); type = Index.FMDTW; end
if ~exist('order', 'var'); order = 8; end
if ~exist('windows', 'var'); windows = 8; end
if ~exist('period', 'var'); period = floor(length(template)/(2*windows)); end
if ~exist('training', 'var'); training = 0.7; end

% Set placeholders
sig = zeros(length(template), signals);
% features = zeros(2*order*windows, signals);
features = zeros(order*windows, signals);
l2 = zeros(signals, 1); tw = zeros(signals, 1); sqi = zeros(signals, 1);
% testFeatures = zeros(2*order*windows, signals);
testFeatures = zeros(order*windows, signals);

% -------------------------------------------------------------------------
% Feature Extraction
% -------------------------------------------------------------------------

% Decompose template signal
% Arguments: signal, order, windows, Fs, varargin
St = cardio.sqi.decompose(template, order, windows, Fs, 'period', period, 'offset', 0);
Rt = cardio.sqi.reconstruct(St);

% Extract the testing features, if signals are provided
% Arguments: template, signal, Fs, order, windows
if exist('testing', 'var')
    for i = 1:size(testing, 2)
        [testFeatures(:, i), l2(i), tw(i)] = ...
            cardio.sqi.extractData(St, testing, Fs, order, windows);
    end
end

% Generate the specified number of synthetic signals
% Arguments: template, Fs, order, windows, varargin
for i = 1:signals
    % Generate synthetic signal
    [sig(:, i), features(:, i)] = cardio.sqi.synthetic(Rt.signal, Fs, order, windows, 'period', period);
    % Get distance
    if type == Index.Norm; l2(i) = norm(template - sig(:, i)); end
    if type == Index.DTW; tw(i) = dtw(template, sig(:, i)); end
    if type == Index.FMDTW; sqi(i) = cardio.sqi.filter(sig(:, i), Fs, inf, 'template', template); end
end
features = features';   % Obtain proper dimensionality for features

% Choose random indices to map
indices = randperm(signals);                    % Random indices
if ~exist('testing', 'var')
    trainCutoff = floor(training*signals);      % Number of training indices
    trainIndices = indices(1:trainCutoff);      % Training indices
    testIndices = indices(trainCutoff + 1:end); % Testing indices
else
    trainIndices = indices;                     % Train on all
    testIndices = randperm(size(testing, 2));   % Test on all
end

% Train on the training indices
if type == Index.Norm
    x = features(trainIndices, :)\l2(trainIndices);
elseif type == Index.DTW
    x = features(trainIndices, :)\tw(trainIndices);
elseif type == Index.FMDTW
    x = features(trainIndices, :)\sqi(trainIndices);
end

% Test on testing indices
if ~exist('testing', 'var')
    y = features(testIndices, :)*x;
else; y = testFeatures*x;
end

% Get distance between predictions and targets
if type == Index.Norm
    performance = cardio.general.rsquared(l2(testIndices), y);
elseif type == Index.DTW
    performance = cardio.general.rsquared(tw(testIndices), y);
elseif type == Index.FMDTW
    performance = cardio.general.rsquared(sqi(testIndices), y);
end

% Get impulse response
response = x;

% Get targets
if type == Index.Norm; targets = l2;
elseif type == Index.DTW; targets = tw;
elseif type == Index.FMDTW; targets = sqi;
end

end