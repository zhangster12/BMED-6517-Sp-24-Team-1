function [SQI, signals_c, signals_r, rejected, features] = ...
    warpMask(signals, numFeatures, Fs, tol, Plot, varargin)

% -------------------------------------------------------------------------
% SUMMARY
% This function computes a signal quality index for a set of signals given a
% template. Signals are first projected into the subspace spanned by the
% Fourier series, with the coefficient weights being modified based on the
% properties of the DTW filter for the template, or applies a custom mask.
% The peaks of the template most consistent in the signal set are then
% identified, and the index of each feature is returned. 
% 
% ARGUMENTS (REQUIRED)
% - signals           [NxM]   Set of signal vectors
% - numFeatures               Number of features to extract
% - Fs                        Sampling frequency
% - tol                       Cutoff for SQI
% - Plot              Bool    Plot results?
% 
% ARGUMENTS (OPTIONAL)
% - varargin
%     - 'M'                   Smoothing factor for EMA
%     - 'windows'             Number of windows for Fourier series
%     - 'order'               Order number for Fourier series
%     - 'period'              Period for Fourier series
%     - 'offset'              Offset for Fourier series
%     - 'maxStretch'          Maximum feature stretching in DTW
%     - 'maxMatch'            Maximum distance for feature matching
%     - 'maskType'            Mask type for coefficient weights
%         - 'uniform'         Equal importance for all coefficients
%         - 'inverse'         Inverse of coefficient importance
%         - 'fronthalf'       Considering front-half coefficients only
%         - 'backhalf'        Considering back-half coefficients only
% -------------------------------------------------------------------------

% Parse input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'windows'); windows = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'order'); order = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'period'); period = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'offset'); offset = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'M'); M = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'maskType'); maskType = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'maxStretch'); maxStretch = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'maxMatch'); maxMatch = varargin{arg + 1};
        end
    end
end

% Set defaults
if ~exist('windows', 'var'); windows = 8; end
if ~exist('order', 'var'); order = 3; end
if ~exist('period', 'var'); period = floor(size(signals, 1)/(2*windows)); end
if ~exist('M', 'var'); M = 1; end
if ~exist('maskType', 'var'); maskType = 'uniform'; end
if ~exist('maxStretch', 'var'); maxStretch = 50; end
if ~exist('maxMatch', 'var'); maxMatch = 25; end

% -------------------------------------------------------------------------
% Generate Template and Align Signals
% -------------------------------------------------------------------------

% Arguments: signals, Fs, varargin
[normSignals, templates] = ...
    cardio.general.template(signals, Fs, 'maxlag', 0.1, 'smooth', true);

% Fetch overall template
normTemplate = templates(:, end);

% -------------------------------------------------------------------------
% Generate Mask
% -------------------------------------------------------------------------

if strcmp(maskType, 'uniform')
    
    % Generate uniform mask
    mask = ones(2*order*windows, 1);
    offsetMask = ones(windows, 1);
    
elseif strcmp(maskType, 'inverse')
    
    % Characterize DTW for template type
    % Arguments: template, signals, Fs, varargin
    [response, ~, ~, ~] = ...
        cardio.sqi.impulse(normTemplate, 10000, Fs, 'order', order, 'windows', windows);
    response = response./min(response); % Relative response
    
    % Generate inverse mask
    temp = response.^-1;
    
    % Expand mask
    mask = zeros(2*length(temp), 1); cc = 1;
    for i = 1:length(temp); mask(cc) = temp(i); mask(cc + 1) = temp(i); ...
            cc = cc + 2; end
    
    % Set offset mask
    offsetMask = zeros(windows, 1);
    
elseif strcmp(maskType, 'fronthalf')
    
    % Generate front-half mask
    mask = ones(2*order*windows, 1); offsetMask = ones(windows, 1);
    mask((order*windows + 1):end) = 0; offsetMask(floor(windows/2)+1:end) = 0;
    
elseif strcmp(maskType, 'backhalf')
    
    % Generate back-half mask
    mask = zeros(2*order*windows, 1); offsetMask = zeros(windows, 1);
    mask((order*windows + 1):end) = 1; offsetMask(floor(windows/2)+1:end) = 1;
    
end

% -------------------------------------------------------------------------
% Get Template Peaks
% -------------------------------------------------------------------------

% Get peaks and valleys from template signal and consolidate into raw ft.
midway = floor(length(normTemplate)/2); % Get midway point in signal
if strcmp(maskType, 'fronthalf')
    [peaks, valleys, ~] = cardio.general.getPeaks(normTemplate(1:midway));
elseif strcmp(maskType, 'backhalf')
    [peaks, valleys, ~] = cardio.general.getPeaks(normTemplate(midway:end));
else
    [peaks, valleys, ~] = cardio.general.getPeaks(normTemplate);
end

% Template features
templateFt = [peaks; valleys];

% Set placeholder for feature and presence vectors
presence = zeros(size(signals, 2), length(templateFt));
features = presence;

% -------------------------------------------------------------------------
% Warp Signals
% -------------------------------------------------------------------------

% Arguments: signals, template, Fs, varargin
if exist('offset', 'var')
    [~, ~, sIdx, tIdx] = ...
        cardio.sqi.warpFeatures(signals, normTemplate, Fs, 'windows', windows, 'order', order, ...
        'period', period, 'M', M, 'mask', mask, 'offsetMask', offsetMask, 'offset', offset);
else
    [~, ~, sIdx, tIdx] = ...
        cardio.sqi.warpFeatures(signals, normTemplate, Fs, 'windows', windows, 'order', order, ...
        'period', period, 'M', M, 'mask', mask, 'offsetMask', offsetMask);
end

% -------------------------------------------------------------------------
% Match Peaks
% -------------------------------------------------------------------------

% For each signal...
for s = 1:size(signals, 2)
    
    % Find the features of the raw signal
    [peaks, valleys, ~] = cardio.general.getPeaks(normSignals(:, s));
    signalFt = [peaks; valleys];
    
    % Map features to warped signals
    warpedTemplateFt = cardio.general.map(templateFt, tIdx{s});
    warpedSignalFt = cardio.general.map(signalFt, sIdx{s});
    
    % For each warped template feature...
    for f1 = 1:length(warpedTemplateFt)
        
        % Find the closest feature in the warped signal
        distance = inf; % Placeholder for closest distance
        idx = 0;        % Placeholder for closest index
        for f2 = 1:length(warpedSignalFt)
            warpedDist = abs(warpedSignalFt(f2) - warpedTemplateFt(f1));
            rawDist = abs(signalFt(f2) - templateFt(f1));
            if warpedDist < distance && rawDist < maxStretch
                distance = abs(warpedSignalFt(f2) - warpedTemplateFt(f1));
                idx = signalFt(f2);
            end
        end
        
        % If the closest feature is close enough, there is a match
        if distance < maxMatch
            presence(s, f1) = 1; features(s, f1) = idx;
        end
        
    end
    
end

% -------------------------------------------------------------------------
% Select Features
% -------------------------------------------------------------------------

% Compute feature quality
featureQuality = mean(presence);

% Select top features
[~, rank] = sort(featureQuality, 'descend');
topFeatures = rank(1:numFeatures);
rejectedFeatures = rank(numFeatures + 1:end);

% -------------------------------------------------------------------------
% Determine SQI
% -------------------------------------------------------------------------

% Set placeholder for return value
SQI = zeros(size(signals, 2), 1);

% For each signal...
for s = 1:size(signals, 2)
    
    % Determine how many of the top features it has
    numFt = 0;
    for f = 1:numFeatures
        numFt = numFt + presence(s, topFeatures(f));
    end
    
    % Write result to SQI
    SQI(s) = numFt/numFeatures;
    
end

% Sort accepted and rejected segments
accepted = find(SQI >= tol); rejected = find(SQI < tol);

% -------------------------------------------------------------------------
% Visualize Results
% -------------------------------------------------------------------------

if Plot
    
    % Plot accepted segments, annotating the key points
    figure; hold on; grid on; title("Accepted Segments")
    for s = 1:length(accepted)
        plot(normSignals(:, accepted(s)))
        for f = 1:numFeatures
            xPoint = features(accepted(s), topFeatures(f));
            if xPoint ~= 0
                plot(xPoint, normSignals(xPoint, accepted(s)), 'ok')
            end
        end
    end

    % Plot rejected segments
    figure; hold on; grid on; title("Rejected Segments")
    for s = 1:length(rejected)
        plot(normSignals(:, rejected(s)))
        for f = 1:numFeatures
            xPoint = features(rejected(s), topFeatures(f));
            if xPoint ~= 0
                plot(xPoint, normSignals(xPoint, rejected(s)), 'ok')
            end
        end
    end
    
end

% -------------------------------------------------------------------------
% Return Values
% -------------------------------------------------------------------------

signals_c = signals(:, accepted);
signals_r = signals(:, rejected);
features(accepted, :) = []; features(:, rejectedFeatures) = [];
features(features == 0) = NaN;

end