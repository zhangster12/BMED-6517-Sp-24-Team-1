function [newSignal, newTemplate, signalMap, templateMap] = ...
    warpFeatures(signals, template, Fs, varargin)

% -------------------------------------------------------------------------
% SUMMARY
% This function adapts an input signal to match a template by first using
% dynamic time warping (DTW) to match peaks in the signals, and then
% projecting the signal into the Fourier series subspace to smooth
% discontinuities. Features extracted from the modified signal may be mapped
% back to the original signal using the generated map.
% 
% ARGUMENTS (REQUIRED)
% - signals   [NxM]	Signals to modify
% - template  [Nx1]	Template signal
% - Fs                Sampling frequency (Hz)
% 
% ARGUMENTS (OPTIONAL)
% - windows           Windows in which to divide the signal
% - order             Order of Fourier series
% - period            Period length for Fourier series
% - offset            Offset for Fourier series
% - M                 Smoothing factor for exponential moving average
% - mask              Weights for Fourier coefficients
% - offsetMask        Weights for Fourier offsets
% -------------------------------------------------------------------------

% Parse input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'windows'); windows = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'order'); order = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'period'); period = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'offset'); offset = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'M'); M = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'mask'); mask = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'offsetMask'); offsetMask = varargin{arg + 1};
        end
    end
end

% Set defaults
if ~exist('windows', 'var'); windows = 8; end
if ~exist('order', 'var'); order = 3; end
if ~exist('period', 'var'); period = floor(length(template)/(2*windows)); end
if ~exist('M', 'var'); M = 1; end
if ~exist('mask', 'var'); mask = ones(2*order*windows); end
if ~exist('offsetMask', 'var'); offsetMask = ones(windows, 1); end

% Set placeholders for return values
signalMap = cell(size(signals, 2), 1); templateMap = signalMap;
newSignal = cell(size(signals, 2), 1); newTemplate = newSignal;

% -------------------------------------------------------------------------
% Signal Pre-Processing
% -------------------------------------------------------------------------

% Normalize amplitude of prototypical and sample signals
template = normalize(template); signals = normalize(signals);

% Filter signals with exponential moving average
sig = cardio.general.ema(signals, M, false);

% Decompose template using Fourier series
% Arguments: signal, order, windows, Fs
if exist('offset', 'var')
    St = cardio.sqi.decompose(template, order, windows, Fs, 'period', period, 'offset', offset);
else
    St = cardio.sqi.decompose(template, order, windows, Fs, 'period', period);
end
Rt = cardio.sqi.reconstruct(St);    % Reconstruct signal

% Get total number of signals
S = size(sig, 2);

% For each signal...
for s = 1:S
    
    % ---------------------------------------------------------------------
    % Reconstruct Signal
    % ---------------------------------------------------------------------
    
    % Decompose the signal using Fourier series
    if exist('offset', 'var')
        Sc = cardio.sqi.decompose(sig(:, s), order, windows, Fs, 'period', period, 'offset', offset);
    else
        Sc = cardio.sqi.decompose(sig(:, s), order, windows, Fs, 'period', period);
    end
    
    % Get the error between the signal and template
    D.offset = offsetMask.*(Sc.offset - St.offset);   % Offset
    for win = 1:windows                 % Coefficients
        D.coeffs{win} = Sc.coeffs{win} - St.coeffs{win};
    end
    
    % Organize the error vector such that a mask can be applied to it
    theta_e = [];   % Initialize placeholder
    for win = 1:windows     % For each window...
        temp = D.coeffs{win}'; temp = temp(:); theta_e = [theta_e; temp];
    end
    
    % Apply mask to error vector
    errorAdj = theta_e.*mask;

    % Re-organize error vector 
    for win = 1:windows     % For each window...
        D.coeffs{win} = reshape(errorAdj(1:2*order), [2, order]);
        D.coeffs{win} = D.coeffs{win}'; errorAdj(1:2*order) = [];
    end

    % Reconstruct the signal from the template using the scaled error
    Rs = cardio.sqi.reconstruct(St, 'noise', D);
    
    % ---------------------------------------------------------------------
    % Dynamic Time Warping
    % ---------------------------------------------------------------------

    % Obtain a map of the original signal to the new signal
    [~, templateMap{s}, signalMap{s}] = dtw(Rt.signal, Rs.signal);

    % Construct the new signal from the map
    warpedSignal = Rs.signal(signalMap{s})';
    warpedTemplate = Rt.signal(templateMap{s})';

    % ---------------------------------------------------------------------
    % Signal reconstruction
    % ---------------------------------------------------------------------
    
    % Remove discontinuities using moving average
    newSignal{s} = movmean(warpedSignal, period/2);
    newTemplate{s} = movmean(warpedTemplate, period/2);
    
end

end

