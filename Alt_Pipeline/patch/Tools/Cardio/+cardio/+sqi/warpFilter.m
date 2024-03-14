function [distances, signal_c, signal_r, outliers] = ...
    warpFilter(prototype, signals, tol, Fs, varargin)

% -------------------------------------------------------------------------
% SUMMARY
% This function uses dynamic time warping (DTW) to return a subset of signals
% representative of the template. The filter first breaks down the template
% and signal using Fourier series to obtain error coefficients. These
% coefficients are then weighted using the noise response of the given
% template to obtain an estimate of the relative DTW distance. The noise
% response is optionally modified by a mask that may be specified by the user
% (for instance, the counter-mask such that all noise components have equal
% influence on the distance, or a selective mask, such that some noise
% components are not considered).
%     
% ARGUMENTS (REQUIRED)
% - prototype     [Nx1]   Prototypical signal for template
% - signals       [NxM]   Raw signal set for comparison
% - tol                   Percentage of signals to retain
% - Fs                    Sampling frequency (Hz)
% 
% ARGUMENTS (OPTIONAL)
% - varargin
%     - 'M'               Coefficient for exponential moving average
%     - 'mask'    [Yx1]   Mask for noise response
%     - 'plot'    Bool    Plot accepted/rejected signals?
%     - 'title'           Title of data plot
%     - 'order'           Order of Fourier series decomposition
%     - 'windows'         Number of windows for Fourier series
%     - 'period'          Period for Fourier series
%     - 'response'    Bool    Characterize noise response of the signal?
%     - 'impulse'     Bool    Generate impulse response of template?
% 
% RETURNS
% - varargout
%     - signal_c  [NxA]   Accepted signal segments given tol
%     - signal_r  [NxR]   Rejected signal segments given tol
%     - distances [Mx2]   Projected and actual distances from template
%     - ranks     [Mx2]   Projected and actual distance ranks
%     - outliers  [Rx1]   Indices of rejected segments
% -------------------------------------------------------------------------

% Parse optional arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'M'); M = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'mask'); mask = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'plot'); Plot = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'title'); Title = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'order'); order = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'windows'); windows = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'period'); period = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'response'); impulse = varargin{arg + 1};
        end
    end
end

% Set defaults
if ~exist('M', 'var'); M = 1; end
if ~exist('Plot', 'var'); Plot = false; end
if ~exist('Title', 'var'); Title = []; end
if ~exist('order', 'var'); order = 8; end
if ~exist('windows', 'var'); windows = 8; end
if ~exist('period', 'var'); period = floor(length(prototype)/(2*windows)); end
if ~exist('impulse', 'var'); impulse = false; end

% -------------------------------------------------------------------------
% Signal Pre-Processing
% -------------------------------------------------------------------------

% Filter signals with exponential moving average
sig = cardio.general.ema(signals, M, false);

% Normalize amplitude of prototypical and sample signals
template = normalize(prototype); sig = normalize(sig);

% Decompose template using Fourier series
% Arguments: signal, order, windows, Fs
St = cardio.sqi.decompose(template, order, windows, Fs, 'period', period);
Rt = cardio.sqi.reconstruct(St);    % Reconstruct signal

% Get total number of signals
S = size(sig, 2);

% Set placeholders
d_hat = zeros(S, 1);    % Estimated relative for each signal
d_star = zeros(S, 1);   % Actual total DTW distance

% -------------------------------------------------------------------------
% Characterize Filter
% -------------------------------------------------------------------------

if impulse
    % Arguments: template, signals, Fs, varargin
    [response, ~, ~, ~] = cardio.sqi.impulse(template, 10000, Fs);
    response = response./min(response); % Relative response

    % From response, generate a mask (inverse mask)
    mask = response.^-1;
end

% -------------------------------------------------------------------------
% Compute Distances
% -------------------------------------------------------------------------

% For each signal for comparison...
for s = 1:S
    
    % Obtain the actual DTW distance for each signal
    d_star(s) = dtw(template, sig(:, s));
    
    if exist('mask', 'var')
        
        % Decompose the signal using Fourier series
        Sc = cardio.sqi.decompose(sig(:, s), order, windows, Fs, 'period', period);
        
        % Get the error between the signal and template
        D.offset = Sc.offset - St.offset;   % Offset
        for win = 1:windows                 % Coefficients
            D.coeffs{win} = Sc.coeffs{win} - St.coeffs{win};
        end
        
        % Organize the error vector such that a mask can be applied to it
        theta_e = getCoeffs(D);
        
        % Expand mask vector
        eMask = zeros(2*length(mask), 1); cc = 1;
        for i = 1:length(mask); eMask(cc) = mask(i); eMask(cc + 1) = mask(i); ...
                cc = cc + 2; end
        
        % Apply mask to error vector
        errorAdj = theta_e.*eMask;
        
        % Re-organize error vector 
        for win = 1:windows     % For each window...
            D.coeffs{win} = reshape(errorAdj(1:2*order), [2, order]);
            D.coeffs{win} = D.coeffs{win}'; errorAdj(1:2*order) = [];
        end
        
        % Reconstruct the signal from the template using the scaled error
        Rs = cardio.sqi.reconstruct(St, 'noise', D);
        
        % Get the adjusted DTW distance between the template and signal
        d_hat(s) = dtw(Rt.signal, Rs.signal);
        
    else
        
        % Set the estimate to d_star
        d_hat(s) = d_star(s);
        
    end
    
end

% -------------------------------------------------------------------------
% Remove Outliers
% -------------------------------------------------------------------------

% Sort the vectors based on their distances
[~, rank_hat] = sort(d_hat);    % Estimated relative distance
% [~, rank_star] = sort(d_star);  % Actual DTW distance

% Get the percentage of segments indicated by tol
outliers = find(rank_hat > round(tol*length(d_hat)));

% Remove outliers from signal array
signal_c = sig;  % Placeholder for return value
signal_c(:, outliers) = []; % Remove outliers

% Return rejected signals
signal_r = sig(:, outliers);

% -------------------------------------------------------------------------
% Visualize Results
% -------------------------------------------------------------------------

if Plot
    
    % Plot outliers
    figure; hold on; grid on;
    title(Title + " Rejected Samples (DTW)");
    xlabel("Timestep"); ylabel("Amplitude");
    for s = 1:length(outliers)
        plot(signal_r(:, s)./range(signal_r(:, s)))
    end
    plot(template, '--k', 'LineWidth', 2)   % Plot baseline
    
    % Plot accepted
    figure; hold on; grid on;
    title(Title + " Accepted Samples (DTW)")
    xlabel("Timestep"); ylabel("Amplitude")
    for s = 1:size(signal_c, 2)
        plot(signal_c(:, s)./range(signal_c(:, s)))
    end
    plot(template, '--k', 'LineWidth', 2)   % Plot baseline
    
    % Plot impulse response
    figure; hold on; grid on;
    title(Title + " Noise Response (DTW)")
    xlabel("Time"); ylabel("Frequency")
    imagesc(reshape(response, [order, windows])); colorbar
    
end

% -------------------------------------------------------------------------
% Organize Output Arguments
% -------------------------------------------------------------------------

% Set distance and rank outputs
distances = [d_hat d_star];

end