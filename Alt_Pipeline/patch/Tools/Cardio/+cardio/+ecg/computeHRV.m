function [hrv, varargout] = computeHRV(indices, Fs, varargin)

% -------------------------------------------------------------------------
% This function computes heart rate variability using one of three methods:
% 1. R-R interval RMS difference ('difference')
% 2. Power Spectral Density - High vs. Low Frequency ('spectral')
% 3. Poincare Method ('poincare')
%
% Arguments (required)
% - indices     [Nx1]   R-peak time indices
% - Fs                  Sampling frequency (Hz)
%
% Arguments (optional)
% - difference  FLAG    Use window method for HRV (default)
% - spectral    FLAG    Use spectral method for HRV
% - poincare    FLAG    Use Poincare method for HRV
% - tol                 Tolerance for outliers (in standard deviations)
%
% Outputs
% - hrv         [Nx1]   Heart rate variability
% - LF/Cov1     [Nx1]   Low frequency component (spectral) or first
%                       eigenvalue (poincare)
% - HF/Cov2     [Nx1]   High frequency component (spectral) or second
%                       eigenvalue (poincare)
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'Fs'); Fs = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'difference'); difference = true;
        elseif strcmp(varargin{arg}, 'spectral'); spectral = true;
        elseif strcmp(varargin{arg}, 'poincare'); poincare = true;
        elseif strcmp(varargin{arg}, 'tol'); tol = varargin{arg + 1};
        end
    end
end

% Set defaults for optional input arguments
if ~exist('window', 'var'); difference = false; end
if ~exist('spectral', 'var'); spectral = false; end
if ~exist('poincare', 'var'); poincare = false; end
if ~difference && ~spectral && ~poincare; difference = true; end
if ~exist('tol', 'var'); tol = inf; end

% If the concurrent ECG and Fs were not provided for the spectral method,
% return an error message to the console
if spectral && ~exist('Fs', 'var')
    disp("-> Error in computeHRV(): Fs must be provided when using spectral method")
    hrv = nan; return
end

% Remove nans
indices = indices(~isnan(indices));

% Get difference between R peak indices
intervals = diff(indices); intervals = intervals./Fs;

% Remove outliers from the intervals, if indicated
if ~isinf(tol)
    meanint = mean(intervals); stdint = std(intervals);
    intervals(intervals < meanint - tol*stdint | intervals > meanint + tol*stdint) = nan;
    intervals = fillmissing(intervals, 'linear');
end
    
% Correct for edge cases
intervals(intervals < 0) = 0;

% -------------------------------------------------------------------------
% RMS Difference
% -------------------------------------------------------------------------

% Compute HRV as the RMS of R-R intervals
if difference; hrv = rms(intervals); end

% Correct for output errors
if difference && nargout > 1; varargout{2} = []; end
if difference && nargout > 2; varargout{2} = []; end

% -------------------------------------------------------------------------
% Power Spectral Density (HF vs. LF)
% -------------------------------------------------------------------------

% RATIONALE
% http://www.medscape.com/viewarticle/536280_3
% 1980s, Akselrod et al.[2] showed that the low frequency (LF) band
% (0.04-0.15 Hz) is related to both sympathetic and parasympathetic
% modulation, and the high frequency(HF) band (0.15-0.40 Hz) is governed
% almost exclusively by parasympathetic effects. The ratio of LF to HF power
% is often used as a metric of sympathetic-parasympathetic balance. It is
% important to note,however, that the main driver of HRV in the HF band is
% respiration,which produces the vagally mediated respiratory sinus
% arrhythmia. The magnitude of HF power is highly dependent on the depth of
% respiration,which often varies greatly from one recording epoch to
% another. Unless the depth of respiration is taken into account, assessment
% of autonomic balance by HRV measurement alone may be quite
% difficult.

if spectral
    
    % Determine the total length of time of recording
    ANN = cumsum(intervals) - intervals(1);
    
    % Compute the next nearest power of 2 of the window length and the
    % resulting updated sampling frequency
    NFFT = 2^nextpow2(length(intervals)); newFs = NFFT/ANN(end);
    
    % Interpolate the RR interval such that its length is L
    RR_interp = interp1(ANN, intervals, linspace(0, ANN(end), NFFT), 'spline');
    
    % Set parameters for power spectral density (PSD) calculation
    x = RR_interp; F = (newFs/2)*linspace(0, 1, NFFT/2+1);
    
    % Mean-center the RR-intervals and detrend
    x = x - mean(x); x = detrend(x); 

    % Obtain PSD of interpolated R-R intervals
    Y = fft(x, NFFT); PSD = 2*abs(Y(1:NFFT/2 + 1));
    
    % Compute the high- and low-frequency components
    VLF = sum(PSD(F < 0.04));
    LF = sum(PSD(F < 0.15)) - VLF;
    HF = sum(PSD(F < 0.40)) - LF - VLF;
    total = sum(PSD(F < 0.40)) - VLF;
    
    % Compute high- and low-frequency power as a percentage of total
    Apsd_LF = LF/total;
    Apsd_HF = HF/total;

    % Divide the PSD area of the low-frequency band by the high-frequency band
    % to obtain the ratio
    hrv = Apsd_LF/Apsd_HF;
    
    % Output HF and LF if indicated
    if nargout > 1; varargout{1} = Apsd_LF; end
    if nargout > 2; varargout{2} = Apsd_HF; end

end

% -------------------------------------------------------------------------
% Poincare Method (w/ Covariance for Quantification)
% -------------------------------------------------------------------------

if poincare

    % Separate points into corresponding x and y values
    x_vals = zeros(length(intervals)-1,1);
    y_vals = zeros(length(intervals)-1,1);

    % x_vals is the current interval (i), y_vals is the interval i+1
    for i = 1:length(intervals)-1
        x_vals(i) = intervals(i); y_vals(i) = intervals(i+1);
    end

    % Obtain the covariance matrix of the data and return the eigenvalues to
    % obtain the eigenvector length for the covariance matrix
    RR_covariance = cov([x_vals y_vals]);
    RR_eigenvalues = eig(RR_covariance);
    Cov1 = RR_eigenvalues(1);   % Return first eigenvalue
    Cov2 = RR_eigenvalues(2);   % Return second eigenvalue
    
    % Return the HRV as the ratio between the eigenvalues
    hrv = Cov1/Cov2;
    
    % Return Cov1 and Cov2 if indicated
    if nargout > 1; varargout{1} = Cov1; end
    if nargout > 2; varargout{2} = Cov2; end

end



end

