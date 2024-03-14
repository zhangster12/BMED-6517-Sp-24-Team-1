function [tHR, HR, threshold] = automatedHR(ecg, varargin)

% -------------------------------------------------------------------------
% This function computes the heart rate of an ECG waveform with the desired
% precision. The function selects a threshold that minimizes the heart rate
% variability of the extracted heart rate.
%
% Arguments (required)
% - ecg             [Nx1]   ECG signal
% Arguments (optional)
% - numThresholds           Number of thresholds to test
% - Fs                      Sampling frequency (Hz)
% - minDist                 Minimum distance between R-peaks
% - plot            FLAG    Plot resutls?
% - corrected       FLAG    Correct heart rate for outliers?
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        switch varargin{arg}
            case 'numThresholds'; numThresholds = varargin{arg + 1};
            case 'Fs'; Fs = varargin{arg + 1};
            case 'minDist'; minDist = varargin{arg + 1};
            case 'plot'; Plot = true;
            case 'corrected'; corrected = true;
        end
    end
end

% Set defaults for optional input arguments
if ~exist('numThresholds', 'var'); numThresholds = 100; end
if ~exist('Fs', 'var'); Fs = 2000; end
if ~exist('minDist', 'var'); minDist = 500; end
if ~exist('Plot', 'var'); Plot = false; end
if ~exist('corrected', 'var'); corrected = false; end

% Set threshold search parameters
thresholds = 0.5*max(ecg):-max(ecg)/numThresholds:0;

% Placeholder for standard deviation and corresponding counter
stddev = zeros(size(thresholds)); counter = 1;

for i = 1:length(thresholds)

    % Set parameters
    threshold = thresholds(i);  % Threshold for detection

    % Compute heart rate (corrected)
    [~, HR] = cardio.ecg.computeHR(ecg, threshold, 'minDist', minDist, 'Fs', Fs, 'corrected');

    % Obtain standard deviation for heart rate and increment counter
    stddev(counter) = std(HR); counter = counter + 1;

end

% Find the index of the minimum standard deviation
% Format vector to find nonzero min
stddev_f = stddev(~isnan(stddev)); stddev_f = stddev_f(stddev_f > 0);
[~, dd_index] = min(stddev_f);

if Plot 
    
    % Plot standard deviation with chosen threshold labeled
    figure; hold on; grid on;
    % Format x axis, plot, and format
    x_axis = length(stddev_f)-1:-1:0;
    plot(x_axis, stddev_f); plot(x_axis(dd_index), stddev_f(dd_index), 'ok')
    title('Standard Deviation Inflection Point')
    xlabel('Threshold'); ylabel('Standard Deviation'); hold off
    
end

% Re-compute the heart rate with the chosen threshold
threshold = thresholds(stddev == stddev_f(dd_index)); threshold = threshold(1);
if corrected
    [tHR, HR] = cardio.ecg.computeHR(ecg, threshold, 'minDist', minDist, 'Fs', Fs, 'corrected');
else
    [tHR, HR] = cardio.ecg.computeHR(ecg, threshold, 'minDist', minDist, 'Fs', Fs);
end

if Plot
    
    % Plot HR for the chosen threshold
    figure; plot(tHR, HR); hold on;
    title('HR for threshold ' + string(x_axis(dd_index)))
    xlabel('Time (s)'); ylabel('HR (bpm)'); hold off
    
end

end

