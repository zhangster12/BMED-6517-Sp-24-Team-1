function [t_HR, RR_bpm] =...
    computeHR(ECG, Threshold, varargin)

% -------------------------------------------------------------------------
% This function computes the heart rate given a threshold.
%
% Arguments (required)
% - ecg         [Nx1]   ECG signal vector
% - threshold           Threshold for R-peak detection
%
% Arguments (optional)
% - minDist             Minimum distance between R-peaks
% - Fs                  Sampling frequency (Hz)
% - plot        FLAG    Plot results?
% - corrected   FLAG    Correct HR (omit outliers and interpolate)
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        switch varargin{arg}
            case 'minDist'; MinDist = varargin{arg + 1};
            case 'Fs'; Fs = varargin{arg + 1};
            case 'plot'; Plot = true;
            case 'corrected'; corrected = true;
        end
    end
end

% Set defaults for optional arguments
if ~exist('MinDist', 'var'); MinDist = 800; end
if ~exist('Fs', 'var'); Fs = 2000; end
if ~exist('Plot', 'var'); Plot = false; end
if ~exist('corrected', 'var'); corrected = false; end

% Identify peaks (indices)
[~,Rpeakslocs] = findpeaks(ECG, 'MinPeakHeight', Threshold, ...
    'MinPeakDistance', MinDist);

% Handle the case where there are no Rpeakslocs
if length(Rpeakslocs) < 2; Rpeakslocs = zeros(1,2); end
                                
% Calculate number of samples
samples = 1:length(ECG);
% Find the locations of the peak indices in seconds
t_Rpeaks = Rpeakslocs/Fs;

% Figure Subplot 1: ECG w/ markers at peaks
if Plot && corrected
    
    figure; subplot(3, 1, 1); hold on; plot(samples, ECG);
    plot(Rpeakslocs, ECG(Rpeakslocs), 'rv', 'MarkerFaceColor', 'r');

    % Format Figure 1
    grid on; xlabel('Sample #'); ylabel('ECG'); legend('ECG', 'R-peaks');
    
end

% Calculate RR interval and corresponding BPM for each interval
for i = 1:(length(Rpeakslocs) - 1)
    RR(i) = (Rpeakslocs(i + 1) - Rpeakslocs(i))*1000/Fs;	%in ms
    RR_bpm(i) = 60000./RR(i);                               %in bpm
end

% Calculate mean HR
RR_bpm_mean = mean(RR_bpm);

% Calculate HR stdev
RR_bpm_std = std(RR_bpm);

% Remove first datapoint from HR indices
t_HR = t_Rpeaks(2:end);

if corrected
    
    % Find outlier indices
    percntiles = prctile(RR_bpm,[5 95]); %5th and 95th percentile
    outlierIndex = RR_bpm < percntiles(1) | RR_bpm > percntiles(2);
    % Remove outlier values
    RR_bpm_corr = RR_bpm; RR_bpm_corr(outlierIndex) = [];
    t_Rpeaks_corr = t_Rpeaks; t_Rpeaks_corr(outlierIndex)=[];
    
    % Set corrected HR BPM
    RR_bpm_mean_corr = mean(RR_bpm_corr);

    % Set corrected HR stdev
    RR_bpm_std_corr = std(RR_bpm_corr);
    
    % Remove first datapoint from HR indices
    t_HR_corr = t_Rpeaks_corr(2:end);

end

if Plot && corrected
    
    % Figure 1 Subplot 2: RR interval vs. time
    subplot(3, 1, 2); plot(t_HR, RR_bpm, '.');
    figTitle = sprintf('HR=%.2f bpm, HR std=%.2f', RR_bpm_mean, RR_bpm_std);
    title(figTitle); xlabel('Time (s)'); ylabel('RR interval (bpm)');

    % Figure 1 Subplot 3: Corrected RR interval vs. time
    grid on; subplot(3, 1, 3); plot(t_HR_corr, RR_bpm_corr, '.');
    figTitle=sprintf(' corrected HR=%.2f bpm, HR std=%.2f', RR_bpm_mean_corr, RR_bpm_std_corr);
    title(figTitle); xlabel('Time (s)'); ylabel('RR interval corrected (bpm)'); grid on;
    
end

% Assign corrected variables to output
if corrected; t_HR = t_HR_corr; RR_bpm = RR_bpm_corr; end

end