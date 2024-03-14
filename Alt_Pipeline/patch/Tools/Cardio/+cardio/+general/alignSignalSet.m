function [ecg1_c, ecg2_c, set1_c, set2_c] = ...
    alignSignalSet(set1, set2, ECG1, ECG2, Fs1, Fs2, Plot)

% -------------------------------------------------------------------------
% This function aligns two separate sets of signals (Set 1 and Set 2) that
% are synced to two separate ECG signals and may have different sampling
% frequencies. The ECG signals are aligned using cross-correlation, and the
% alignment settings are then used to align the rest of the signals. The
% input arguments are as follows: 
% - set1:   % A cell array of all signals synced to ECG1
% - set2:   % A cell array of all signals synced to ECG2
% - ECG1:   % The ECG signal corresponding to Set 1
% - ECG2:   % The ECG signal corresponding to Set 2
% - Fs1:    % The sampling frequency of Set 1
% - Fs2:    % The sampling frequency of Set 2
% - Plot:   % Plot results?
% -------------------------------------------------------------------------

% Obtain number of signals in each set
numSignalsSet1 = length(set1);
numSignalsSet2 = length(set2);

% -------------------------------------------------------------------------
% Step 1: Resample the Higher-Frequency Signal Set
% -------------------------------------------------------------------------
% Initialize placeholders
set1_r = cell(size(set1));
set2_r = cell(size(set2));
Fs = 0;

% If Fs1 > Fs2, resample Set 1 at Fs2; else, resample Set 2 at Fs1
if Fs1 > Fs2
    % Set the global sampling frequency
    Fs = Fs2;
    % For each signal in the original set, resample at Fs2
    for i = 1:numSignalsSet1
        set1_r{i} = resample(set1{i}, Fs2, Fs1);
    end
elseif Fs2 > Fs1
    % Set the global sampling frequency
    Fs = Fs1;
    % For each signal in the original set, resample at Fs1
    for i = 1:numSignalsSet2
        set2_r{i} = resample(set2{i}, Fs1, Fs2);
    end
else
    % Handle the case where Fs1 = Fs2
    set1_r = set1; set2_r = set2;
end

% -------------------------------------------------------------------------
% Step 2: Align ECG1 and ECG2
% -------------------------------------------------------------------------
% Align these two signals using cross-correlation and obtain the offset
[ecg1_c, ecg2_c, lagSamples, ~] = cardio.general.alignSignals(ECG1 ,ECG2, Fs, Plot);

% -------------------------------------------------------------------------
% Step 3: Offset Signal Sets 1 and 2
% -------------------------------------------------------------------------
% Initialize placeholders
set1_c = cell(size(set1_r)); set2_c = cell(size(set2_r));

% If lagSamples < 0, Set 2 lags Set 1; thus, Set 1 will be truncated at the
% end and Set 2 will be truncated at the beginning by lagSamples. If
% lagSamples > 0, the opposite is true.
if lagSamples < 0
    for i = 1:numSignalsSet1
        set1_c{i} = set1_r{i}(1:end-abs(lagSamples));	% Truncate at end
    end
    for i = 1:numSignalsSet2
        set2_c{i} = set2_r{i}(abs(lagSamples):end);     % Truncate at beginning
    end
elseif lagSamples > 0
    for i = 1:numSignalsSet1
        set1_c{i} = set1_r{i}(lagSamples:end);          % Truncate at beginning
    end
    for i = 1:numSignalsSet2
        set2_c{i} = set2_r{i}(1:end-lagSamples);        % Truncate at end
    end
end

end