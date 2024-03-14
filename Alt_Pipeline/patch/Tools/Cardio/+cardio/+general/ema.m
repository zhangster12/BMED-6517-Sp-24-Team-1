function signalOutput = ema(signal, M, plotResults)

% -------------------------------------------------------------------------
% EXPONENTIAL MOVING AVERAGE
%
% This function takes the exponential moving average of a signal with a
% specified smoothing factor M. The signal should be an nxm matrix where m
% is the number of signal slices and n is the length of each slice. This
% function accepts the following arguments:
% - signal:         Signal slices with size [ n x m ]
% - M:              Smoothing factor (1 = no smoothing, > 1 = more smoothing)
% - plotResults:    BOOL indicating whether to plot results
% -------------------------------------------------------------------------

% Initialize matrices to hold exponential moving average of signal amplitude
signalInput = zeros(size(signal));
signalOutput = zeros(size(signal));

% Set constant for moving average calculation
alpha = 2/(M+1);

% For each signal segment...
for i = 1:size(signal,2)
    
    % Set placeholder for ICG input
    signalInput(:,i) = signal(:,i);
    
    % Map the first input directly to the output
    % For all other inputs, the output is a*input(t) + (1-a)output(t-1)
    if i==1
        signalOutput(:,i)=signalInput(:,i);
    else
        signalOutput(:,i)=alpha*signalInput(:,i)+(1-alpha)*signalOutput(:,i-1);
    end
    
end

if plotResults
    % Figure Subplot 1: Plot averaged signal
    figure
    subplot(211);
    hold on;
    for i=1:size(signalOutput,2)
        plot(signalOutput(:,i));
    end
    % Format figure
    xlabel('Timestamp');
    title('EMA Signal');
    ylabel('EMA Signal Amplitude');
    grid on;

    % Figure Subplot 2: Plot raw signal (not averaged)
    subplot(212);
    hold on;
    for i=1:size(signal,2)
        plot(signal(:,i));
    end
    % Format figure
    xlabel('Timestamp');
    title('Raw Signal Amplitude'); 
    ylabel('Raw Signal');
    grid on;
end

end