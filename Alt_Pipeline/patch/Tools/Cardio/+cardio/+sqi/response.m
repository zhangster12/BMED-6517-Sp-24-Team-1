function [dist, dSignal, d] = response(signal, D, Fs, type)

% -------------------------------------------------------------------------
% SUMMARY
% This function computes the distance between a template signal and the same
% signal with a sinusoidal disturbance.
% 
% ARGUMENTS
% - signal    [Nx1]   Signal vector
% - D         Struct  Disturbance
%     - a             Coefficient of cosine sinusoidal disturbance
%     - b             Coefficient of sine sinusoidal disturbance
%     - period        Period of sinusoidal disturbance
%     - order         Order or sinusoidal disturbance (harmonic)
%     - window        Window number of disturbance
%     - numWindows    Number of windows in signal
%     - offset        Offset of signal disturbance
% - Fs                Sampling frequency (Hz)
% - type              Distance type
%     - "norm"        L2-norm
%     - "dtw"         Dynamic time warping
% -------------------------------------------------------------------------
    
% Generate time vector
time = (0:(length(signal) - 1))/Fs;

% Generate sinusoidal disturbance
d = zeros(size(signal));                        % Set placeholder
winLength = fix(length(signal)/D.numWindows);   % Length of each window

% Generate sinusoid to place in appropriate window
t = time((1 + (winLength*(D.window - 1))):winLength*D.window);
offset = D.offset*ones(size(t));
s = offset + D.a*cos(D.order*t*D.period) + D.b*sin(D.order*t*D.period);

% Place sinusoid in desired window
d((1 + (winLength*(D.window - 1))):winLength*D.window) = s;

% Get distance between signal and signal + disturbance
dSignal = signal + d;
if strcmp(type, "norm"); dist = norm(signal - dSignal); 
else; dist = dtw(signal, dSignal);
end

end

