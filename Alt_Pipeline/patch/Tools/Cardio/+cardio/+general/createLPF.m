function Hd = createLPF(Fpass, Fstop, Fs)

% -------------------------------------------------------------------------
% This function creates a low pass filter with the given input parameters:
% Fs:       Sampling frequency
% Fpass1:   Pass band frequency
% Fstop1:   Stop band frequency
% -------------------------------------------------------------------------

% Equiripple Lowpass filter designed using the FIRPM function.

Dpass = 0.1;    % Passband Ripple
Dstop = 0.001;  % Stopband Attenuation
dens  = 20;     % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);

% [EOF]
