function Hd = createHPF(varargin)

% -------------------------------------------------------------------------
% This function creates a high pass filter with Kaiser window.
%
% Arguments (optional)
% - 'Fs'    Sampling Frequency (Hz)
% - 'Fc'    Cutoff Frequency (Hz)
% - 'order' Filter order (Hz)
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'Fs'); Fs = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'Fc'); Fc = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'order'); order = varargin{arg + 1};
        end
    end
end

% Set default values for optional inputs
if ~exist('Fs', 'var'); Fs = 2000; end
if ~exist('Fc', 'var'); Fc = 50; end
if ~exist('order', 'var'); order = 50; end

flag = 'scale';  % Sampling Flag
Beta = 0.5;      % Window Parameter

% Create the window vector for the design algorithm.
win = kaiser(order + 1, Beta);

% Calculate the coefficients using the FIR1 function.
b  = fir1(order, Fc/(Fs/2), 'high', win, flag);
Hd = dfilt.dffir(b);

% [EOF]
