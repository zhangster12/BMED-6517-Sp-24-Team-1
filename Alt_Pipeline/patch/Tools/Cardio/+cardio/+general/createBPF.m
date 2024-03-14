function Hd = createBPF(Fstop1, Fstop2, Fs, varargin)

% -------------------------------------------------------------------------
% This function creates a band pass filter with the given input parameters:
% Fs:       Sampling frequency
% Fpass1:   First pass band frequency
% Fstop1:   First stop band frequency
% Fpass2:   Second pass band frequency
% Fstop2:   Second stop band frequency
%
% ARGS (OPT'L)
% 'kaiser'  FLAG    Kaiser Window
% 'order'           Order of Kaiser window
% 'Fpass1'          First pass band frequency
% 'Fpass2'          Second pass band frequency
% -------------------------------------------------------------------------

% Parse optional arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'kaiser'); Kaiser = true;
        elseif strcmp(varargin{arg}, 'order'); order = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'Fpass1'); Fpass1 = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'Fpass2'); Fpass2 = varargin{arg + 1};
        end
    end
end

% Set defaults
if ~exist('Kaiser', 'var'); Kaiser = false; end
if ~exist('order', 'var'); order = 0; end
if ~exist('Fpass1', 'var'); Fpass1 = Fstop1; end
if ~exist('Fpass2', 'var'); Fpass2 = Fstop2; end

% Set other parameters
Dstop1 = 0.001;           % First Stopband Attenuation
Dpass  = 0.05;            % Passband Ripple
Dstop2 = 0.0001;          % Second Stopband Attenuation
dens   = 20;              % Density Factor

if ~Kaiser
    % Calculate the order from the parameters using FIRPMORD.
    [N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
                              0], [Dstop1 Dpass Dstop2]);

    % Calculate the coefficients using the FIRPM function.
    b  = firpm(N, Fo, Ao, W, {dens});
    Hd = dfilt.dffir(b);
elseif order ~= 0
    flag = 'scale';  % Sampling Flag
    Beta = 0.5;      % Window Parameter
    % Create the window vector for the design algorithm.
    win = kaiser(order+1, Beta);

    % Calculate the coefficients using the FIR1 function.
    b  = fir1(order, [Fstop1 Fstop2]/(Fs/2), 'bandpass', win, flag);
    Hd = dfilt.dffir(b);
else
    flag = 'scale';  % Sampling Flag
    % Calculate the order from the parameters using KAISERORD.
    [N,Wn,BETA,TYPE] = kaiserord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 ...
                                 1 0], [Dstop1 Dpass Dstop2]);

    % Calculate the coefficients using the FIR1 function.
    b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
    Hd = dfilt.dffir(b);
end

end

