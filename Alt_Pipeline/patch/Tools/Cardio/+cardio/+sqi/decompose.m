function S = decompose(signal, order, windows, Fs, varargin)

% -------------------------------------------------------------------------
% SUMMARY
% This function decomposes a signal using a windowed Fourier series, or
% short-time Fourier series (STFS). Each window is described as a linear
% combination of sinusoids of varying frequency and amplitude.
% 
% ARGUMENTS (REQUIRED)
% - signal    [Nx1]   Signal vector
% - order             Number of terms for Fourier series (max 8)
% - windows           Number of windows in which to partition signal
% - Fs                Sampling frequency (Hz)
% 
% ARGUMENTS (OPTIONAL)
% - period            Fix period to static value
% - offset            Fix offset to static value
% - limit             Limit coefficient values
% 
% RETURNS
% - S     Struct
%     - S.time        [Tx1]   Time vector for signal
%     - S.wintime     [Mx1]   Time vector for window
%     - S.coeffs      [Nx2]   Coefficients an and bn for each order N
%     - S.offset              Coefficient a0 in Fourier series
%     - S.period              Coefficient w in Fourier series
% 
% USAGE
% S = decompose(signal, order, windows, Fs);
% S = decompose(signal, order, windows, Fs, 'period', 50);
% -------------------------------------------------------------------------

% Parse optional input arguments
if nargin > 4
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'period'); period = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'offset'); offset = varargin{arg + 1}; 
        elseif strcmp(varargin{arg}, 'limit'); limit = varargin{arg + 1};
        end
    end
end

% Get length of signal and window
sigLen = length(signal); winLen = fix(sigLen/windows);

% Set time vector for signal
time = (0:(sigLen - 1))/Fs;

% Set placeholders for return values
S.wintime = cell(windows, 1); S.coeffs = cell(windows, 1);
S.offset = zeros(windows, 1); S.period = zeros(windows, 1);

% For each window, obtain the Fourier coefficients of the signal
initial = 1;    % Set initial point for windowing
for win = 1:windows
    
    % Window the signal
    sig = signal(initial:(initial + winLen - 1));   % Signal vector 
    t = time(initial:(initial + winLen - 1));       % Time vector
    
    % Use first iteration to set Fourier fit options if period provided
    if win == 1 && (exist('period', 'var') || exist('offset', 'var')) || exist('limit', 'var')
        
        % Obtain Fourier coefficients
        f = fit(t', sig, "fourier" + string(order));

        % Get fourier field options
        options = fitoptions("fourier" + string(order));
        % Get names of fields
        names = coeffnames(f); index = 1;
        % Find the field which contains the period/offset and update
        if ~exist('limit', 'var')
            upper = inf*ones(1, length(names)); lower = -inf*ones(1, length(names));
        else
            upper = limit*ones(1, length(names)); lower = -limit*ones(1, length(names));
        end
        if exist('period', 'var')   % If the period is specified
            % Get the index corresponding to the period
            for n = 1:length(names); if names{n}(1) == 'w'; index = n; break; end; end
            % Set upper and lower bounds for period
            upper(index) = period; lower(index) = period;
        end
        if exist('offset', 'var')   % If the offset is specified
            % Get the index corresponding to the offset
            for n = 1:length(names); if strcmp(names{n}, 'a0'); index = n; break; end; end
            % Set upper and lower bounds for offset
            upper(index) = offset; lower(index) = offset;
        end
        % Update Fourier coefficient parameters
        options.Upper = upper; options.Lower = lower;
        
    elseif win == 1
        
        % Se default options if no period provided
        options = fitoptions("fourier" + string(order));
        
    end
    
    % Obtain fit with specified options
    f = fit(t', sig, "fourier" + string(order), options);
    
    % Extract coefficients for return value
    coeffs = zeros(order, 2);	% Set placeholder for coefficients
    for ord = 1:order
        coeffs(ord, :) = [f.("a" + string(ord)) f.("b" + string(ord))];
    end
    
    % Organize segment information
    S.wintime{win} = t;     % Time vector for window
    S.coeffs{win} = coeffs; % Coefficients ax and bx
    S.offset(win) = f.a0;   % Coefficient a0 (offset)
    S.period(win) = f.w;    % Coefficient w (period)
    
    % Update initial point
    initial = initial + winLen;
    
end

% Set time vector for signal
S.time = time;

end

