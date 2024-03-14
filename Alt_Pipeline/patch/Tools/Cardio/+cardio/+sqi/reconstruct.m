function R = reconstruct(S, varargin)

% -------------------------------------------------------------------------
% SUMMARY
% This function reconstructs a signal from its windowed Fourier series, or
% short-time Fourier series (STFT). Each window is described as a linear
% combination of sinusoids of varying frequency and amplitude.
% 
% ARGUMENTS (REQUIRED)
% - S     Struct  Data structure for short-time Fourier series
%     - S.time    [Tx1]   Time vector for signal
%     - S.wintime [Mx1]   Time vector for window
%     - S.coeffs  [Nx2]   Coefficients an and bn for each order N
%     - S.offset          Coefficient a0 in Fourier series
%     - S.period          Coefficient w in Fourier series
% 
% ARGUMENTS (OPTIONAL)
% - 'noise'   Struct  Data structure for additional noise components
% 
% RETURNS
% - R     Struct
%     - R.signal  [Tx1]   Reconstructed signal
%     - R.time    [Tx1]   Time vector for reconstructed signal
% -------------------------------------------------------------------------

% Parse optional agruments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'noise'); D = varargin{arg + 1}; end
    end
end

% Obtain signal parameters
order = size(S.coeffs{1}, 1);           % Order of STFS
windows = length(S.coeffs);             % Number of windows
winLen = length(S.wintime{1});          % Window length

% Set placeholder for return value
R.time = S.time; R.signal = zeros(size(R.time));

% For each window, reconstruct the signal using the coefficients
initial = 1;    % Set initial point for windowing
for win = 1:windows
    
    % Set time vector
    t = S.wintime{win};
    
    % Set placeholder for segment result
    temp = S.offset(win); if exist('D', 'var'); temp = temp + D.offset(win); end
    
    % For each order, add a term to the above
    for ord = 1:order
        
        % Adjust coefficients if noise is provided
        if exist('D', 'var')
            S.coeffs{win}(ord, :) = S.coeffs{win}(ord, :) + D.coeffs{win}(ord, :);
        end
        
        % Generate signal component
        temp = temp + S.coeffs{win}(ord, 1)*cos(ord*t*S.period(win)) + ...
            S.coeffs{win}(ord, 2)*sin(ord*t*S.period(win));
    end
    
    % Append temp to signal
    R.signal(initial:(initial + winLen - 1)) = temp;
    
    % Update initial point
    initial = initial + winLen;
    
end

end