function [signal, theta_e] = synthetic(template, Fs, order, windows, varargin)

% -------------------------------------------------------------------------
% SUMMARY
% This function generates a synthetic SCG signal based on a template. The
% template is first projected in to the subspace defined by the Fourier
% series. Subsequently up to N sinusoids of order up to 'order' are
% generated, with parameters drawn from a bounded distribution.
% 
% ARGUMENTS (REQ'D)
% - template  [Nx1]   Prototypical signal vector in Fourier domain
% - Fs                Sampling frequency
% - order             Maximum order of noise components
% - windows           Number of windows in which to divide signal segments
% 
% ARGUMENTS (OPTIONAL)
% - varargin
%     - 'N'           Maximum number of noise components
%     - 'lower'       Lower bound on noise amplitude
%     - 'upper'       Upper bound on noise amplitude
%     - 'period'      Period of Fourier series
% -------------------------------------------------------------------------

% Parse inputs
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'lower'); lower = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'upper'); upper = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'period'); period = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'N'); N = varargin{arg + 1};
        end
    end
end

% Set defaults
if ~exist('lower', 'var'); lower = -10; end
if ~exist('upper', 'var'); upper = 10; end
if ~exist('period', 'var'); period = floor(length(template)/(2*windows)); end
if ~exist('N', 'var'); N = 10; end

% Generate time vector
time = (0:(length(template) - 1))/Fs;

% Get the length of each window
winLength = floor(length(template)/windows);

% Determine the number of noise components
N = randi([1, N]);

% For each noise component, determine the order
R = randi([1, order], N, 1);

% For each noise component, determine the window
W = randi([1, windows], N, 1);

% For each noise component, determine the coefficients
C = zeros(N, 2);
for n = 1:N
    C(n, :) = (upper - lower)*rand(2, 1) + lower;
end

% Initialize the synthetic signal vector as the template
signal = template;

% Generate the noisy signal and add it to the synthetic signal
for n = 1:N
    
    % Generate noise and place the noise in the proper window
    noise = zeros(size(signal));
    t = time((1 + (winLength*(W(n) - 1))):winLength*W(n));
    temp = C(n, 1)*cos(R(n)*t*period) + C(n, 2)*sin(R(n)*t*period);
    noise((1 + (winLength*(W(n) - 1))):winLength*W(n)) = temp;
    
    % Add the noise to the signal
    signal = signal + noise;
    
end

% Organize theta
% theta = zeros(2, order, windows);
% for n = 1:N; theta(:, R(n), W(n)) = C(n, :); end
theta = zeros(order, windows);
for n = 1:N; theta(R(n), W(n)) = 1;

% Return theta_e
% theta_e = {a1_w1, b1_w1, a2_w1, b2_w1 ... aO_wW, b0_wW]
theta_e = theta(:);

end

