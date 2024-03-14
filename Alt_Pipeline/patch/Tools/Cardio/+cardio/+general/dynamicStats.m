function [statVector, error, varargout] = ...
    dynamicStats(signal, order, Plot, varargin)

% -------------------------------------------------------------------------
% This function returns vectors of statistics adjusted for trends in the
% signal. This starts by fitting the signal to a polynomial of the
% specified order, then returning the values of that polynomial over the
% input range. Further processing may use these values as a basis for
% future resampling decisions, when determining outliers. Input parameters
% of the function include:
% - signal: [Nx1]   Signal vector
% - order:  Int     Order of polynomial to which data will be fit
% - Plot:   Bool    Plot results?
% - indices [Nx1]   (OPTIONAL) Indices on which to fit signal
%
% Returns:
% - statsVector:    [N x 1] Expected value at each sample
% - error:          Estimate of the standard deviation of the error in
%                   predicting a future observation at X by P(X)
% - varargout
%   - polySig       Polynomial coefficients
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Process Signal
% -------------------------------------------------------------------------

% Generate X and Y coordinates
X = (1:length(signal))';            % Set X coordinate as indices
Y = fillmissing(signal, 'linear');  % Fill missing values in Y coordinate
len = length(X);                    % Original signal length

% Remove unnecessary indices
if nargin > 3
    
    % Get vector of indices to keep
    indices = varargin{1};
    
    extraIndices = ~ismember(X, indices);           % Indices to remove
    X(extraIndices) = []; Y(extraIndices) = [];     % Remove indices
    
end

% -------------------------------------------------------------------------
% Get Polynomial Fit
% -------------------------------------------------------------------------

% Fit signal to polynomial
[polySig, S] = polyfit(X, Y, order);

% Return polynomial fit and error
[statVector, error] = polyval(polySig, 1:len, S);

% Assign return value
varargout{1} = polySig;

% -------------------------------------------------------------------------
% Plot results
% -------------------------------------------------------------------------

if Plot
    
    figure; hold on; grid on;               % Initialize figure
    plot(X, signal); plot(X, statVector, '--k');  % Plot signal and best fit
    title("Generating Best-Fit for Signal (Delta: " + string(error) + ")")
    xlabel('Timestep'); ylabel('Expected Value');
    
end

end

