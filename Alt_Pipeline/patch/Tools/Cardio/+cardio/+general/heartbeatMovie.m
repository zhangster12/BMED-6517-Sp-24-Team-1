function heartbeatMovie(signals, varargin)

% -------------------------------------------------------------------------
% This function creates an animation of a heartbeat-separated signal with
% optional feature points overlaid.
%
% Arguments (required)
% - signals     [MxN]   N signal vectors of length M for visualization
%
% Arguments (optional)
% - 'features'  [NxF]   Feature points for visualization (F per heartbeat)
% - 'start'             Sample at which to start animation
% - 'stop'              Sample at which to stop animation
% - 'speed'             Animation speed ('slow', 'med', or 'fast')
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'features'); features = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'start'); start = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'stop'); stop = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'speed'); speed = varargin{arg + 1};
        end
    end
end

% Get the number of signals and features
numSignals = size(signals, 2);
if ~exist('features', 'var'); features = []; end
if ~isempty(features); numFeatures = size(features, 2); end

% Set defaults for optional arguments
if ~exist('start', 'var'); start = 1; end
if ~exist('stop', 'var'); stop = numSignals; end
if ~exist('speed', 'var'); speed = 'med'; end

% Set the speed parameter tau
switch speed
    case 'slow'
        tau = 1;
    case 'med'
        tau = 0.1;
    case 'fast'
        tau = 0.01;
    otherwise
        tau = 0.1;
end

% Plot the signals
figure
for i = start:stop
    
    % Plot the current signal
    plot(signals(:, i), 'k'); hold on; grid on;
    
    % Display figure title
    title("Segment " + string(i) + " of " + string(stop - start))
    
    % Plot the signal feature(s), if indicated
    if ~isempty(features) && isempty(find(isnan(features(i, :)), 1))
        for j = 1:numFeatures
            plot(features(i, j), signals(features(i, j), i), 'o')
        end
    end
    
    % Pause and reset figure
    pause(tau); clf('reset');
    
end

end

