function ICG = processICG(signals, varargin)

% -------------------------------------------------------------------------
% SUMMARY:
% This function extracts the B and X points from an ICG signal. The B point
% is identified as the maximum second derivative before the global maximum
% of the signal. The X point is identified as the signal minimum after the
% global maximum. Some signals require special configurations of the
% cardio.general.pulseTime() functions for consistent feature extraction;
% in these cases, the pulseTime function may be implemented directly.
%
% Arguments (required):
% - signals     [NxM]       M vectors of ICG segments of length N
%
% Arguments (optional)
% - bptMethod               Method for selecting B points
% - xptMethod               Method for selecting X points
% - bptMin                  Minimum value for B point
% - bptMax                  Maximum value for B point
% - xptMin                  Minimum value for X point
% - xptMax                  Maximum value for X point
% - memory                  Numer of samples to store in memory
% - smooth                  Smoothing factor applied after differentiation
% - plot        FLAG        Plot results?
%
% OUTPUT:
% - ICG             Struct
%   .bpt            Double vector of processed B point indices
%   .xpt            Double vector of processed X point indices
% -------------------------------------------------------------------------

% Extract optional arguments if necessary
if ~isempty(varargin)
    for arg = 1:length(varargin)
        switch varargin{arg}
            case 'bptMethod'; bptMethod = varargin{arg + 1};
            case 'xptMethod'; xptMethod = varargin{arg + 1};
            case 'bptMin'; bptMin = varargin{arg + 1};
            case 'bptMax'; bptMax = varargin{arg + 1};
            case 'xptMin'; xptMin = varargin{arg + 1};
            case 'xptMax'; xptMax = varargin{arg + 1};
            case 'memory'; memory = varargin{arg + 1};
            case 'smooth'; smoothing = vararging{arg + 1};
            case 'plot'; Plot = true;
        end
    end
end

% Set optional arguments if they were not set by the user
if ~exist('bptMethod', 'var'); bptMethod = 'diff2peak'; end
if ~exist('xptMethod', 'var'); xptMethod = 'globalMin'; end
if ~exist('bptMin', 'var'); bptMin = 1; end
if ~exist('xptMin', 'var'); xptMin = 1; end
if ~exist('bptMax', 'var'); bptMax = size(signals, 1); end
if ~exist('xptMax', 'var'); xptMax = size(signals, 1); end
if ~exist('memory', 'var'); memory = 30 ;end
if ~exist('smooth', 'var'); smoothing = 0; end
if ~exist('Plot', 'var'); Plot = false; end

% -------------------------------------------------------------------------
% Compute B-point
% -------------------------------------------------------------------------

if smoothing == 0
    
    % Compute the B-point index using cardio.general.pulseTime()
    ICG.bpt = cardio.general.pulseTime(signals, bptMethod, 'min', bptMin, ...
        'max', bptMax, 'memory', memory);
    
else
    
    % Use the smoothing parameter
    ICG.bpt = cardio.general.pulseTime(signals, bptMethod, 'min', bptMin, ...
        'max', bptMax, 'memory', memory, 'smooth', smoothing);
    
end

% -------------------------------------------------------------------------
% Compute X-point
% -------------------------------------------------------------------------

% Flipd the signals
flipped = flipud(signals);

if smoothing == 0
    
    % Compute the B-point index using cardio.general.pulseTime()
    xpt = cardio.general.pulseTime(flipped, xptMethod, 'min', xptMin, ...
        'max', xptMax, 'memory', memory); ICG.xpt = size(signals, 1) - xpt;
    
else
    
    % Use the smoothing parameter
    xpt = cardio.general.pulseTime(flipped, xptMethod, 'min', xptMin, ...
        'max', xptMax, 'memory', memory, 'smooth', smoothing);
    ICG.xpt = size(signals, 1) - xpt;
    
end

% -------------------------------------------------------------------------
% Plot Results
% -------------------------------------------------------------------------

if Plot
    
    % Plot results as an imagesc overlaid with features
    title("ICG B-point and X-point")
    figure; imagesc(normalize(signals)); colormap(gray); hold on
    plot(ICG.bpt); plot(ICG.xpt); legend("B-point", "X-point");
    
end

end

