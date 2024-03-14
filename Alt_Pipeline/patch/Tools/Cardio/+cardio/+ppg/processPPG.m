function pat = processPPG(signals, varargin)

% -------------------------------------------------------------------------
% SUMMARY:
% This function extracts the PAT from PPG signals. The PAT may be found
% using a variety of methods defined in cardio.general.pulseTime(). Some
% signals may requrie more complex analysis, in which case
% cardio.general.pulseTime() should be implemented directly
%
% Arguments (required):
% - signals     [NxM]       M vectors of PPG segments of length N
%
% Arguments (optional)
% - method                  Method for selecting pulse arrival
% - min                     Minimum value for X point
% - max                     Maximum value for X point
% - memory                  Numer of samples to store in memory
% - smooth                  Smoothing factor applied after differentiation
% - plot        FLAG        Plot results?
% -------------------------------------------------------------------------

% Extract optional arguments if necessary
if ~isempty(varargin)
    for arg = 1:length(varargin)
        switch varargin{arg}
            case 'method'; method = varargin{arg + 1};
            case 'min'; min = varargin{arg + 1};
            case 'max'; max = varargin{arg + 1};
            case 'memory'; memory = varargin{arg + 1};
            case 'smooth'; smooth = varargin{arg + 1};
            case 'plot'; Plot = varargin{arg + 1};
        end
    end
end

% Set optional arguments if they were not set by the user
if ~exist('method', 'var'); method = 'tangents'; end
if ~exist('min', 'var'); minVal = 1; end
if ~exist('max', 'var'); maxVal = size(signals, 1); end
if ~exist('memory', 'var'); memory = 30; end
if ~exist('smooth', 'var'); smoothing = 0; end
if ~exist('plot', 'var'); Plot = False; end

% -------------------------------------------------------------------------
% Compute PAT
% -------------------------------------------------------------------------

if smoothing == 0
    
    % Compute the PAT cardio.general.pulseTime()
    pat = cardio.general.pulseTime(signals, method, 'min', minVal, ...
        'max', maxVal, 'memory', memory);
    
else
    
    % Use the smoothing parameter
    pat = cardio.general.pulseTime(signals, method, 'min', minVal, ...
        'max', maxVal, 'memory', memory, 'smooth', smoothing);
    
end

% -------------------------------------------------------------------------
% Plot Results
% -------------------------------------------------------------------------

if Plot
    
    % Plot results as an imagesc overlaid with features
    title("PPG-derived PAT") figure; imagesc(normalize(signals)); 
    colormap(gray); hold on; plot(pat)
    
end

end

