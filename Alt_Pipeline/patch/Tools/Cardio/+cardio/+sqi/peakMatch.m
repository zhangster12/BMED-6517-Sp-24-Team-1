function [distances, varargout] = peakMatch(signals, template, varargin)

% -------------------------------------------------------------------------
% SUMMARY
% This function implements a modified dynamic time warping (DTW) algorithm
% optimized for peak-matching in addition to waveform morphology matching.
% The function first extracts the signal and template features and uses it to
% construct a warping grid that determines anchor points in the warping path.
% The distance is then computed between signals along each warping path, the
% optimal path being the one that minimizes the final distance.
% 
% ARGUMENTS (REQ'D)
% - signals         [NxM]   Signal segments
% - template        [Nx1]   Template vector
% 
% ARGUMENTS (OPT'L)
% - 'maxDist'               Maximum distance between candidate features
%
% Outputs
% - distances       [Mx1]   Distance between each signal and the template
% - varargout{1}    {Mx1}   Cell vector of warped signals
% - varargout{2}    {Mx1}   Cell vector of warped templates
% - varargout{3}    {Mx1}   Cell vector of signal paths
% - varargout{4}    {Mx1}   Cell vector of template paths
% - varargout{5}    [Mx1]   Percent of matched peaks in template
% - varargout{6}    [Mx1]   Number of excess signal peaks
% -------------------------------------------------------------------------

% Parse optional arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'maxDist'); maxDist = varargin{arg + 1}; end
    end
end

% Set defaults
if ~exist('maxDist', 'var'); maxDist = 50; end

% Initialize placeholders for return values
matched = zeros(size(signals, 2), 1);   % Percent matched peaks in template
total = zeros(size(signals, 2), 1);     % Number of excess signal peaks
distances = zeros(size(signals, 2), 1); % Minimum matched distances (length-normalized)
warped_s = cell(size(signals, 2), 1);   % Cell vector of warped signals
warped_t = cell(size(signals, 2), 1);   % Cell vector of warped templates
sigPath = cell(size(signals, 2), 1);	% Cell vector of signal paths
temPath = cell(size(signals, 2), 1);	% Cell vector of template paths

% Perform analysis separately for each signal
for s = 1:size(signals, 2)
    
    % Set counters
    matchCounter = 0;
    
    % ---------------------------------------------------------------------
    % Extract Features
    % ---------------------------------------------------------------------
    
    % For the signal segments, get the peaks and valleys
    [peaks, valleys, ~] = cardio.general.getPeaks(signals(:, s));
    fs = [peaks; valleys]; [fs, idx] = sort(fs, 'ascend'); ts = cell(size(fs));
    % Determine labels for combined feature vector
    ts(1:length(peaks)) = {'p'}; ts(length(peaks) + 1:end) = {'v'}; ts = ts(idx);
    
    % For the template, get the peaks and valleys
    [peaks, valleys, ~] = cardio.general.getPeaks(template);
    ft = [peaks; valleys]; [ft, idx] = sort(ft, 'ascend'); tt = cell(size(ft));
    % Determine labels for combined feature vector
    tt(1:length(peaks)) = {'p'}; tt(length(peaks) + 1:end) = {'v'}; tt = tt(idx);
    
    % Set parameters for algorithm
    M = length(fs); N = length(ft);
    
    % ---------------------------------------------------------------------
    % Construct Grid
    % ---------------------------------------------------------------------
    
    % Initialize grid to 0
    A = zeros(M, N);
    
    % Where the labels match, set A to 1
    for n = 1:N; for m = 1:M; if strcmp(ts{m},tt{n}); A(m,n) = 1; end; end; end
    
    % ---------------------------------------------------------------------
    % Remove Outliers
    % ---------------------------------------------------------------------
    
    % Where the features are too far apart, set A to 0
    for n = 1:N; for m = 1:M; if abs(fs(m)-ft(n)) > maxDist; A(m,n) = 0; end; end; end
    
    % ---------------------------------------------------------------------
    % Obtain Losses for each Path
    % ---------------------------------------------------------------------
    
    % Set initial values
    idx_t = 1; idx_s = 1;       % Start indices
    path_t = []; path_s = [];   % Warp paths for signal and template
    
    % For each template feature...
    for n = 1:N
        
        % If there is a possible path...
        if sum(A(:, n)) > 0
            
            % Set placeholder for possible distances and warp paths
            dist = zeros(sum(A(:, n)), 1);
            wp_t = cell(sum(A(:, n)), 1); wp_s = cell(sum(A(:, n)), 1);
            coord = cell(sum(A(:, n)), 1);
            
            % Get indices
            indices = find(A(:, n) == 1);
            
            % For each possible path
            for p = 1:length(indices)
                
                % Set signal feature
                m = indices(p);
                
                % Segment signal
                tmp = template(idx_t:ft(n)); sig = signals(idx_s:fs(m), s);
                
                % Warp signals
                [dist(p), wp_t{p}, wp_s{p}] = dtw(tmp, sig);
                
                % Correct distance for length of warping path
                dist(p) = dist(p)/length(wp_t{p});
                
                % Save coordinates and update counter
                coord{p} = [m, n];
                
            end
            
            % Select path with shortest distance
            [~, op_path] = min(dist);
            
            % Remove invalid paths from grid
            A(:, 1:n) = 0; A(1:coord{op_path}(1), :) = 0;
            
            % Save path
            path_t = [path_t; (idx_t - 1) + wp_t{op_path}(:)];
            path_s = [path_s; (idx_s - 1) + wp_s{op_path}(:)];
            
            % Update start indices
            idx_t = ft(n) + 1; idx_s = fs(coord{op_path}(1)) + 1;
            
            % Indicate matched peak
            matchCounter = matchCounter + 1;
            
        end
        
    end
    
    % Get the euclidian distances of the warped signals
    minDist = norm(signals(path_s, s) - template(path_t));
    
    % Update counters
    matched(s) = matchCounter; total(s) = M; distances(s) = minDist/length(path_t);

    % Capture warped signals
    warped_s{s} = signals(path_s, s); warped_t{s} = template(path_t);
    sigPath{s} = path_s; temPath{s} = path_t;
    
%     % Plot orignal signal, template, and points
%     figure; hold on; grid on;
%     plot(template); for f = ft; plot(f, template(f), 'ob'); end
%     plot(signals(:, s)); for f = fs; plot(f, signals(f, s), 'or'); end
%     
%     % Produce warped signal and template
%     warped_t = template(path_t); warped_s = signals(path_s, s);
%     figure; hold on; grid on;
%     plot(warped_t)
%     for f = ft
%         newf = cardio.general.map(f, path_t);
%         plot(newf, warped_t(newf), 'ob')
%     end
%     plot(warped_s)
%     for f = fs
%         newf = cardio.general.map(f, path_s);
%         plot(newf, warped_s(newf), 'or')
%     end
%     
%     % Plot the DTW signals
%     [~, ix, iy] = dtw(signals(:, s), template);
%     figure; hold on; grid on;
%     warped_s = signals(ix, s); warped_t = template(iy);
%     plot(warped_t)
%     for f = ft
%         newf = cardio.general.map(f, iy);
%         plot(newf, warped_t(newf), 'ob')
%     end
%     plot(warped_s)
%     for f = fs
%         newf = cardio.general.map(f, ix);
%         plot(newf, warped_s(newf), 'or')
%     end

end

% Return output arguments
if nargout > 1; varargout{1} = warped_s; end
if nargout > 2; varargout{2} = warped_t; end
if nargout > 3; varargout{3} = sigPath; end
if nargout > 4; varargout{4} = temPath; end
if nargout > 5; varargout{5} = matched; end
if nargout > 6; varargout{6} = total; end

end