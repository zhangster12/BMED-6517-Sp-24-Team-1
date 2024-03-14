function SCG = ...
    processSCG(signalSegments, numFeatures, M, tol, Plot, varargin)

% -------------------------------------------------------------------------
% SUMMARY:
% This function extracts the time-domain features of an SCG signal with
% potential transient corruptions and discontinuities. To remove
% discontinuities, the following steps are performed. Once a feature is
% extracted from the signal using the traditional peak-counting method,
% gaussian clustering is used to separate signal segments for which there
% are discontinuities. A line of best fit is generated for each segment,
% and the signal is resampled such that the chosen peaks are closest to the
% line of best fit. The line which minimizes the detrended standard
% deviation of the final signal is determined as the line of best fit, and
% the features are reselected along this line. Subsequent features are
% chosen to follow the indices of the preceding features. If there still
% exists low signal quality due to highly complex transient disruptions,
% the genetic algorithm library (Smart.genetic) may be used to properly
% sort the signal vectors.
%
% ARGUMENTS (MANDATORY):
% - signalSegments  [NxM]   M vectors of SCG segments of length N
% - numFeatures             Number of time-domain features to extract
%                           - 1: first minimum
%                           - 2: first maximum
%                           - 3: second minimum ...
%                           - 6: third maximum
% - M                       Smoothing factor for exponential moving average
% - tol                     Tolerance for outliers (in standard deviations)
% - Plot            Bool    Plot results?
%
% ARGUMENTS (OPTIONAL)
% - K                       Maximum cluster size for mixed gaussian model
% - iterations              Number of iterations for mixed gaussian model
% - order                   Order number for lines of best fit
% - ao              [Nx1]   AO points for comparison
% - removed         {2x1}   Vectors of removed indices for AO point and signal
%                   {1}:    Removed indices for AO point
%                   {2}:    Removed indices for SCG signal
% - rawOnly         Bool    Raw indices only?
% - inflect         FLAG    Consider inflection points for resampling
% - startModel      {GMM}   Starting GMM for analysis for each feature
% - min             FLAG    Select the first local minimum as the first feature (default)
% - max             FLAG    Select the first local maximum as the first feature
%
% OUTPUT:
% - SCG             Struct
%   .features       {numFeatures}   Cell array of processed feature vectors
%   .startModel     {GMM}           Optimal GMM parameters for each feature
%
% USAGE:
% scg = cardio.scg.processSCG(signalSegments, numFeatures, M, tol, Plot)
% scg = cardio.scg.processSCG(signalSegments, numFeatures, M, tol, Plot, 'K', 3)
% scg = cardio.scg.processSCG(signalSegments, numFeatures, M, tol, Plot, 'iterations', 100)
% scg = cardio.scg.processSCG(signalSegments, numFeatures, M, tol, Plot, 'order', 5)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Part 0: Setup
% -------------------------------------------------------------------------

% Extract optional arguments (if necessary)
if nargin > 5
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'K'); K = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'iterations'); iterations = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'order'); order = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'ao'); ao = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'removed'); removed = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'rawOnly'); rawOnly = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'inflect'); inflect = true;
        elseif strcmp(varargin{arg}, 'startModel'); startModel = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'min'); firstMin = true;
        elseif strcmp(varargin{arg}, 'max'); firstMin = false;
        end
    end
end

% Set variables to defaults if not set by user
if ~exist('K', 'var'); K = 2; end
if ~exist('iterations', 'var'); iterations = 100; end
if ~exist('order', 'var'); order = 1; end
if ~exist('rawOnly', 'var'); rawOnly = false; end
if ~exist('inflect', 'var'); inflect = false; end
if ~exist('firstMin', 'var'); firstMin = true; end

% Initialize global arrays
indices = cell(numFeatures, 1);         % Raw indices derived from peak-counting (overwritten)
origIdx = cell(numFeatures, 1);         % Raw indices derived from peak-counting (not overwritten)
newIdx = cell(numFeatures, K, order);	% Indices derived from resampling clusters
baselineIdx = cell(numFeatures, order); % Indices derived from resampling baseline
resIdx = cell(numFeatures, 1);          % Indices derived from optimal resampling
finalIdx = cell(numFeatures, 1);        % Indices with outliers removed
remIdx = cell(numFeatures, 1);          % Removed indices

% -------------------------------------------------------------------------
% Part 1: Exponential Moving Average
% -------------------------------------------------------------------------

% Arguments: signal, M, plotResults
scg = cardio.general.ema(signalSegments, M, false);

% -------------------------------------------------------------------------
% Part 2: Extract Features
% -------------------------------------------------------------------------

% Set starting index to search for features
startIdx = ones(size(scg, 2), 1);

for feature = 1:numFeatures
    
    % ---------------------------------------------------------------------
    % Part 2(a): Extract Features with Peak-Counting
    % ---------------------------------------------------------------------
    
    if firstMin
        % Even features correspond to maxima -> odd features correspond to minima
        if mod(feature,2) == 0; findMax = true; else; findMax = false; end
    else
        % Even features correspond to minima -> odd features correspond to maxima
        if mod(feature,2) == 0; findMax = false; else; findMax = true; end
    end
    
    % Set placeholder for feature indices
    indices{feature} = zeros(size(scg, 2), 1);
    
    % Set placeholders for peaks and valleys in each segment
    peaks = cell(size(scg, 2), 1); valleys = cell(size(scg, 2), 1);
    % Inflection points (resampling)
    inflections = cell(size(scg, 2), 1);
    
    % Start index must  be >= 0
    startIdx(startIdx < 1) = 1;
    
    % Find desired feature type along each signal segment
    for segment = 1:size(scg, 2)
        
        % Extract segment
        if startIdx(segment) <= 0 % || (startIdx(segment)>length(scg(:,segment))))
            disp('error')
        end; sig = scg(startIdx(segment):end, segment);
       
        % Find all features of the given type
        [peaks{segment}, valleys{segment}, inflections{segment}] = ...
            cardio.general.getPeaks(sig);
        
        % Add back start indices
        peaks{segment} = peaks{segment} + startIdx(segment);
        valleys{segment} = valleys{segment} + startIdx(segment);
        
        % Replace empty vectors with NaN
        if isempty(peaks{segment}); peaks{segment} = NaN; end
        if isempty(valleys{segment}); valleys{segment} = NaN; end
        
        % Return the first feature within the search range
        if findMax; indices{feature}(segment) = min(peaks{segment}); else; ...
                indices{feature}(segment) = min(valleys{segment}); end
        
        % Save indices
        origIdx{feature}(segment) = indices{feature}(segment);
        
        % Get raw features
        % Find all features of the given type without correcting for segment
        [rawPeaks, rawValleys, ~] = ...
            cardio.general.getPeaks(scg(:,segment));
        % Write appropriate feature to output value
        if feature == 1
            if(isempty(rawValleys) && exist('SCG','var'))
                SCG.raw{feature}(segment) = SCG.raw{feature}(segment-1); 
            elseif(exist('SCG','var'))
                SCG.raw{feature}(segment) = rawValleys(1);
            else
                SCG.raw{feature}(segment) = 100;    
            end
        elseif feature == 2
            SCG.raw{feature}(segment) = rawPeaks(1);
        elseif feature == 3
            SCG.raw{feature}(segment) = rawValleys(2);
        else
            SCG.raw{feature}(segment) = rawPeaks(2);
        end
        
        
    end
    
    if ~rawOnly
    
    % Set initial performance
    PERFORMANCE = 0;
    
    % ---------------------------------------------------------------------
    % Iterate Until Best-Fit Converges
    % ---------------------------------------------------------------------
    while true
    
    % ---------------------------------------------------------------------
    % Part 2(b): Cluster Features with Gaussian Mixture Models (GMMs)
    % ---------------------------------------------------------------------
    
    % Perform GMM Clustering
    % Arguments: signal, Plot, (maxClusters)/(startModel, numClusters), (iterations)
    if exist('startModel', 'var')
        gmm = smart.cluster.gaussian(indices{feature}, Plot, 'startModel', startModel{feature}, ...
            'numClusters', startModel{feature}.numClusters, 'iterations', iterations);
        clusters = gmm.indices; SCG.startModel{feature} = gmm.model;
    else
        gmm = smart.cluster.gaussian(indices{feature}, Plot, 'maxClusters', K, ...
            'iterations', iterations);
        clusters = gmm.indices; SCG.startModel{feature} = gmm.model;
    end
    % Handle case where there are no outliers (Cluster 1)
    % In either case, shift cluster numbers down
    if ~isempty(find(clusters == 1, 1))
        numClusters = length(unique(clusters)) - 1; clusters = clusters - 1;
    else
        numClusters = length(unique(clusters)); clusters = clusters - 1;
    end
    
    % ---------------------------------------------------------------------
    % Part 2(c): Determine Optimal Line of Best Fit
    % ---------------------------------------------------------------------
    
    % Perform analysis for each order of best fit line
    for ord = 1:order
    
    % Perform analysis for each cluster
    for cluster = 1:max(1,numClusters)
        
        % Get indices of features in cluster (cluster 1 is for outliers only)
        allIndices = 1:length(indices{feature}); featIdx = allIndices(clusters == cluster);
       
        % Arguments: signal, order, Plot, varargin (indices)
        [E, error] = cardio.general.dynamicStats(indices{feature}, ord, false, featIdx);
        
        % Set placeholders for recording updated features
        newIdx{feature, cluster, ord} = zeros(size(indices{feature}));
        newIdx{feature, cluster, ord}(:) = NaN;
        
        % Resample features
        % Arguments: originalIndex, peakArray, infArray, expectedValue
        for segment = 1:length(indices{feature})
            if ~isnan(indices{feature}(segment)) % Resample iff the index is not NaN
                if findMax
                    if inflect
                        newIdx{feature, cluster, ord}(segment) = ...
                            cardio.scg.resamplePeaks(indices{feature}(segment), ...
                            [peaks{segment}; inflections{segment}], [], E(segment));
                    else
                        newIdx{feature, cluster, ord}(segment) = ...
                            cardio.scg.resamplePeaks(indices{feature}(segment), ...
                            peaks{segment}, [], E(segment));
                    end
                else
                    if inflect
                        newIdx{feature, cluster, ord}(segment) = ...
                            cardio.scg.resamplePeaks(indices{feature}(segment), ...
                            [valleys{segment}; inflections{segment}], [], E(segment));
                    else
                        newIdx{feature, cluster, ord}(segment) = ...
                            cardio.scg.resamplePeaks(indices{feature}(segment), ...
                            valleys{segment}, [], E(segment));
                    end
                end
            end
        end
        
        % Record best fit line
        bestFit.line{feature, cluster, ord} = E; bestFit.error{feature, cluster, ord} = error;

        % Get the performance as the correlation coefficient between the SCG
        % and AO gold standard signals (if provided), else the inverse of the 
        % detrended standard deviation of the signal
        [~, error] = cardio.general.dynamicStats(newIdx{feature, cluster, ord}, 1, false);
        bestFit.performance(feature, cluster, ord) = 1/max(error);
        if exist('ao', 'var') && exist('removed', 'var')
            bestFit.performance(feature, cluster, ord) = cardio.general.rsquared(ao, ...
                newIdx{feature, cluster, ord}, removed{1}, removed{2})/mean(error);
        elseif exist('ao', 'var')
            bestFit.performance(feature, cluster, ord) = cardio.general.rsquared(ao, newIdx{feature, cluster, ord})/mean(error);
        end
        
    end
    
    % The segments obtained during clustering must outperform
    % resampling against a line of best fit (guard against overfitting)
    % Arguments: signal, order, Plot, varargin (indices)
    [E, error] = cardio.general.dynamicStats(indices{feature}, ord, false);
    
    % Set placeholders for recording updated features
    baselineIdx{feature, ord} = zeros(size(indices{feature}));
    baselineIdx{feature, ord}(:) = NaN;
    
    % Resample features
    % Arguments: originalIndex, peakArray, infArray, expectedValue
    for segment = 1:length(indices{feature})
        if ~isnan(indices{feature}(segment)) % Resample iff the index is not NaN
            if findMax
                if inflect
                    baselineIdx{feature, ord}(segment) = cardio.scg.resamplePeaks(indices{feature}(segment), ...
                        [peaks{segment}; inflections{segment}], [], E(segment));
                else
                    baselineIdx{feature, ord}(segment) = cardio.scg.resamplePeaks(indices{feature}(segment), ...
                        peaks{segment}, [], E(segment));
                end
            else
                if inflect
                    baselineIdx{feature, ord}(segment) = cardio.scg.resamplePeaks(indices{feature}(segment), ...
                        [valleys{segment}; inflections{segment}], [], E(segment));
                else
                    baselineIdx{feature, ord}(segment) = cardio.scg.resamplePeaks(indices{feature}(segment), ...
                        valleys{segment}, [], E(segment));
                end
            end
        end
    end
    
    % Record best fit line
    baseline.line{ord} = E; baseline.error{ord} = error;
    
    % Get the performance as the correlation coefficient between the SCG
    % and AO gold standard signals (if provided), else the inverse of the 
    % detrended standard deviation of the signal
    [~, error] = cardio.general.dynamicStats(baselineIdx{feature, ord}, 1, false);
    baseline.performance(ord) = 1/max(error);
    if exist('ao', 'var') && exist('removed', 'var')
        baseline.performance(ord) = cardio.general.rsquared(ao, ...
            baselineIdx{feature, ord}, removed{1}, removed{2})/mean(error);
    elseif exist('ao', 'var')
        baseline.performance(ord) = cardio.general.rsquared(ao, baselineIdx{feature, ord})/mean(error);
    end
    
    end
    
    % Get best fit line (clusters must exceed performance of best fit)
    % 1. Find the baseline fit with minimal error
    [basePerformance, baseOrder] = max(baseline.performance);
    % 2. Find the resampled cluster and order with minimal error
    temp = squeeze(bestFit.performance(feature, :, :));
    % (Fix random error)
    if size(feature, 1) == 1 && order == 1; temp = temp(:); end
    [performance, clusterOrder] = max(temp(:));
    [bestCluster, bestOrder] = ind2sub(size(temp), clusterOrder);
    % 3. Select the best fit
    if performance > basePerformance
        bestFitLine = bestFit.line{feature, bestCluster, bestOrder};
        bestFitError = bestFit.error{feature, bestCluster, bestOrder};
    else
        bestFitLine = baseline.line{baseOrder};
        bestFitError = baseline.error{baseOrder};
    end
    
    % ---------------------------------------------------------------------
    % Part 2(d): Resample Features with Optimal Line
    % ---------------------------------------------------------------------
    
    % Set placeholder for resampled indices
    % disp("Feature: " + string(feature)); disp("---");
    if performance > basePerformance
        resIdx{feature} = newIdx{feature, bestCluster, bestOrder};
        % disp("Selected Cluster " + string(bestCluster) + " with Order " + string(bestOrder));
    else; resIdx{feature} = baselineIdx{feature, baseOrder};
        % disp("Selected Baseline Cluster with Order " + string(baseOrder));
    end
    
    % If the performance has stopped increasing, break the loop now
    if performance <= PERFORMANCE && basePerformance <= PERFORMANCE; break
    else
        % Else, update the performance
        if performance > basePerformance; PERFORMANCE = performance;
        else; PERFORMANCE = basePerformance;
        end
    end
    
    % Update current indices
    indices{feature} = resIdx{feature};
    
    end
    
    % ---------------------------------------------------------------------
    % Part 3: Remove Outliers
    % ---------------------------------------------------------------------
    
    % Set placeholder for return values
    finalIdx{feature} = resIdx{feature}; remIdx{feature} = [];
    
    % Replace an index with NaN if it lies outside the tolerance
    for segment = 1:length(finalIdx{feature})
       
        % Calculate tolerance zone
        minVal = bestFitLine(segment) - tol*mean(bestFitError);   % Min value
        maxVal = bestFitLine(segment) + tol*mean(bestFitError);   % Max value
        
        % Replace outliers
        if finalIdx{feature}(segment) < minVal || finalIdx{feature}(segment) > maxVal
            finalIdx{feature}(segment) = NaN;
            remIdx{feature} = [remIdx{feature} segment];
        end
        
    end
    
    % Write finalized values to output array
    SCG.features{feature} = finalIdx{feature};
    
    % Update start index for next feature (fill NaNs with nearest neighbor)
    recon = fillmissing(finalIdx{feature},'nearest');
    startIdx = round(cardio.general.dynamicStats(recon, 1, false));
    
    else
        
        % Set empty placeholder for return value
        SCG.features = [];
        
    end
    
end

% -------------------------------------------------------------------------
% Part 4: Calculate Statistics
% -------------------------------------------------------------------------

if ~rawOnly
% Initialize placeholder for statistics
statistics = zeros(numFeatures, 2); % [Mean] [Standard Deviation]
norms = cell((numFeatures), 1);     % Cell array for norm PDFs
x = 0:size(SCG.features{1}, 1);     % Define x-coordinate
% Return statistics for each feature
for feature = 1:numFeatures
    statistics(feature, 1) = nanmean(SCG.features{feature});    % Mean
    statistics(feature, 2) = nanstd(SCG.features{feature});     % Stddev
    % Fit norm PDF
    norms{feature} = normpdf(x, statistics(feature,1), statistics(feature,2));
end
end

% -------------------------------------------------------------------------
% Part 5: Visualization
% -------------------------------------------------------------------------

if Plot && ~rawOnly
    
    % Set graph colors
    colors = {[0, 0, 0], [1, 0, 0], [0, 0.4470, 0.7410], [0.4660, 0.6740, 0.1880]};
    
    % FIGURE: Annotated Signal Segments
    figure; hold on;
    grid on
    subplot(2, 1, 1); grid on; title("Annotated SCG Segments - Raw Features")
    xlabel("Timestep"); ylabel("Amplitude"); zlabel("Signal")
    subplot(2, 1, 2); grid on; title("Annotated SCG Segments - Resampled Features")
    xlabel("Timestep"); ylabel("Amplitude"); zlabel("Signal")
    for signal = 1:size(scg, 2)
        offset = signal*ones(length(scg(:, signal)), 1);
        subplot(2, 1, 1); hold on; 
        p = plot3(1:length(scg(:, signal)), scg(:, signal), offset, 'LineWidth', 1);
        p.Color(4) = 0.25; hold off
        subplot(2, 1, 2); hold on; 
        p = plot3(1:length(scg(:, signal)), scg(:, signal), offset, 'LineWidth', 1);
        p.Color(4) = 0.25; hold off
        % Plot features
        for feature = 1:numFeatures
            if ~isnan(SCG.raw{feature}(signal))
                subplot(2, 1, 1); hold on
                plot3(SCG.raw{feature}(signal), ...
                    scg(SCG.raw{feature}(signal), signal), signal, 'o', 'Color', colors{feature})
            end
            if ~isnan(finalIdx{feature}(signal))
                subplot(2, 1, 2); hold on
                plot3(finalIdx{feature}(signal), ...
                    scg(finalIdx{feature}(signal), signal), signal, 'o', 'Color', colors{feature})
            end
        end
    end
    % FIGURE: Feature Distribution
    figure; hold on; grid on; title("Feature Distribution")
    xlabel("Timestep"); ylabel("Probability")
    for feature = 1:numFeatures; area(x, norms{feature}, 'FaceAlpha', 0.3); end
    
end

end

