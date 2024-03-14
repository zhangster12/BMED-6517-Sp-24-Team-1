function SCG = ...
    geneticSCG(signals, rawFeatures, numFeatures, M, tol, Plot, Detrend, varargin)

% -------------------------------------------------------------------------
% SUMMARY
% This function extracts the time-domain features of an SCG signal with
% potential transient corruptions and discontinuities. To remove
% discontinuities, the following steps are performed. Once a feature is
% extracted from the signal using the traditional peak-counting method,
% k-means clustering is used to separate signal segments for which there
% are discontinuities. The signals are then split at these indices and
% recombined using a genetic algorithm such that the resulting feature
% vectors have the lowest detrended variability and are thereby properly
% sorted to remove discontinuities and outliers.
%
% ARGUMENTS (REQUIRED)
% - signals     [NxM]   M signal vectors of N timepoints
% - rawFeatures         Number of time-domain features to extract
% - numFeatures         Number of feature vectors to sort ( < rawFeatures)
% - M                   Smoothing factor for exponential moving average
% - tol                 Tolerance for outliers (in standard deviations)
% - Plot        Bool    Plot results?
% - Detrend     Bool    Detrend features?
%
% ARGUMENTS (OPTIONAL)
% - K               Maximum cluster size for K means
% - S               Scaling factor
%
% OUTPUTS
% SCG           Struct
%   .features       Cell array of feature vectors
%
% USAGE
% scg = cardio.scg.geneticSCG(signals, rawFt, numFt, M, tol, Plot, Detrend, K)
% scg = cardio.scg.geneticSCG(signals, rawFt, numFt, M, tol, Plot, Detrend, K, S)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Part 0: Setup
% -------------------------------------------------------------------------

% Parse optional arguments if necessary
if nargin > 7
    for i = 1:length(varargin)
        if ~exist('K', 'var'); K = varargin{i}; 
        elseif ~exist('S', 'var'); S = varargin{i}; 
        end
    end
end

% Set optional arguments if not provided
if ~exist('K', 'var'); K = 2; end; if ~exist('S', 'var'); S = 1; end

% Get toal number of signals
numSignals = size(signals, 2);

% Initialize global arrays
indices = cell(rawFeatures, 1);     % Raw indices derived from peak-counting
breakpoints = cell(rawFeatures, 1);	% Indices of discontinuities in each signal

% -------------------------------------------------------------------------
% Part 1: Exponential Moving Average
% -------------------------------------------------------------------------

% Arguments: signal, M, plotResults
scg = cardio.general.ema(signals, M, false);

% -------------------------------------------------------------------------
% Part 2: Extract Features
% -------------------------------------------------------------------------

% Set starting index to search for features
startIdx = ones(size(scg, 2), 1);

% For each feature...
for feature = 1:rawFeatures
    
    % ---------------------------------------------------------------------
    % Part 2(a): Extract Features with Peak-Counting
    % ---------------------------------------------------------------------
    
    % Even features correspond to maxima -> odd features correspond to minima
    if mod(feature,2) == 0; findMax = true; else; findMax = false; end
    
    % Set placehodler for feature indices
    indices{feature} = zeros(numSignals, 1);
    
    % Set placeholders for peaks and valleys in each segment
    peaks = cell(numSignals, 1); valleys = cell(numSignals, 1);
    
    % Find desired feature type along each signal segment
    for signal = 1:numSignals
        
        % Extract segment
        sig = scg(startIdx(signal):end, signal);
       
        % Find all features of the given type
        [peaks{signal}, valleys{signal}, ~] = cardio.general.getPeaks(sig);
        
        % Add back start indices
        peaks{signal} = peaks{signal} + startIdx(signal);
        valleys{signal} = valleys{signal} + startIdx(signal);
        
        % Repalce empty vectors with NaN
        if isempty(peaks{signal}); peaks{signal} = NaN; end
        if isempty(valleys{signal}); valleys{signal} = NaN; end
        
        % Return the first feature within the search range
        if findMax; indices{feature}(signal) = min(peaks{signal}); else; ...
                indices{feature}(signal) = min(valleys{signal}); end
        
    end
    
    % ---------------------------------------------------------------------
    % Part 2(b): Cluster Features with K-Means
    % ---------------------------------------------------------------------
    
    % Detrend and scale before clustering
    if Detrend
        dt_indices = cardio.general.detrendFeatures({indices{feature}});
        dt_indices = dt_indices{1}; % Extract from cell array
    else
        dt_indices = indices{feature};
    end
    
    % Scale indices by optional scaling factor (separation for K-means)
    scaled_indices = dt_indices.^S;
    
    % Perform K-Means Clustering
    clusters = kmeans(scaled_indices, K); numClusters = length(unique(clusters));
    
    % ---------------------------------------------------------------------
    % Part 2(c): Obtain Breakpoints in Signal
    % ---------------------------------------------------------------------
    
    % Differentiate cluster signal to determine breakpoint locations
    breaks = diff(clusters); breakpoints{feature} = find(breaks ~= 0);
    
    % ---------------------------------------------------------------------
    % Part 2(d): Update Start Indices
    % ---------------------------------------------------------------------
    
    % Update start index for next feature (fill NaNs with nearest neighbor)
    startIdx = fillmissing(indices{feature},'nearest');
    
end


% ---------------------------------------------------------------------
% Part 3: Consolidate Breakpoints
% ---------------------------------------------------------------------

% Concatenate cell array and keep only unique values
rawIndices = smart.general.plotCell(breakpoints, false);
brkIndices = unique(rawIndices);

% ---------------------------------------------------------------------
% Part 4: Sort Features with Genetic Algorithm
% ---------------------------------------------------------------------

% Create population(s)
% Arguments: mixedVectors, B, R, Mu, C, L, P, N, Mi, F, V
Population1 = ...
    smart.genetic.definePopulation([indices{:}], ...
    brkIndices, 12, 15, 15, 10, 2, 2, 0.01, 'scgmean', 'none');
Population2 = ...
    smart.genetic.definePopulation([indices{:}], ...
    brkIndices, 12, 20, 20, 10, 2, 2, 0.01, 'scgmax', 'none');

% Arguments: patience, maxGen, Plot, Population1, Population2, ...
stdIdx = smart.genetic.geneticSort(5, 1000, false, Population1, Population2);

SCG.features = stdIdx;

% ---------------------------------------------------------------------
% Part 5: Remove Outliers
% ---------------------------------------------------------------------

% Set placeholder for return value
finalIdx = stdIdx;

% For each feature vector
for feature = 1:numFeatures
    
    % Get the feature vector
    ftVector = stdIdx(:, feature);
    
    % Detrend feature vector
    dtVector = cardio.general.detrendFeatures({ftVector}); dtVector = dtVector{1};
    
    % Compute bounds
    minVal = nanmean(dtVector) - tol*nanstd(dtVector);
    maxVal = nanmean(dtVector) + tol*nanstd(dtVector);
    
    % Replace an index with NaN if it lies outside the tolerance
    ftVector(ftVector < minVal) = NaN; ftVector(ftVector > maxVal) = NaN;
    
    % Write result to finalIdx
    finalIdx(:, feature) = ftVector;

end

% Return final features
SCG.features = finalIdx;

% -------------------------------------------------------------------------
% Part 6: Calculate Statistics
% -------------------------------------------------------------------------

% Initialize placeholder for statistics
statistics = zeros(numFeatures, 2); % [Mean] [Standard Deviation]
norms = cell((numFeatures), 1);     % Cell array for norm PDFs
x = 0:size(SCG.features, 1);        % Define x-coordinate
% Return statistics for each feature
for feature = 1:numFeatures
    statistics(feature, 1) = nanmean(SCG.features(:, feature));    % Mean
    statistics(feature, 2) = nanstd(SCG.features(:, feature));     % Stddev
    % Fit norm PDF
    norms{feature} = normpdf(x, statistics(feature,1), statistics(feature,2));
end

% ---------------------------------------------------------------------
% Part 7: Visualization
% ---------------------------------------------------------------------

if Plot

    % FIGURE: Annotated Signal Segments
    figure; hold on; grid on; title("Annotated SCG Segments")
    xlabel("Timestep"); ylabel("Amplitude")
    for signal = 1:numSignals
        plot(scg(:,signal)) % Plot signal
        % Plot features
        for feature = 1:numFeatures
            if ~isnan(indices{feature}(signal))
                plot(indices{feature}(signal), ...
                    scg(indices{feature}(signal), signal), 'or')
            end
            if ~isnan(finalIdx(signal, feature))
                plot(finalIdx(signal, feature), ...
                    scg(finalIdx(signal, feature), signal), 'ok')
            end
        end
    end

    % FIGURE: Feature Distribution
    figure; hold on; grid on; title("Feature Distribution")
    xlabel("Timestep"); ylabel("Probability")
    for feature = 1:numFeatures; area(x, norms{feature}, 'FaceAlpha', 0.3); end
    
end

end

