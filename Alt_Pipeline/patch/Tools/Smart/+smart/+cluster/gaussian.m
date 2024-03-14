function output = gaussian(signal, Plot, varargin)

% -------------------------------------------------------------------------
% SUMMARY
% This function clusters an input signal using expectation maximization and
% Gaussian Mixture Models (GMMs) for classification. Currently, the
% function builds and tests up to a threshold number of clusters,
% evaluating each model with a silhouette score to determine the optimal
% number of clusters.
%
% ARGS (REQ'D)
% - signal      [Nx1]   Signal vector
% - Plot        Bool    Plot results?
%
% ARGS (OPT)
% - varargin
%   - maxClusters           Maximum number of clusters (> 2, < 3)
%   - numClusters           Specify a fixed number of clusters (> 2)
%   - iterations            Algorithm iterations per cluster
%   - startModel    GMM     Starting GMM from previous iteration
%
% OUTPUTS (OPT)
% - indices
% - startModel
%   - means             Mean values of Gaussian distributions (optimal)
%   - sigma             Covariance matrices of Gaussian distributions (optimal)
%   - prior             Prior probabilities of Gaussian distributions (optimal)
%   - opClusters        Optimal number of clusters
%
% NOTES
% Note that the current function is limited to 2-3 clusters due to custom
% initiation of priors, however this script may be altered to accept any
% number of clusters by automating prior and mean selection.
%
% USAGE
% idx = smart.cluster.gaussian(signal, Plot);
% idx = smart.cluster.gaussian(signal, Plot, 'maxClusters', 3);
% idx = smart.cluster.gaussian(signal, Plot, 'iterations', 100);
% idx = smart.cluster.gaussian(signal, Plot, 'maxClusters', 3, 'iterations', 100);
% -------------------------------------------------------------------------

% Parse optional arguments
if nargin > 2
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'maxClusters')
            maxClusters = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'iterations')
            iterations = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'startModel')
            startModel = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'numClusters')
            numClusters = varargin{arg + 1};
        end
    end
end

% Set defaults
if ~exist('maxClusters', 'var'); maxClusters = 3; end   % Default 3 clusters
if ~exist('iterations', 'var'); iterations = 100; end   % Default 100 iterations
if size(signal, 2) == 2; twoDim = true; else; twoDim = false; end

% Adjust for errors
if maxClusters < 2 || maxClusters > 3; maxClusters = 3; end

% Set placeholders
ind = cell(maxClusters, 1);     % Indices for each model
scores = zeros(maxClusters, 1); % Scores for each model

% Set placeholders
c_means = cell(maxClusters - 1, 1); % Placeholder for Gaussian means
c_sigma = cell(maxClusters - 1, 1); % Placeholder for Gaussian covariances
c_prior = cell(maxClusters - 1, 1); % Placeholder for Gaussian priors

counter = 1;  % Initialize counter

% Set range of clusters over which to test
if exist('numClusters', 'var'); range = numClusters;
else; range = 2:maxClusters; numClusters = maxClusters;
end

for clusters = range
    % -------------------------------------------------------------------------
    % Process Signal
    % -------------------------------------------------------------------------

    if ~twoDim
        
        % Normalize and detrend signal
        nsig = normalize(signal);

        % Remove large signal outliers
        nsig(nsig < -3) = NaN; nsig(nsig > 3) = NaN; sig = fillmissing(nsig, 'nearest');

        % Format signal as two-column vector (time, signal)
        x = [sig'; linspace(-2, 2, length(sig))]; N = length(sig);
        
    else
        
        % Extract dimensions
        dim1 = signal(:, 1); dim2 = signal(:, 2);
        % ndim1 = normalize(signal(:, 1)); ndim2 = normalize(signal(:, 2));
        
        % Format signal as two-column vector
        x = [dim1 dim2]'; N = length(dim1);
        
    end

    % -------------------------------------------------------------------------
    % Define Multivariate Gaussian Distribution
    % -------------------------------------------------------------------------

    % Gaussian Distribution (phi)
    % Arguments:
    % - x:      Index for evaluation
    % - mu:     Mean of distribution
    % - var:    Variance of distribution
    phi = @(x,mu,var) (1/(2*pi*sqrt(det(var))))*exp(-0.5*(x-mu).'*(var\(x-mu)));

    % -------------------------------------------------------------------------
    % Initialize Parameters
    % -------------------------------------------------------------------------

    % Initialize means, if necessary
    if ~exist('startModel', 'var')

        % Set placeholders for distribution means
        means = zeros(2, clusters);  % [mean(x), mean(y)]

        % Initialize distribution centers for each cluster
        if clusters == 1; means(:, 1) = [0, 0];
        elseif clusters == 2; means(:, 1) = [2 0]; means(:, 2) = [-2 0];
        else; means(:, 1) = [2 0]; means(:, 2) = [-2 0]; means(:, 3) = [0 0];
        end

    else
        means = startModel.means;    % Set initial values if provided
    end

    % Distribute priors for initial guesses, if necessary
    if ~exist('startModel', 'var')
        if clusters == 3; prior = [1/4; 1/4; 1/2];
        elseif clusters == 2; prior = [1/2; 1/2];
        else; prior = 1;
        end
    else
        prior = startModel.prior;   % Set initial values if provided
    end

    % Initialize placeholders for variance of each distribution (sigma)
    if ~exist('startModel', 'var'); sigma = cell(clusters, 1); sigma(:, 1) = {eye(2)};
    else
        sigma = startModel.sigma;   % Set initial values if provided
    end

    % Visualize initial distributions
    if Plot; visualize("Initial Distributions"); end

    % -------------------------------------------------------------------------
    % Update Distributions
    % -------------------------------------------------------------------------

    % Set placeholder for probability that each point falls in each distribution
    posterior = zeros(N, clusters);

    % Perform the update K times
    for loop = 1:iterations

        % Show distribution at each loop
        % visualize("Distribution at Loop " + string(loop))

        % For each point...
        for p = 1:N

            % Set placeholder for normalizer (sum of all distributions)
            den = 0;

            % Find the total distribution area scaled by the priors
            for c = 1:clusters
                den = den + prior(c)*phi(x(:, p), means(:, c), sigma{c});
            end

            % For each cluster, find the posterior distribution
            for c = 1:clusters
                if clusters > 1
                    post = prior(c)*phi(x(:, p), means(:, c), sigma{c});
                    posterior(p, c) = post/den;
                else; posterior(p, c) = 1;
                end
            end

        end

        % Update the priors and means
        post_sum = sum(posterior);  % Get sum of posteriors for each cluster
        if clusters > 1
            prior = (1/N)*post_sum; means = (x*posterior)./post_sum;
        else; prior = 1;
        end

        % Update covariance
        temp = zeros(2);    % Initialize placeholder
        for c = 1:clusters  % For each cluster...
            for p = 1:N     % For each point in cluster...
                % Update placeholder by covariance for each point
                temp = temp + posterior(p, c)*(x(:, p) - means(:, c))*...
                    (x(:, p) - means(:, c)).';
            end
            % Handle case where cov is 0
            if temp(1, 1) == 0; temp(1, 1) = eps; end
            if temp(2, 2) == 0; temp(2, 2) = eps; end
            % Normalize covariance and update placeholder
            temp = temp./post_sum(c); sigma{c} = temp;
        end

    end
    
    % Classify points
    ind{clusters} = classify();
    
    % Get score for the current model
    if numClusters > 1
        scores(clusters) = smart.cluster.evaluate(x, ind{clusters}, means);
    else
        scores(clusters) = 1;
    end
    
    % Visualize results
    if isnan(scores(clusters)); scores(clusters) = 0; end
    if Plot; visualize("Model: " + string(clusters - 1) + ", Score: " + ...
            string(scores(clusters))); end
    
    % Save means, covariances, and priors
    c_means{counter} = means; c_sigma{counter} = sigma; c_prior{counter} = prior;
    
    % Increment counter
    counter = counter + 1;

end

% Return indices of optimal model
if ~isempty(scores(scores ~= 0))
    [~, op] = max(scores(scores ~= 0)); indices = ind{op + 1};
else
    op = 1; indices = ind{clusters};
end

% Print result
if Plot; disp("Optimal Model: " + string(op) + ", Score: " + ...
        string(scores(op + 1))); end

% Return configuration of optimal model
output.means = c_means{op};
output.sigma = c_sigma{op};
output.prior = c_prior{op};
output.numClusters = length(output.prior);
output.model = output;
output.indices = indices;

% -------------------------------------------------------------------------
% FUNCTION: Classify Points
% -------------------------------------------------------------------------

    function [idx, prob] = classify()
        
        % Determine the probability that each point falls within each cluster
        prob = zeros(N, clusters + 1);  % Set placeholder for probabilities
        
        % For each point...
        for pt = 1:N
            
            % If the point is an outlier, automatically asign it to class 0
            if ~twoDim && isnan(nsig(pt)); prob(pt, 1) = 1; continue; end
            if twoDim && (isnan(dim1(pt)) || isnan(dim2(pt))); prob(pt, 1) = 1; continue; end
            
            % For each cluster
            for cl = 1:clusters
                % Get probability
                prob(pt, cl + 1) = prior(cl)*phi(x(:, pt), means(:, cl), sigma{cl});
            end
        end
        
        % For each point, find the cluster for which the probability is max
        [~, idx] = max(prob, [], 2);
        
    end

% -------------------------------------------------------------------------
% FUNCTION: Visualize Results
% -------------------------------------------------------------------------

    function visualize(figTitle)
        
        if ~twoDim
            % Define test region
            x1 = linspace(-2, 2, length(sig));
            x2 = linspace(min(sig), max(sig), length(sig));
        else
            x1 = linspace(min(dim1), max(dim1), length(dim1));
            x2 = linspace(min(dim2), max(dim2), length(dim2));
        end

        % Placeholders for contours
        X = cell(clusters, 1); X(:) = {zeros(length(x1), length(x2))};
        
        c1 = 1;         % Loop counter
        for i = x1      % y-coordinate
            c2 = 1;     % Loop counter
            for j = x2  % x-coordinate
                test = [j; i];    % Set test point
                % Populate contours
                for cl = 1:clusters
                    X{cl}(c1, c2) = prior(cl)*phi(test, means(:, cl), sigma{cl});
                end
                c2 = c2 + 1;    % Increment counter
            end
            c1 = c1 + 1;        % Increment counter
        end
        
        % Classify the signal points
        idx = classify();
        
        % Plot the signals and contours
        figure; hold on; grid on; title(figTitle)
        xlabel('Signal Segment'); ylabel('Amplitude');
        colors = [".k"; ".r"; ".b"; ".g"; ".y"];
        for cl = 1:clusters
            plot(x(2, idx == cl + 1), x(1, idx == cl + 1), colors(cl))
            contour(x1, x2, X{cl}')
        end
        
    end

end
