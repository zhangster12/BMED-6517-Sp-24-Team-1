function score = evaluate(signal, idx, means)

% -------------------------------------------------------------------------
% This function computes a performance score for how well the given signal
% was clustered based on the provided indices. Input arguments inclue:
% - signal  [2xN]   Vector of points to cluster
% - idx     [Nx1]   Cluster of each point
% - means   [2xM]   Mean points for each cluster
%
% The current metric is a silhouette score, defined as:
%
% s = (b-a)/max(a,b)
%
% Where a is the mean distance between a sample and all other points in the
% same class and b is the mean distance between a sample and all other
% points in the next nearest cluster.
% -------------------------------------------------------------------------

% Remove outliers and decrement indices
signal(:, idx == 1) = []; idx(idx == 1) = []; idx = idx - 1;

% Get number of clusters and points
C = max(idx); N = length(signal);

% Set placeholder for silhouette coefficient for each sample
s = zeros(N, 1); a = zeros(N, 1); b = zeros(N, 1);

% For each cluster, find all points in cluster
points = cell(C, 1);    % Set placeholder
for c = 1:C
    points{c} = find(idx == c);
end

% For each cluster, find all points in next nearest cluster
% First, find next-nearest cluster for each cluster
nextNearest = zeros(C, 1);  % Set placeholder for nearest neighbor
for c = 1:C
    temp = zeros(1, C); % Set placeholder for distances
    % Find distance between current mean and all other means
    for c2 = 1:C
        temp(:, c2) = norm(means(:, c) - means(:, c2));
    end
    % Find minimum distance, accounting for nonzero index
    [~, index] = sort(temp); nextNearest(c) = index(2);
end

% Find all points in next nearest cluster
notPoints = cell(C, 1);     % Set placeholder
for c = 1:C
    notPoints{c} = find(idx == nextNearest(c));
end

% For each point...
for point = 1:N
    
    % Determine class of current point
    class = idx(point);
    
    % Find value of (a)
    % 1. Find the total distance between this point and others within class
    dist = 0;   % Initialize counter
    for point2 = 1:length(points{class})
        temp = abs(signal(:, point) - signal(:, points{class}(point2)));
        dist = temp(1);	% Get vertical distance only (for timeseries)
        % dist = dist + norm(signal(:, point) - signal(:, points{class}(point2)));
    end
    
    % 2. Normalize by number of points and write to (a)
    a(point) = dist/length(points{class});
    
    % Find value of (b)
    % 1. Find the total distance between this point and points in other classes
    dist = 0;   % Initialize counter
    for point2 = 1:length(notPoints{class})
        temp = abs(signal(:, point) - signal(:, notPoints{class}(point2)));
        dist = temp(1);     % Get vertical distance only (for timeseries)
        % dist = dist + norm(signal(:, point) - signal(:, notPoints{class}(point2)));
    end
    
    % 2. Normalize by number of points and write to (b)
    b(point) = dist/length(notPoints{class});
    
    % Compute silhouette score
    s(point) = (b(point) - a(point))/max(a(point),b(point));
    
end

% Return the overall score as the mean silhouette score
score = nanmean(s);

end
