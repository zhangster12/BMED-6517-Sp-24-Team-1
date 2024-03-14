function circlePoints = project(data, radius, varargin)

% -------------------------------------------------------------------------
% This function projects a 2-D dataset onto a circular function defined by
% its radius, centered at the origin.
%
% Arguments (required)
% - data        [MxN]       N data vectors of length M
% - radius                  Radius of circle
%
% Arguments (optional)
% - start                   Angle (in radians) to start circle
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'start'); start = varargin{arg + 1};
        end
    end
end

% Set defaults for optional input arguments
if ~exist('start', 'var'); start = 0; end

% Get the number of points to evaluate
numPoints = size(data, 1);

% Set placeholder for return values
delta = zeros(numPoints, 1);
circlePoints = zeros(numPoints, 1);

% Generate (x, y) coordinates
angles = linspace(start, 2*pi + start, 1000);
x = radius*cos(angles);
y = radius*sin(angles);
circumference = 2*pi*radius;
displacement = linspace(0, circumference, 1000);

% For each point...
for i = 1:numPoints
    
    % Set placeholder for distances/angles
    distances = zeros(length(x), 1);
    
    % Get coordinate of current point
    coord = [data(i, 1); data(i, 2)];
    
    % For each circle point...
    for j = 1:1000
        
        % Compute the distance
        distances(j) = norm(coord - [x(j); y(j)]);
        
    end
    
    % Get the index of the point with the lowest distance
    [~, idx] = min(distances);
    
    % Get the displacement
    delta(i) = displacement(idx);
    
    % Set initial point
    if i == 1; circlePoints(i) = 0; continue; end
    
    % Get displacement
    if i > 1
        
        % Compute displacement
        circlePoints(i) = delta(i) - delta(i - 1);
        
        % If there is an error, fix it
        if delta(i) < 0.1*circumference && delta(i - 1) > 0.9*circumference
            circlePoints(i) = circumference - (delta(i - 1) - delta(i));
        elseif delta(i - 1) < 0.1*circumference && delta(i) > 0.9*circumference
            circlePoints(i) = (delta(i) - delta(i - 1)) - circumference;
        end
            
        
    end
    
end

% Return the cumulative sum of the offsets
circlePoints = cumsum(circlePoints);

end

