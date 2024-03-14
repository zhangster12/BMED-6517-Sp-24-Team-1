function parent = getParent(~, level)

% -------------------------------------------------------------------------
% This function returns the parent for the specified level.
% -------------------------------------------------------------------------

% If "level" is not a list, return a single parent

if length(level) == 1

    % Return the parent
    parent = assignParent(level);
    
else
    
    % If "level" is a list, return a list of parents
    
    % Initialize placeholder for return value
    parent = level;
    
    % For each level...
    for l = 1:length(level)
        
        % Assign the parent
        parent(l) = assignParent(level(l));
        
    end

end

% Sub-function for assigning a parent
    function parent = assignParent(level)
        
        if level == Level.allRelative || level == Level.relBaseline1 || ...
                level == Level.relBaseline2 || level == Level.relative5 || ...
                level == Level.relative10 || level == Level.relative20
            parent = Level.allRelative;
        elseif level ~= Level.all
            parent = Level.allAbsolute;
        else
            parent = Level.all;
        end
        
    end

end