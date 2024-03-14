function newIdx = map(oldIdx, path)

% -------------------------------------------------------------------------
% This function maps indices from a signal into its warped signal, given the
% old index and the warping path. If an old index maps to multiple values on
% the path, the middle value is taken.
% -------------------------------------------------------------------------

% Set placeholder for return value
newIdx = ones(size(oldIdx));

% For each old index
for idx = 1:length(oldIdx)
    
    % Find its corresponding indices in the warp path
    candidates = find(path == oldIdx(idx));
    
    if isempty(candidates); break; end
    
    % If there is more than one candidate, take the mean value
    if length(candidates) > 1
        newIdx(idx) = round(mean(candidates));
    else
        newIdx(idx) = candidates;
    end
    
end

end