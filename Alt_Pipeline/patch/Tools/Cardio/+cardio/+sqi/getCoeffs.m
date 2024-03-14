function coeffs = getCoeffs(S, varargin)

% This function returns the fourier coefficients of S in standardized form.
% - 'combine'   Combine coefficients in same window and frequency

% Parse optional arguments
if ~isempty(varargin) && strcmp(varargin{1}, 'combine'); combine = true; else; ...
        combine = false; end

% Organize the error vector such that a mask can be applied to it
coeffs = [];   % Initialize placeholder
for win = 1:length(S.coeffs)	% For each window...
    temp = S.coeffs{win}'; temp = temp(:); coeffs = [coeffs; temp];
end

% If 'combine', then average every two components
if combine
    % Combine components
    temp = zeros(length(coeffs)/2, 1); c = 1;
    for i = 1:2:length(coeffs)
        temp(c) = (coeffs(i) + coeffs(i + 1))/2; c = c + 1;
    end
    % Update coefficients
    coeffs = temp;
end
    
end

