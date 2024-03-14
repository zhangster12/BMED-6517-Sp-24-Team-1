function beat = sampleToBeat(~, samples, beats)

% Converting samples to beats indices. Samples is a vector of the sample
% numbers to convert and beats is a vector of samples containing R-peaks.

% Set placeholder for return value
beat = zeros(size(samples));

for i = 1:length(samples)
    [~, beat(i)] = min(abs(beats - samples(i)));
end

end