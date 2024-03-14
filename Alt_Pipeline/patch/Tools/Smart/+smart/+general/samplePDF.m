function samples = samplePDF(userPDF, numSamples, Plot)

% -------------------------------------------------------------------------
% This function generates pseudorandom samples based on a user-specified
% probability density function (PDF). Input arguments are as follows:
% - userPDF     [Nx1]   A normalized vector containing index probabilities
% - numSamples          Number of samples to draw from distribution
% - Plot        Bool    Plot results?
% -------------------------------------------------------------------------

% To protect against NaNs perform the following in a while loop
samples = nan;

while ~isempty(find(isnan(samples), 1)) || ~isempty(find(samples == 0, 1))

    x = 1:length(userPDF);  % Define x-axis
    
    % Plot PDF
    if Plot
        figure; subplot(3,1,1); hold on; grid on;
        title('Probability Density Function')
        plot(x, userPDF); xlabel('Value'); ylabel('Probability');
    end

    % Integrate PDF to obtain CDF
    CDF = cumsum(userPDF);

    % Plot CDF
    if Plot
        subplot(3,1,2); hold on; grid on;
        title('Cumulative Distribution Function')
        plot(x, CDF); xlabel('Value'); ylabel('Probability');
    end

    % Generating random numbers between 0 and 1 uniformly
    uniformSamples = rand(numSamples, 1);
    
    % If CDF is less than 1, ensure a sample isn't drawn from that range
    % (Accounting for lack of accuracy)
    if max(CDF) < 1; uniformSamples(uniformSamples > max(CDF)) = max(CDF); end

    % Find indices at which CDF is non-zero and less than 1
    indices = [false; not(diff(CDF) == 0)];  % Get proper indices
    CDF = CDF(indices); x = x(indices);     % Update vectors

    % Obtain samples by mapping uniformSamples to values using CDF
    % Find value in CDF that the uniformSample intersects (ceil)
    samples = zeros(numSamples, 1);     % Placeholder for PDF samples
    for i = 1:numSamples
        [~, temp] = min(abs(CDF - uniformSamples(i))); % Index of CDF value
        % Round the value up to the nearest CDF value
        if CDF(temp) < uniformSamples(i); temp = temp + 1; end
        samples(i) = x(temp);   % Get index of CDF value closest to sample
    end

    % Plot histogram
    if Plot
        subplot(3,1,3); hold on; grid on;
        title('Sample Histogram')
        histogram(samples)
        xlabel('Value'); ylabel('Frequency');

        % Set X and Y limits
        subplot(3,1,1); xlim([x(1), x(end)]); ylim([0, 1]);
        subplot(3,1,2); xlim([x(1), x(end)]); ylim([0, 1]);
        subplot(3,1,3); xlim([x(1), x(end)]);
    end
    
end

end

