function [rawScores, normScores] = similarityScore(generation, fitness, fit_rank, L, Plot)

% -------------------------------------------------------------------------
% This function produces a symmetrical matrix containing the average
% Euclidian distance (and normalized distance) between organisms in each
% generation. Input parameters are as follows:
% - generation  {Nx1}   Cell array containing each organism
% - fitness     [Nx1]   Fitness score vector for N organisms
% - fit_rank    [Nx1]   Fitness ranks for each organism
% - L                   Progenitor limitation coefficient
% - Plot        Bool    Plot results?
% -------------------------------------------------------------------------

% Get number of organisms
pop = length(generation);

% Set placeholder for vectorized alleles
vectors = cell(size(generation));

% Initialize placeholder for return value
rawScores = zeros(pop);

% Limit computation of similarity score to top L organisms
if ~strcmp(L, 'all')
    fitness(fit_rank(L+1:end)) = 0;     % Zero out unfit organisms
    fitness = fitness/norm(fitness, 1); % Normalize fitness
end

% For each organism, get its distance from all other organisms
for i = 1:pop
    
    % Extract bases for first organism as vectors
    % If the organism's fitness is 0, do not perform this operation
    if i == 1
        organism1 = smart.general.plotCell(generation{i}, false);
        vectors{i} = organism1; else; organism1 = vectors{i}; 
    end
    
    for j = (i+1):pop
        
        % For each chromosome, extract the bases as vectors and compute the
        % distance as the L2 norm (Euclidian distance).
        if i == 1
            organism2 = smart.general.plotCell(generation{j}, false);
            vectors{j} = organism2; else; organism2 = vectors{j};
        end
        
        if fitness(i) == 0 || fitness(j) == 0
            distance = 0; else; distance = norm(organism1 - organism2);
        end
        
        % Write distance to return value (matrix is symmetrical)
        rawScores(i,j) = distance; rawScores(j,i) = distance;
        
    end
    
end

% Normalize the matrix
normScores = rawScores/norm(rawScores);

% Plot results
if Plot; figure; title("Similarity Scores"); imagesc(rawScores); colorbar; end

end

