function offspring = shuffleAlleles(template)

% -------------------------------------------------------------------------
% This function generates a new organism from a template organism by
% shuffling the alleles at each locus. Input arguments include:
% - template:   {NxM}   Cell array containing N loci on M chromosomes
% -------------------------------------------------------------------------

% Set placeholder for return value
offspring = cell(size(template));

% Get number of loci and chromosomes
loci = size(template, 1); chromosomes = size(template, 2);

% for each locus...
for locus = 1:loci
    
    % Shuffle chromosomes at random
    shuffledAlleles = randperm(chromosomes);
    
    % Write alleles to new organism
    for i = 1:chromosomes
        offspring{locus, i} = template{locus, shuffledAlleles(i)};
    end
    
end

end

