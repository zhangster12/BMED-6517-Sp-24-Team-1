function [infectedPopulation, ViralSignature, figNum] = ...
    infect(population, penetration, vLoad, gen, pop, Plot)

% -------------------------------------------------------------------------
% This function generates a viral vector and introduces it in a population.
% Viral vectors may be introduced to increase genetic diversity to combat
% decreasing genetic diversity when local optima in the fitness functions
% are converged upon. Input parameters are as follows:
% - population  {NxM}   Cell array containing population data
% - penetration         Maximum penetration of viral vector
% - vLoad               Viral load (number of loci to affect)
%                       - 'min': Viral load is limited to one locus
% - gen                 Generation counter (from geneticSort.m)
% - pop                 Population number (from geneticSort.m)
% - Plot        Bool    Plot results?
% -------------------------------------------------------------------------

% Set placeholder for return value
infectedPopulation = population;

% Determine the loci to alter with the viral vector
numLoci = size(population{1}, 1);       % Get number of loci
numAlleles = size(population{1}, 2);    % Get number of alleles
maxLoci = floor(vLoad*numLoci);         % Maximum number of loci to infect
if maxLoci == 0; maxLoci = 1; end       % Possibility to infect at least one locus
% Total loci to infect
if ~strcmp(vLoad, 'min'); viralLoad = randi(maxLoci); else; viralLoad = 1; end  
% Determine loci to modify
modLoci = datasample(1:numLoci, viralLoad, 'Replace', false);

% Generate modified loci
modAlleles = cell(viralLoad, numAlleles);	% Placeholder for modified alleles
alleleOrder = cell(viralLoad, 1);           % Placeholder for allele order
for i = 1:viralLoad
    
    % Determine modified allele order
    alleleOrder{i} = randperm(numAlleles);
    
    % Generate modified alleles at each locus
    for j = 1:numAlleles
        modAlleles{i, j} = population{1}{modLoci(i), alleleOrder{i}(j)};
    end
    
end

% Determine number of organisms to infect
maxInfections = floor(penetration*length(population));
numInfections = randi(maxInfections);

% Determine specific organisms to infect
victims = datasample(1:length(population), numInfections, 'Replace', false);

% Infect organisms
for i = 1:numInfections
    for j = 1:viralLoad
        for k = 1:numAlleles
            infectedPopulation{victims(i)}{modLoci(j), k} = modAlleles{j, k};
        end
    end
end

% Record infections
infected = zeros(length(population), 1); infected(victims) = 1;

% Plot results
if Plot
    figure; subplot(1, 2, 1); hold on; figNum = get(gcf, 'Number');
    title("Infected Organisms for Population " + string(pop) + " at generation " + string(gen))
    imagesc(reshape(infected, [sqrt(length(population)), sqrt(length(population))]))
    xlim([0, sqrt(length(population)) + 1]); ylim([0, sqrt(length(population)) + 1]); 
    colorbar; subplot(1, 2, 2);
    for i = 1:viralLoad
        area(alleleOrder{i}, 'FaceAlpha', 0.3); xlabel("Chromosome"); ylabel("Allele");
        hold on
    end
    grid on; title("Viral Vector"); hold off
end

% Handle case where no figure is generated
if ~Plot; figNum = []; end

% Return viral signature as a struct
ViralSignature.loci = modLoci;          % List of modified loci
ViralSignature.alleles = modAlleles;    % Alleleic combination at each locus
ViralSignature.initial = victims;       % Initial set of victims

end

