% -------------------------------------------------------------------------
% Genetic Sorting Algorithm
% Created By: Jonathan Zia
% Date Created: 12/28/2018
% Inan Research Lab @ Gatech
% -------------------------------------------------------------------------

function sortedVectors = ...
    geneticSort(patience, maxGen, Plot, varargin)

% -------------------------------------------------------------------------
% SUMMARY:
% This function implements a genetic algorithm to sort a set of time-series
% vectors, the elements of which may be scrambled due to changes in the
% morphology of the underlying signal from which the feature vectors were
% extracted. Input parameters are as follows:
% - mixedVectors    [NxM]   An array of columnwise feature vectors
% - B               HYPERPARAMETER: Bases per allele
%                   - Scalar or vector of indices
% - R               HYPERPARAMETER: Population size (R^2)
% - Mu              HYPERPARAMETER: Maximum mutations per organism
% - C               HYPERPARAMETER: Maximum crossover points per generation
% - L               HYPERPARAMETER: Limiting breeding to top L organisms
%                   - 'all': All organisms have opportunity to breed
% - P               HYPERPARAMETER: Propensity score (bias toward selecting
%                                   parents with relatively higher fitness)
% - N               HYPERPARAMETER: Niche penalty (bias toward diversity)
% - Mi              HYPERPARAMETER: Migration coefficient
% - Pp              HYPERPARAMETER: Number of parallel populations
% - F               HYPERPARAMETER: Method of calculating fitness
%                   - 'avgstd': Mean standard deviation
%                   - 'maxstd': Maximum standard deviation
% - V               HYPERPARAMETER: Number of generations between viral pandemics
% - patience        Stop patience (number of generations)
% - maxGen          Maximum number of generations
% - Plot            Plot migration / pandemics? 

% These input parameters are organized as part of the Population struct and
% are parsed accordingly. This struct may be constructed with the
% smart.structs.definePopulation() method.

% DESCRIPTION:
% This function first accepts a columnwise array of time-series vectors and
% divides them into cell arrays in which multiple bases are grouped
% together into alleles (cells) based on the hyperparameter B, which
% specifies the number of bases per allele. This step increases the
% efficiency of the algorithm by reducing the total number of possible
% combinations, while conversely decreasing the accuracy of the result.
% This array is then recombined 2^R times to produce 2^R distinct vector
% sets ("organisms"). The fitness of each organism is determined by the
% detrended standard deviation, based on the assumption that the signals
% will minimize the standard deviation when properly sorted. The organisms
% are then mutated based on the parameter M and recombined with other
% organisms based on their fitness and the crossover parameter C, and the
% process is repeated. The process is repeated until the peak fitness is no
% longer increasing (1/max_std), at which point the paramter M is decreased
% and the process is repeated. This continues until M has been decayed to
% its minimum point, and the results are displayed.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Parse and Initialize Variables
% -------------------------------------------------------------------------

Pp = nargin - 3;    % Number Population objects -> number of populations

% Initialize global variables
mixedVectors = cell(Pp, 1); B = cell(Pp, 1); R = cell(Pp, 1); 
Mu = cell(Pp, 1); C = cell(Pp, 1); L = cell(Pp, 1); P = cell(Pp, 1);
N = cell(Pp, 1); Mi = cell(Pp, 1); F = cell(Pp, 1); V = cell(Pp, 1);

for pop = 1:Pp
   
    % Set variables for each population
    mixedVectors{pop} = varargin{pop}.mixedVectors;
    B{pop} = varargin{pop}.B; R{pop} = varargin{pop}.R;
    Mu{pop} = varargin{pop}.Mu; C{pop} = varargin{pop}.C;
    L{pop} = varargin{pop}.L; P{pop} = varargin{pop}.P;
    N{pop} = varargin{pop}.N; Mi{pop} = varargin{pop}.Mi;
    F{pop} = varargin{pop}.F; V{pop} = varargin{pop}.V;
    
end

% Set placeholders for necessary universal variables
generation = cell(Pp, 1);   % Placeholder for each generation
fitness = cell(Pp, 1);      % Normalized fitness fore each individual
fit_rank = cell(Pp, 1);     % Fitness rankings for each individual
fitness_op = cell(Pp, 1);   % Optimal fitness for each population
fitOrganism = cell(Pp, 1);  % Optimal organism in each population
kingOfTheHill = cell(Pp, 1);% King of the hill in each population
optimal_fitness = cell(Pp, 1);  % Optimal fitness in each population
average_fitness = cell(Pp, 1);  % Average fitness in each population
loopFLAG = true;            % Loop flag for each population
loopCounter = 1;            % Loop counter for generations
unimproved = cell(Pp, 1);   % Initializing stop parameter
numChromosomes = size(mixedVectors{1}, 2); % Number of vectors (chromosomes)
similarity = cell(Pp, 1);   % Similarity scores for each population
diversity = cell(Pp, 1);    % Diversity of each population
samples = cell(Pp, 1);      % Sampling frequency per parent
virusCounter = cell(Pp, 1); % Flag indicating how many viruses have been released
Fig = cell(Pp, 1);          % Figure numbers fore each population
ViralSig = cell(Pp, 1);     % Viral signatures for each population

% Initialize placeholders
unimproved(:) = {0}; virusCounter(:) = {0};

% Perform initial processing for each community
parfor pop = 1:Pp

    % ---------------------------------------------------------------------
    % Part 1: Process Original Vectors
    % ---------------------------------------------------------------------

    % Organize the chromosomes into cell arrays divided into alleles
    % Arguments: inputArray, B
    origChromosomes = ...
        smart.genetic.generateAlleles(mixedVectors{pop}, B{pop});

    % ---------------------------------------------------------------------
    % Part 2: Generate Initial Generation
    % ---------------------------------------------------------------------

    % Arguments: old_gen, fitness, fit_rank, similarity, R, M, C, L, N
    [generation{pop}, ~] = smart.genetic.produceGeneration(origChromosomes, ...
        [], [], [], R{pop}, Mu{pop}, C{pop}, L{pop}, N{pop});

    % ---------------------------------------------------------------------
    % Part 3: Evaluate Initial Generation
    % ---------------------------------------------------------------------

    % For each organism, obtain fitness, normalized fitness score, and rank
    [fit_raw, fitness{pop}, fit_rank{pop}] = ...
        smart.genetic.getFitness(generation{pop}, P{pop}, F{pop});

    % Record optimal fitness
    fitness_op{pop} = fit_raw(fit_rank{pop}(1));

    % Set placeholders for training data
    optimal_fitness{pop} = []; average_fitness{pop} = []; diversity{pop} = [];
    optimal_fitness{pop} = [optimal_fitness{pop} fit_raw(fit_rank{pop}(1))];
    average_fitness{pop} = [average_fitness{pop} mean(fit_raw)];
    
    % Generate similarity matrix if N is provided
    if N{pop} > 0
        [similarity{pop}, ~] = smart.general.similarityScore(generation{pop}, ...
            fitness{pop}, fit_rank{pop}, L{pop}, false);
    else; similarity{pop} = ones(2^R{pop}, 2^R{pop});
    end

end

% ---------------------------------------------------------------------
% Plot Results
% ---------------------------------------------------------------------

% Generate plots for each population
for pop = 1:Pp
    
    % Generate a plot to display original vectors
    figure(pop); subplot(2,3,1); hold on; grid on;   % Format figure
    title("Original Vectors"); xlabel("Timestep"); ylabel("Amplitude")
    for i = 1:numChromosomes
        area(mixedVectors{pop}(:,i), 'FaceAlpha', 0.3)
    end; xlim([0, length(mixedVectors{pop})]); hold off

    % Display heatmap of fitness
    fitMatrix = reshape(fitness{pop}, [sqrt(2^R{pop}), sqrt(2^R{pop})]);
    figure(pop); subplot(2,3,2); hold on
    imagesc(fitMatrix); title("Initial Fitness per Organism"); colorbar; 
    xlim([0,sqrt(2^R{pop})+1]); ylim([0,sqrt(2^R{pop})+1]); hold off

    % Plot organism with highest fitness
    figure(pop); subplot(2,3,4); grid on;
    title("Fittest Organism - Generation 0")
    xlabel('Timestep'); ylabel('Amplitude');
    fitOrganism{pop} = smart.general.plotCell(generation{pop}{fit_rank{pop}(1)}, false);
    hold on
    for i = 1:numChromosomes
        area(fitOrganism{pop}(:,i), 'FaceAlpha', 0.3)
    end; xlim([0, length(mixedVectors{pop})]); hold off
    
    % Plot optimal fitness and average fitness of population
    figure(pop); subplot(2,3,3)
    area(0, optimal_fitness{pop}, 'FaceAlpha', 0.3); hold on; grid on
    area(0, average_fitness{pop}, 'FaceAlpha', 0.3)
    title("Population Fitness for Generation 0/" + ...
        string(maxGen) + " - Patience Counter: " + string(patience))
    xlabel("Generation"); ylabel("Fitness")
    legend("Optimal Fitness", "Overall Fitness"); hold off
    
end


while loopFLAG
    
    parfor pop = 1:Pp
    
        % -----------------------------------------------------------------
        % Part 4: Generate Subsequent Generation
        % -----------------------------------------------------------------

        % Arguments: old_gen, fitness, fit_rank, similarity, R, M, C, L, N
        [new_generation, samples{pop}] = smart.genetic.produceGeneration(generation{pop}, ...
            fitness{pop}, fit_rank{pop}, similarity{pop}, R{pop}, Mu{pop}, C{pop},...
            L{pop}, N{pop});

        % -----------------------------------------------------------------
        % Part 5: Evaluate Subsequent Generation
        % -----------------------------------------------------------------
        
        % For each organism, obtain fitness, normalized fitness score, and rank
        [fit_raw, fitness{pop}, fit_rank{pop}] = ...
            smart.genetic.getFitness(new_generation, P{pop}, F{pop});
        
        % Generate similarity matrix if N is provided
        if N{pop} > 0
            [similarity{pop}, ~] = ...
                smart.general.similarityScore(generation{pop}, fitness{pop}, ...
                fit_rank{pop}, L{pop}, false);
        else; similarity{pop} = ones(2^R{pop}, 2^R{pop});
        end

        % -----------------------------------------------------------------
        % Part 6: Assess Progress
        % -----------------------------------------------------------------

        % Record change in optimal fitness
        % Get change in fitness
        fitness_delta = fit_raw(fit_rank{pop}(1)) - fitness_op{pop};
        if fitness_delta > 0    % If there has been an improvement
            % Record new optimal value
            fitness_op{pop} = fit_raw(fit_rank{pop}(1));
            unimproved{pop} = 0;                         % Reset counter
            kingOfTheHill{pop} = fitOrganism{pop};	% Save the top organism
        else
            unimproved{pop} = unimproved{pop} + 1;        % Increment counter
        end

        % Record training progress
        optimal_fitness{pop} = [optimal_fitness{pop} fit_raw(fit_rank{pop}(1))];
        average_fitness{pop} = [average_fitness{pop} mean(fit_raw)];
        diversity{pop} = [diversity{pop} mean(similarity{pop}(similarity{pop} ~= 0))];
        
        % Update parameters
        generation{pop} = new_generation;
    
    end
    
    % ---------------------------------------------------------------------
    % Plot Results
    % ---------------------------------------------------------------------
    
    for pop = 1:Pp
        
        % Display heatmap of fitness
        fitMatrix = reshape(fitness{pop}, [sqrt(2^R{pop}), sqrt(2^R{pop})]);
        figure(pop); subplot(2,3,2); hold on
        imagesc(fitMatrix); title("Fitness at Loop " + string(loopCounter));
        xlim([0,sqrt(2^R{pop})+1]); ylim([0,sqrt(2^R{pop})+1]); colorbar; hold off
        
        % Plot organism with highest fitness
        figure(pop); subplot(2,3,4);
        fitOrganism{pop} = ...
            smart.general.plotCell(generation{pop}{fit_rank{pop}(1)}, false);
        for i = 1:numChromosomes
            area(fitOrganism{pop}(:,i), 'FaceAlpha', 0.3)
            hold on
        end; xlim([0, length(mixedVectors{pop})])
        title("Fittest Organism - Generation " + string(loopCounter))
        xlabel('Timestep'); ylabel('Amplitude'); grid on; hold off
        
        % Plot sampling frequency per parent
        samplingMatrix = reshape(samples{pop}, [sqrt(2^R{pop}), sqrt(2^R{pop})]);
        figure(pop); subplot(2,3,5); hold on
        imagesc(samplingMatrix); title("Sampling Frequency of Parents");
        xlim([0,sqrt(2^R{pop})+1]); ylim([0,sqrt(2^R{pop})+1]); colorbar; hold off
        
        % Plot training progress
        figure(pop); subplot(2,3,3)
        max_fitness = fitness_op{pop}*ones(length(optimal_fitness{pop}),1);
        area(0:loopCounter, optimal_fitness{pop}, 'FaceAlpha', 0.3); hold on; grid on
        area(0:loopCounter, average_fitness{pop}, 'FaceAlpha', 0.3)
        plot(0:loopCounter, max_fitness, '--k')
        title("Population Fitness for Generation " + string(loopCounter) + "/" + ...
            string(maxGen) + " - Patience Counter: " + string(max(patience-unimproved{pop},0)))
        legend("Optimal Fitness", "Overall Fitness", "Maximum Fitness"); 
        xlabel("Generation"); ylabel("Fitness")
        xlim([0, length(optimal_fitness{pop})-1]); hold off
        
        % Plot change in population diversity
        figure(pop); subplot(2,3,6);
        area(0:loopCounter-1, diversity{pop}, 'FaceAlpha', 0.3)
        hold on; grid on; title("Population Diversity for Generation " + string(loopCounter))
        xlabel("Generation"); ylabel("Diversity");
        if loopCounter > 1; xlim([0, length(diversity{pop})-1]); end; hold off
        
    end
    
    % ---------------------------------------------------------------------
    % Part 7: Migration
    % ---------------------------------------------------------------------
    
    for pop = 1:Pp
        % Arguments: old_populations, pop, Mi, Plot
        generation = smart.genetic.migrate(generation, pop, Mi{pop}, Plot);

        % Re-calculate fitness of population
        [~, fitness{pop}, fit_rank{pop}] = ...
            smart.genetic.getFitness(generation{pop}, P{pop}, F{pop});
    end
    
    % ---------------------------------------------------------------------
    % Part 8: Pandemic
    % ---------------------------------------------------------------------
    
    % Arguments: population, penetration, vLoad, gen, pop, Plot
    % TODO: Conditional release of pandemic
    for pop = 1:Pp
        if ~strcmp(V{pop}, 'none') && mod(loopCounter, V{pop}) == 0
            % Increment virus counter
            virusCounter{pop} = virusCounter{pop} + 1;
            % Introduce viral vector
            [generation{pop}, ViralSig{pop}{virusCounter{pop}}, Fig{pop}(virusCounter{pop})] = ...
                smart.genetic.infect(generation{pop}, 0.05, 0.05, loopCounter, pop, Plot);
        end

        % Post-infection monitoring of viral signature
        % Post-monitoring once virus has been introduced
        if virusCounter{pop} > 0 && Plot
            % For each population's viruses...
            for virus = 1:virusCounter{pop}
                
                % Set placeholder for vector presence in population
                infectionFLAG = zeros(length(generation{pop}), 1);
                % Inspect each organism for allelic combination at each locus
                for org = 1:length(generation{pop})
                   % Set infection flag if all allelic combinations are present
                   if isequal(generation{pop}{org}(ViralSig{pop}{virus}.loci, :), ...
                           ViralSig{pop}{virus}.alleles)
                       infectionFLAG(org) = 1;  % Set the infection flag
                   end
                end

                % Plot results
                figure(Fig{pop}(virus)); subplot(1,2,1); 
                imagesc(reshape(infectionFLAG, ...
                    [sqrt(length(generation{pop})), sqrt(length(generation{pop}))]))
                hold on; colorbar; 
                xlim([0, sqrt(length(generation{pop}))+1])
                ylim([0, sqrt(length(generation{pop}))+1])
            
            end
        end
    end
    
    % Update loop flag in response to stop criteria
    if loopCounter > maxGen || length(find([unimproved{:}] > patience)) == Pp
        loopFLAG = false;
    end
    
    % Increment loop counter
    loopCounter = loopCounter + 1;
    
end

% -------------------------------------------------------------------------
% Part 9: Format Data for Completion
% -------------------------------------------------------------------------

% Return sorted vector
[~, koth] = max([fitness_op{:}]);    % Get the population with the KotH
sortedVectors = kingOfTheHill{koth};

end

