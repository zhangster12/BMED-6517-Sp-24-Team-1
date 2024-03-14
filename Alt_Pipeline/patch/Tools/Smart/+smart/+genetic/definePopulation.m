function Population = ...
    definePopulation(mixedVectors, B, R, Mu, C, L, P, N, Mi, F, V)

% -------------------------------------------------------------------------
% This function creates a Population struct from a set of input parameters.
% These structs may be input as arguments to geneticSort.m. Input
% parameters include the following:
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
% - F               HYPERPARAMETER: Method of calculating fitness
%                   - 'avgstd': Mean standard deviation
%                   - 'maxstd': Maximum standard deviation
% - V               HYPERPARAMETER: Number of generations between viral pandemics
%                   - 'none': No viral pandemics
% -------------------------------------------------------------------------

% Create Population struct
Population.mixedVectors = mixedVectors; Population.B = B; Population.R = R; 
Population.Mu = Mu; Population.C = C; Population.L = L; Population.P = P;
Population.N = N; Population.Mi = Mi; Population.F = F; Population.V = V;

end

