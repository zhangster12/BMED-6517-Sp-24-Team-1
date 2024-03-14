% -------------------------------------------------------------------------
% Genetic Algorithm Implementation Script
% Created By: Jonathan Zia
% Date Created: 12/28/2018
% Inan Research Lab @ Gatech
% -------------------------------------------------------------------------

close all; clear; clc;

% Add library to path
addpath('/Users/jonathanzia/Github/Smart')

% Generate data
n = 500; m = 4; mixedVectors = zeros(n,m); % Set parameters
for i = 1:m
    defaultVal = i; % Set default value for vector;
    mixedVectors(:,i) = defaultVal;
end

% Introduce errors
mixedArray1 = smart.genetic.generateAlleles(mixedVectors, 10);
mutatedArray1 = smart.genetic.shuffleAlleles(mixedArray1);
mixedVectors1 = smart.general.plotCell(mutatedArray1, false);
mixedArray2 = smart.genetic.generateAlleles(mixedVectors, 10);
mutatedArray2 = smart.genetic.shuffleAlleles(mixedArray2);
mixedVectors2 = smart.general.plotCell(mutatedArray2, false);

% Create Population structs
% Arguments: mixedVectors, B, R, Mu, C, L, P, N, Mi, F, V
Population1 = ...
    smart.genetic.definePopulation(mixedVectors1, 10, 10, 10, 10, 10, 4, 0, 0.02, 'avgstd', 3);
Population2 = ...
    smart.genetic.definePopulation(mixedVectors2, 10, 10, 10, 10, 10, 4, 0, 0.02, 'maxstd', 'none');

% Call function
% Arguments: patience, maxGen, Population1, Population2, ...
sortedVectors = smart.genetic.geneticSort(5, 100, Population1, Population2);