# Smart (Genetic)

Functions for implementing genetic algorithms.

## Table of Contents

1. [Summary](#summary)
2. [Functions](#functions)

## Summary

This sub-package implements genetic algorithms to sort data that may have been scrambled among different vectors. An example implementation of genetic algorithms is shown in `Smart/+tutorials/genetic.m`. For additional help for each of the following functions, type `help Smart.genetic.<function>` in the command window. Most functions are called by `smart.genetic.geneticSort()`.

## Functions

| Function | Purpose |
| --- | --- |
| `definePopulation()` | Define a population for use in `geneticSort()` |
| `generateAlleles()` | Divide a vector into a cell vector (alleles) |
| `geneticSort()` | Sort a set of scrambled vectors using genetic algorithms |
| `getFitness()` | Obtain the fitness for each organism |
| `infect()` | Introduce genetic variability via viral vectors |
| `migrate()` | Perform genetic crossover between two populations |
| `mutate()` | Introduce genetic variability via mutation |
| `produceGeneration()` | Produce a generation of offspring |
| `recombine()` | Perform genetic recombination |
| `shuffleAlleles()` | Generate a new organism via genetic shuffling | 
