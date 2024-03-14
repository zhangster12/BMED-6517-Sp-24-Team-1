# Cardio (SQI)

Sub-package for signal quality index development and implementation

## Table of Contents

1. [Summary](#summary)
2. [Functions](#functions)

## Summary

This subpackage contains code used in the 2019 paper "A Unified Framework for Quality Indexing and Classification of Seismocardiogram Signals" in [IEEE JBHI](https://ieeexplore.ieee.org/document/8777167). These functions were used to analyze the response of the signal quality index (SQI) using different distance metrics. The function `cardio.sqi.filter()` computes the SQI as defined in the paper. 

## Functions

| Function | Purpose |
| --- | --- |
| `decompose()` | Decomposes a signal using a windowed Fourier series (STFS) |
| `filter()` | Computes the SQI |
| `getCoeffs()` | Returns Fourier series coefficients |
| `impulse()` | Models the impulse response of distance filters for SCG |
| `peakMatch()` | Implements DTW optimized for peak-matching |
| `reconstruct()` | Reconstructs a signal from its windowed Fourier series |
| `response()` | Computes changes in distance given windowed sinusoidal noise |
| `synthetic()` | Generates a synthetic SCG signal based on a template |
| `warpFeatures()` | Extracts features using dynamic time warping (DTW) and Fourier series subspaces |
| `warpFilter()` | Uses DTW to return a subset of representative signals |
| `warpGenerate()` | Generate a template from a randomly-generate DTW warp path |
| `warpMask()` | Computes the SQI given a template and Fourier series coefficient mask |
| `wasserstein()` | Computes the wasserstein (earth-mover) distance between two signals |
