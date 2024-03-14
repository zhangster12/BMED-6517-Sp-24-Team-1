# Cardio (ECG)

Functions for processing electrocardiogram data.

## Table of Contents

1. [Functions](#functions)
2. [Example](#example)

## Functions

| Function | Purpose |
| --- | --- |
| `automatedHR()` | Computes heart rate (HR) with automated R-peak threshold calibration |
| `cleanECG()` | Cleaning noisy ECG signals (used with `Oink`) | 
| `computeHR()` | Computes HR given an R-peak threshold |
| `computeHRV()` | Compute heart rate variability (HRV) via several methods |
| `slidingWindowHRV()` | Computes HRV over a recording session using a sliding window | 

## Example

Cardiac signals synced to a concurrent ECG may be segmented with a simple automated script. Note that in this script, the user sets the parameters `minDist` and `maxLen` to indicate the minimum distance between R-peaks and maximum length of the signal segments.

```matlab
% Compute threshold for R-peak detection
[~, ~, threshold] = cardio.ecg.automatedHR(ecg, <varargin>);
% Segment ECG and concurrent signal based on R-peaks
[segmented_ecg, segmented_signal] = cardio.general.separateBeats(signal, 'ecg', ecg, <varargin>);
```
