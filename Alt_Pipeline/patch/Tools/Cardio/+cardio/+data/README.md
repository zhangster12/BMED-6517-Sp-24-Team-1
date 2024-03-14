# Cardio (Data)

Extracting, and formatting data for use with the Cardio package.

## Table of Contents

1. [Summary](#summary)
2. [Exercise Study](#exercise-study)
3. [Placement Study](#placement-study)
4. [Stress Study](#stress-study)

## Summary

This is a MATLAB package for processing cardiovascular signals of the following modalities: EEG, ICG, SCG, and PPG. This README file contains general information for using the package; see each sub-folder for a sub-package-specific README file.

## Exercise Study

**Associated Papers:** [Shandhi 2019](https://ieeexplore.ieee.org/document/8627923), [Hersek 2019](https://ieeexplore.ieee.org/abstract/document/8779609)
**Path:** `+cardio/+data/exerciseStudy.m`

> This dataset includes ECG, ICG, and SCG (accelerometer and gyroscope) data recorded during a rest -> exercise -> recovery protocol (26 subjects). The `Intervals` enumeration specifies the different activity levels.

## Placement Study

**Associated Paper:** [Ashouri 2017](https://ieeexplore.ieee.org/abstract/document/7919162)
**Path:** `+cardio/+data/placementStudy.m`

> This dataset includes ECG, ICG, and SCG (accelerometer) data recorded during a rest -> exercise -> recovery protocol with sensors in five different sensor locations. The `Positions` enumeration specifies the different activity levels.

## Stress Study

**Associated Paper:** [Gurel 2018](https://ieeexplore.ieee.org/document/8478178)
**Path:** `+cardio/+data/stressStudy.m`

> This dataset includes ECG, ICG, PPG, and SCG (accelerometer) data recorded during a variety of mental tasks designed to induce stress.
