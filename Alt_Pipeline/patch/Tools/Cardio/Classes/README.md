# Cardio (Classes)
Supporting classes for the Cardio package.

## Summary

This folder contains various enumerations which may be used with functions in the Cardio package, as well as classes for creating useful objects.

## Table of Contents

1. [Classes](#classes)
2. [Enumerations](#enumerations)

## Classes

### Oink

**Path:** `Cardio/Classes/@Oink`

The `Oink` class extracts, organizes, and stores data collected as part of the Hypovolemia Study. To use the `Oink` class, first instantiate an object as follows, where `subjects` is an integer lists of subjects to process.

```matlab
oink = Oink(subjects, "Description");
```

To extract data, use the member function `oink.extractData()`, where `levels` is a vector of type [`Level`](#level) of blood volume levels to extract. For additional options, type `help Oink.extractData` in the command window.

```matlab
oink.extractData('levels', levels);
```

Within the oink object, data can be accessed as `oink.dataset{<subject>}.<level>.<modality>`. From here, the member functions of the object can be used to compute hemodynamic parameters from the dataset and organize the features into formats for training and testing analytical models. The enumerations [`Feature`](#feature), [`Label`](#label), [`Mode`](#mode), and [`Sample`](#sample) are useful when interfacing with this class, and are described below.

Though data from multiple subjects may be combined in a single object, this is not recommended due to the large file size required for this class. Therefore, each object should typically contain data from one subject. The member functions of the `Oink` class are as follows.

| Function | Purpose |
| --- | --- |
| `Oink()` | Class constructor |
| `Oink.amplitudes()` | Computes mean amplitudes for beat-separated waveforms |
| `Oink.createTemplateSet()` | Creates TemplateSet objects for SQI generation |
| `Oink.ecgHRV()` | Computes heart rate variability using ECG R-peak indices |
| `Oink.extractData()` | Extracts, processes, and organizes the raw data |
| `Oink.getParent()` | Returns the parent `Level` for the specified sub-level |
| `Oink.getSubinterval()` | Pulls data from a sub-interval of `Level.allAbsolute` or `Level.allRelative` |
| `Oink.sampleToBeat()` | Converts sample indices to beat indices |
| `Oink.scgAortic()` | Computes the rAO/C interval from SCG waveforms |
| `Oink.scgPTT()` | Computes the PAT/PTT from SCG and PPG waveforms |
| `Oink.score()` | Generates SQI scores using requisite TemplateSet objects |
| `Oink.timeToSample()` | Converts military time to sample number |
| `Oink.trainingData()` | Creates a feature matrix and label vector from the data |
| `Oink.trueAortic()` | Computes the rAO/C interval from pressure waveforms |
| `Oink.trueMAP()` | Computes the mean arterial pressure from pressure waveforms |
| `Oink.truePTT()` | Computes the PAT/PTT from aortic and femoral pressure waveforms |

## Enumerations

### Feature

**Path:** `Cardio/Classes/Feature.m`

Enumeration describing feature types for the Hypovolemia Study dataset:

| Value | Meaning |
| --- | --- |
| scgPEP | SCG-derived pre-ejection period (PEP)  |
| scgLVET | SCG-derived left ventricular ejection time (LVET) |
| truePEP | Aortic pressure-derived PEP |
| trueLVET | Aortic pressure-derived LVET |
| scgPAT | PPG-dervied pulse arrival time (PAT) |
| scgPTT | SCG/PPG-derived pulse transit time (PTT) |
| truePAT | Femoral pressure-derived PAT |
| truePTT | Aortic/Femoral pressure-derived PTT |
| hr | Heart rate (Biopac) |
| hrvDifference | Heart rate variability (HRV) using difference method |
| hrvSpectral | HRV using spectral method |
| hrvPoincare | HRV using Poincare method |
| scgAmplitudeApex | Amplitude of SCG signal (z-axis) at apex |
| scgAmplitudeSternum | Amplitude of SCG signal (z-axis) at sternum |
| ppgAmplitudeProximal | PPG amplitude at apex |
| ppgAmplitudeDistal | PPG amplitude at femoral artery |
| scgPEPOverLVET | SCG-derived PEP/LVET |
| truePEPOverLVET | Aortic pressure-derived PEP/LVET |
| proximalSPO2 | Oxygen saturation (SPO2) using apex PPG |
| distalSPO2 | SPO2 using femoral artery PPG |
| scgSQI | Signal quality index (SQI) of SCG signal (z-axis) at sternum |
| ppgSQI | SQI of PPG signal at femoral artery |
| meanAorticPressure | Mean aortic pressure |
| meanFemoralPressure | Mean femoral pressure |
| meanWedgePressure | Mean pulmonary capillary wedge pressure |
| meanRAPressure | Mean right atrial pressure |
| aorticPressureSQI | SQI for aortic pressure waveform |
| femoralPressureSQI | SQI for femoral pressure waveform |
| all | The set of all features |

### Index

**Path:** `Cardio/Classes/Index.m`

Enumeration describing types of distance metrics:

| Value | Meaning |
| --- | --- |
| FMDTW | Dynamic-Time Feature Matching |
| DTW | Dynamic Time Warping |
| Model | Trained ML Model |
| Norm | L2-norm |

### Intervals

**Path:** `Cardio/Classes/Intervals.m`

Enumeration describing valid intervals for the Exercise Study dataset:

| Value | Meaning |
| --- | --- |
| Rest | Resting Interval |
| Recovery | Recovery Interval |
| Exercise | Exercise Interval |
| Squat | Squatting Interval |
| Walk | Walking Interval |
| All | All Intervals |

### Label

**Path:** `Cardio/Classes/Label.m`

Enumeration describing labeling methods for Hypovolemia Study dataset:

| Value | Meaning |
| --- | --- |
| bloodVolume | Raw blood volume levels |
| relativeVolume | Adjusted blood volume levels (w.r.t. estimated compensatory reserve) |

### Level

**Path:** `Cardio/Classes/Level.m`

Enumeration describing blood volume levels for the Hypovolemia Study dataset:

| Value | Meaning |
| --- | --- |
| all | All datapoints |
| allRelative | All relative hypovolemia datapoints |
| allAbsolute | All absolute hypovolemia datapoints |
| relBaseline1 | First baseline during relative hypovolemia | 
| relative5 | 5% drop in BP during relative hpovolemia |
| relative10 | 10% drop in BP during relative hypovolemia |
| relative20 | 20% drop in BP during relative hypovolemia |
| relBaseline2 | Second baseline during relative hypovolemia |
| absBasleine1 | First baseline during absolute hypovolemia | 
| absDecrease7 | 7% absolute blood volume loss |
| absDecrease14 | 14% absolute blood volume loss |
| absDecrease21 | 21% absolute blood volume loss |
| absDecrease28 | 28% absolute blood volume loss |
| absIncrease28 | 28% blood volume loss to 21% blood volume loss |
| absIncrease21 | 21% blood volume loss to 14% blood volume loss |
| absIncrease14 | 14% blood volume loss to 7% blood volume loss |
| absIncrease7 | 7% blood volume loss to 0% blood volume loss |
| absBaseline2 | Second baseline during absolute hypovolemia |

### Mode

**Path:** `Cardio/Classes/Mode.m`

Enumeration describing sensor modalities for the Hypovolemia Study dataset:

| Value | Meaning |
| --- | --- |
| biopacECG | Biopac ECG  |
| t3ECG | T3 ECG |
| sternumSCG_x | Sternum SCG (x-axis) |
| sternumSCG_y | Sternum SCG (y-axis) |
| sternumSCG_z | Sternum SCG (z_axis) |
| apexSCG_x | Apex SCG (x-axis) |
| apexSCG_y | Apex SCG (y-axis) |
| apexSCG_z | Apex SCG (z-axis) |
| apexSPO2 | PPG-derived SPO2 (apex) |
| femoralSPO2 | PPG-derived SPO2 (femoral) |
| apexPPG | Apex PPG |
| femoralPPG | Femoral PPG |
| aorticPressure | Aortic pressure waveform (catheter) |
| femoralPressure | Femoral pressure waveform (catheter) |
| wedgePressure | Pulmonary capillary wedge pressure (catheter) |
| rightAtrialPressure | Right atrial pressure (catheter) |

### Positions

**Path:** `Cardio/Classes/Positions.m`

Enumeration describing valid positions for the Placement Study dataset:

| Value | Meaning |
| --- | --- |
| Center | Mid-Sternum Position |
| Above | Upper-Sternum Position |
| Below | Lower-Sternum Position |
| Left | Left of Center Position |
| Right | Right of Center Position |
| All | All Positions |

### Sample

Enumeration describing sampling method for [`Oink.trainingData()`](#oink):

**Path:** `Cardio/Classes/Sample.m`

| Value | Meaning |
| --- | --- |
| allSamples | Extract all datapoints |
| fixedInterval | Extract datapoints at fixed intervals |
| randomized | Extract datapoints at random |
