# Cardio (General)

General scripts for cardiovascular signal processing.

## Functions

| Function | Purpose |
| --- | --- |
| `alignSignals()` | Aligns two signals using cross-correlation |
| `alignSignalSet()` | Align two separate sets of signals synched to two separate ECGs  |
| `alignVectors()` | Align two vectors that have had indices removed |
| `confusion()` | Create a confusion matrix |
| `createBPF()` | Create a band-pass filter |
| `createHPF()` | Create a high-pass filter |
| `createLPF()` | Create a low-pass filter |
| `detrendFeatures()` | Detrends feature vectors |
| `dynamicStats()` | Fits polynomial to feature vector |
| `ema()` | Computes the exponential moving average (EMA) |
| `ensembleAvg()` | Computes the ensemble average |
| `getOutliers()` | Returns the outlier indices for each feature vector |
| `getPeaks()` | Returns the peaks, valleys, and inflections in a data vector |
| `heartbeatMovie()` | Creates an animation of a heartbeat-separated signal |
| `labelFeatures()` | Manual label features in time-series data |
| `labelSegmentFeatures()` | Label features in a single signal segment |
| `map()` | Maps indices from a signal into its warped signal, given a warp path |
| `parsave()` | Saving workspace variables in parfor loops |
| `project()` | Map scatterpoints to the nearest point on a circle |
| `pulseTime()` | Infers the arrival time of pulsatile signals |
| `removeRegion()` | Remove signal segments based on `imagesc` plot of signal segments |
| `rollingFilter()` | Filter a signal using a rolling window (mean and std. dev.) |
| `rsquared()` | Computes the coefficient of determination between two signals |
| `selectTemplate()` | Select a template by finding the most representative signal |
| `separateBeats()` | Performs beat-separation of a signal based on ECG |
| `slideWindow()` | Creates a sliding window of a 1D input signal |
| `slidingMove()` | Visualize sets of heartbeat-separated signals |
| `template()` | Generates template signals using Woody's method |

## Scripts

| Script | Purpose |
| --- | --- |
| `visual_scoring.mlx` | Rapidly generate SCG signals and score quality |
