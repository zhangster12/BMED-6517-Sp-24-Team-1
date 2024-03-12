%% HRV Extraction over Two-Minute Intervals: Demo/Example

close all;
clear all;
clc;

%% Extracting all HRV Features
% Load an RRintervals file. The first column will be the times, second
% column will be the RRinterval values
filename = 'RRintervals/sub3128/stim0_RRintervals.csv';
RRintervals = readmatrix(filename);

Fs = 500; % sampling freq. of ECG in Hz

% Choose a 2-min interval over which we will compute HRV features
RRseg = RRintervals(find(80 <= RRintervals & RRintervals <= 200), :);

% Run RRtoHRVDemo. This will output a struct with 26 total features. Refer
% to "RRtoHRVDemo.m" for the names of the different features
HRV_feat = RRtoHRVDemo(RRseg(:, 2), Fs);

%% HRV Feature Consolidation
% Here, I will consolidate HRV_feat into a set of features that are
% relevant for stress detection

% Time-Domain HRV 
SDNN = HRV_feat(9);
RMMSD = HRV_feat(10);
pnn50 = HRV_feat(11);
nn50 = HRV_feat(12);

% Nonlinear HRV
SD1 = HRV_feat(13);
SD2 = HRV_feat(14);
SDratio = HRV_feat(15);

% Freq-Domain HRV
VLF = HRV_feat(20);
LF = HRV_feat(21);
HF = HRV_feat(22);
LFHF = HRV_feat(23);
LFn = HRV_feat(24);
HFn = HRV_feat(25);
POW = HRV_feat(26);