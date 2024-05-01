# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 22:49:51 2024

@author: jha91
"""

## Parameter initialization (clear all)
from IPython import get_ipython
#get_ipython().magic('reset -sf')

Dir1 = 'Consolidated_Features'
import os
import pandas as pd
import scipy.io as sio
import numpy as np
from scipy.signal import butter, filtfilt, resample 

Subject = os.listdir(Dir1) # Post_MI: 35, Healthy: 9
NumSub = np.size(Subject)
subj_list = [
    "3128", "3129", "3130", "3131", "3132", "3133", "3136", "3137", "3138", "3139", 
    "3140", "3141", "3142", "3143", "3147", "3148", "3149", "3150", "3151", "3152", 
    "3153", "3154", "3155", "3156", "3158", "3159", "3160", "3162", "6037", "6038", 
    "6043", "6044", "6045", "6046", "6047", "6048", "6049"]

NumStim = 5
Stimulus = ['Rest', 'Reading', 'SpeechPrep', 'Speech', 'Recovery']
NumFeat = 5
Feature = ['HR', 'PAT', 'PEP', 'PPGamp', 'PTTrecip']

# 1st column: time, second column: feature values
Subj = np.zeros(len(Subject))

for i in range(len(Subject)):
    Subj[i] = int(Subject[i][3:])

Subj.sort()

dataframes_MI = []
dataframes_Ht = []
sub_MI = -1
sub_Ht = -1

for sub in range(NumSub):

    if str(round(Subj[sub])) in subj_list:

        if str(Subj[sub]).startswith('3'):
            sub_MI = sub_MI + 1
            dataframes_MI.append([]) 

            for stim in range(NumStim):
                dataframes_MI[sub_MI].append([])

                for feat in range(NumFeat):
                    dataframes_MI[sub_MI][stim].append([]) 
                    Feat_load = os.path.join(Dir1, 'sub' + str(int(Subj[sub])), 'stim' + str(stim) + '_' + Feature[feat] + '.csv')
                    data = pd.read_csv(Feat_load)
                    dataframes_MI[sub_MI][stim][feat] = data.values

        if str(Subj[sub]).startswith('6'):
            sub_Ht = sub_Ht + 1
            dataframes_Ht.append([])

            for stim in range(NumStim):
                dataframes_Ht[sub_Ht].append([])
                
                for feat in range(NumFeat):
                    dataframes_Ht[sub_Ht][stim].append([])
                    Feat_load = os.path.join(Dir1, 'sub' + str(int(Subj[sub])), 'stim' + str(stim) + '_' + Feature[feat] + '.csv')
                    data = pd.read_csv(Feat_load)
                    dataframes_Ht[sub_Ht][stim][feat] = data.values

# %% Data visualize
import matplotlib.pyplot as plt

def access_df(df, participant_index, stimulus_index, feature_index):
    data = df[participant_index][stimulus_index][feature_index]
    time = data[:, 0]
    feat_values = data[:, 1]
    return time, feat_values

def join_stimulus(df, participant_index, feature_index):
    time = []
    feat_values = []

    for stimulus_index in range(NumStim):
        data = df[participant_index][stimulus_index][feature_index]
        time_temp = data[:, 0]
        values = data[:, 1]
        time.append(time_temp)
        feat_values.append(values)
    
    time = np.concatenate(time)
    feat_values = np.concatenate(feat_values)

    return time, feat_values

# Compare healthy and unhealthy
fig, axes = plt.subplots(len(Feature), 2, sharex=True, sharey='row')

axes[0, 0].set_title('Healthy')
axes[0, 1].set_title('Post-MI')

for i, feature in enumerate(Feature):
    time_Ht, feat_Ht = join_stimulus(dataframes_Ht, 1, i)
    axes[i, 0].plot(time_Ht, feat_Ht)
    axes[i, 0].set_ylabel(feature)

    time_MI, feat_MI = join_stimulus(dataframes_MI, 3, i)
    axes[i, 1].plot(time_MI, feat_MI)

fig.supxlabel('Time (s)')
plt.tight_layout()
plt.savefig('Figure_1.png')
plt.show()

# Individual stimulus

fig, axes = plt.subplots(len(Feature), len(Stimulus), sharex='col', sharey='row')

for i, feature in enumerate(Feature):
    axes[i, 0].set_ylabel(feature)

    for j, stim in enumerate(Stimulus):
        time_Ht, feat_Ht = access_df(dataframes_Ht, 1, j, i)
        time_MI, feat_MI = access_df(dataframes_MI, 3, j, i)

        axes[i, j].plot(time_Ht, feat_Ht)
        axes[i, j].plot(time_MI, feat_MI)

        axes[len(Feature) - 1, j].set_xlabel(stim)

fig.supxlabel('Time (s)')

fig.legend(['Healthy', 'Post-MI'], loc='upper center', ncols=2, framealpha=1)

plt.tight_layout()
plt.savefig('Figure_2.png')
plt.show()
