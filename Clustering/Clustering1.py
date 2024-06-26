# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 22:49:51 2024

@author: jha91
"""

## Parameter initialization (clear all)
from IPython import get_ipython
get_ipython().magic('reset -sf')

# Dir1 = r'C:\Users\jha91\OneDrive - Georgia Institute of Technology\Year1 2024 Spring\BMED 6517 ML Bio\Group project\Data\Consolidated_Features'
Dir1 = '..\Consolidated_Features'
import os
import pandas as pd
import numpy as np
import seaborn as sns
import scipy.io as sio
from scipy import stats
from scipy.signal import butter, filtfilt, resample 

Subject = os.listdir(Dir1) # Post_MI: 28, Healthy: 9
NumSub = np.size(Subject)
subj_list = [
    "3128", "3129", "3130", "3131", "3132", "3133", "3136", "3137", "3138", "3139", 
    "3140", "3141", "3142", "3143", "3147", "3148", "3149", "3150", "3151", "3152", 
    "3153", "3154", "3155", "3156", "3158", "3159", "3160", "3162", "6037", "6038", 
    "6043", "6044", "6045", "6046", "6047", "6048", "6049"]
NumStim = 5
Stimulus = ['Rest', 'Reading', 'SpeechPrep', 'Speech', 'Recovery'] # rest vs. speechprep -> stim0 vs. stim2
NumFeat = 5
Feature = ['HR', 'PAT', 'PEP', 'PPGamp', 'PTTrecip']
# 1st column: time, second column: feature values
Subj = np.zeros(len(Subject))
for i in range(len(Subject)):
    Subj[i] = int(Subject[i][3:])  
Subj.sort()


dataframes_MI = []
dataframes_Ht = []
dataframes_all = []
sub_MI = -1
sub_Ht = -1
sub_all = -1
for sub in range(NumSub):
    if str(round(Subj[sub])) in subj_list:
        if str(Subj[sub]).startswith('3'):
            sub_MI = sub_MI+1;
            dataframes_MI.append([]) 
            for stim in range(NumStim):
                dataframes_MI[sub_MI].append([]) 
                for feat in range(NumFeat):
                    dataframes_MI[sub_MI][stim].append([]) 
                    Feat_load = os.path.join(Dir1, 'sub' + str(int(Subj[sub])), 'stim' + str(stim) + '_' + Feature[feat] + '.csv');
                    data = pd.read_csv(Feat_load)
                    dataframes_MI[sub_MI][stim][feat] = data.values
                    # dataframes_MI[sub_Ht][stim][feat].append(np.zeros(len(data.values[:,2])))
        if str(Subj[sub]).startswith('6'):
            sub_Ht = sub_Ht+1;
            dataframes_Ht.append([])
            for stim in range(NumStim):
                dataframes_Ht[sub_Ht].append([])
                for feat in range(NumFeat):
                    dataframes_Ht[sub_Ht][stim].append([])
                    Feat_load = os.path.join(Dir1, 'sub' + str(int(Subj[sub])), 'stim' + str(stim) + '_' + Feature[feat] + '.csv');
                    data = pd.read_csv(Feat_load)
                    dataframes_Ht[sub_Ht][stim][feat] = data.values
                    # dataframes_Ht[sub_Ht][stim][feat].append(np.zeros(len(data.values[:,2])))
                    
for sub in range(NumSub):
    if str(round(Subj[sub])) in subj_list:
        sub_all = sub_all+1;
        dataframes_all.append([]) 
        for stim in range(NumStim):
            dataframes_all[sub_all].append([]) 
            for feat in range(NumFeat):
                dataframes_all[sub_all][stim].append([]) 
                Feat_load = os.path.join(Dir1, 'sub' + str(int(Subj[sub])), 'stim' + str(stim) + '_' + Feature[feat] + '.csv');
                data = pd.read_csv(Feat_load)
                dataframes_all[sub_all][stim][feat] = data.values

# %% Data visualize - pairplot using seaborn
# % All

NumSub_all = len(dataframes_all)
feat_values = np.zeros((NumSub, NumFeat)) #np.zeros((NumSub_MI,1))
feat_values_total = []
y_val_total = []
all_feat_avg_total = []

for stim in range(NumStim):
    for feat in range(NumFeat):
        for sub in range(NumSub_all):
            feat_values[sub, feat] = np.mean(dataframes_all[sub][stim][feat][:,1])
        y_val = np.ones(len(feat_values))*stim
        if len(feat_values_total) < NumSub_all:
            feat_values_total = feat_values[:, feat:feat+1]
        else:
            feat_values_total = np.concatenate((feat_values_total, feat_values[:, feat:feat+1]), axis = 1)
    y_val_total = np.concatenate((y_val_total, y_val))
    if len(all_feat_avg_total) < NumSub_all:
        all_feat_avg_total = feat_values_total
    else:
        all_feat_avg_total = np.concatenate((all_feat_avg_total, feat_values_total))
    feat_values_total = []

# Make it into pandas form

# Sort states: rest (0) & stress (2)
target_values = [0, 2]

# Use list comprehension to get the indices of elements with target values
indices = [i for i, value in enumerate(y_val_total) if value in target_values]
y_val_total_1 = y_val_total[indices]
all_feat_avg_total_1 = all_feat_avg_total[indices,:]
  
# dictionary of lists 
dict = {'States': y_val_total_1, Feature[0]: all_feat_avg_total_1[:,0], 
        Feature[1]: all_feat_avg_total_1[:,1],
        Feature[2]: all_feat_avg_total_1[:,2],
        Feature[3]: all_feat_avg_total_1[:,3],
        Feature[4]: all_feat_avg_total_1[:,4]} 
    
df_all = pd.DataFrame(dict) 
# print(df)
sns.pairplot(df_all, hue="States")


# %% Baseline correction (normalization)
mean_all = []
data_all = []
data_y_all = []
num_data = 80
dp = round(num_data/2)

for subj in range(len(subj_list)):
    data_all.append([]) 
    data_y_all.append([]) 
    for stim in range(2):
        data_all[subj].append([]) 
        data_y_all[subj].append([]) 
        for feat in range(NumFeat):
            data_all[subj][stim].append([]) 
            data_y_all[subj][stim].append([]) 
            # sort 200 points each for state 0 (baseline) and state 2 (feature)
            baseline = dataframes_all[subj][0][feat][:, 1]
            baseline_mid = baseline[round(len(baseline)/2)-dp:round(len(baseline)/2)+dp] 
            feature = dataframes_all[subj][2][feat][:, 1]
            feature_mid = feature[round(len(feature)/2)-dp:round(len(feature)/2)+dp]
            mean_all.append((np.mean(feature_mid) - np.mean(baseline_mid))/np.mean(baseline_mid))
            if stim == 0:
                data_all[subj][stim][feat] = (baseline_mid - np.mean(baseline_mid))/np.mean(baseline_mid)
                data_y_all[subj][stim][feat] = np.zeros([len(data_all[subj][stim][feat]),1])
            if stim == 1:
                data_all[subj][stim][feat] = (feature_mid - np.mean(baseline_mid))/np.mean(baseline_mid)
                data_y_all[subj][stim][feat] = np.ones([len(data_all[subj][stim][feat]),1])

# reshape [subj][stim][feat] -> [stim][feat][:]
data_all_1 = []
data_y_all_1 = []

for feat in range(NumFeat):
    temp_data_stim = []
    temp_labels_stim = []
    for stim in range(2):
        temp_data = []
        temp_labels = []
        for subj in range(len(subj_list)):
            temp_data.extend(data_all[subj][stim][feat])
            temp_labels.extend(data_y_all[subj][stim][feat])
        temp_data_stim.extend(temp_data)
        temp_labels_stim.extend(temp_labels)
    data_all_1.append(temp_data_stim)
    data_y_all_1.append(temp_labels_stim)
    
# %% list into array

min_length = min(len(lst) for lst in data_all_1) # Ensure all lists have the same length (2480)
data_all_1_truncated = [lst[:min_length] for lst in data_all_1]
data_all_1_array = np.array(data_all_1_truncated) # Convert the truncated list to a numpy array
data_all_1_array = data_all_1_array.reshape(-1, 5) # Reshape the array to have 5 columns
# print(data_all_1_array.shape) 

# same for y
min_length_y = min(len(lst) for lst in data_y_all_1)
data_y_all_1_truncated = [lst[:min_length_y] for lst in data_y_all_1]
data_y_all_1_array = np.array(data_y_all_1_truncated)
data_y_all_1_array = data_y_all_1_array.reshape(-1, 5)

# %% Dimensionality reduction using PCA
import numpy as np
from sklearn.datasets import load_digits
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import load_wine
from sklearn.decomposition import PCA

feats = data_all_1_array[:,:]
X = feats
labels = data_y_all_1_array[:,0]
y = labels

# feats_upc = PCA(n_components=2).fit_transform(feats)
# feats_spc = PCA(n_components=2).fit_transform(StandardScaler().fit_transform(feats))

# fig, axes = plt.subplots(1, 2, figsize=(9, 5))
# axes[0].scatter(feats_upc[:, 0], feats_upc[:, 1], c=labels)
# axes[0].set_title('Unscaled PCA')
# axes[1].scatter(feats_spc[:, 0], feats_spc[:, 1], c=labels)
# axes[1].set_title('Scaled PCA')
# plt.show()

# %%
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.decomposition import SparsePCA, KernelPCA
from sklearn.decomposition import PCA

k=2

lPCA = PCA(n_components=k)
lPCA.fit(X)
X_PCA = lPCA.transform(X)

sPCA = SparsePCA(n_components=k)
sPCA.fit(X)
X_sPCA = sPCA.transform(X)

# lda = LinearDiscriminantAnalysis(n_components=k)
# lda.fit(X,y)
# X_lda = lda.transform(X)

fig, axes = plt.subplots(1,2,figsize=(16,4))

axes[0].scatter(X_PCA[:,0], X_PCA[:,1],c=y)
axes[0].set_title('PCA')

axes[1].scatter(X_sPCA[:,0], X_sPCA[:,1],c=y)
axes[1].set_title('Sparse PCA')

# axes[2].scatter(X_lda[:,0], X_lda[:,1],c=y)
# axes[2].set_title('LDA')

plt.show()

# %% TSNE 
from sklearn.manifold import TSNE

tsne = TSNE(n_components=2, perplexity=20.0,
            early_exaggeration=12.0,
            learning_rate=200.0,
            n_iter=1000,
            init='random')

X_tsne = tsne.fit_transform(feats)

fig, ax = plt.subplots()
ax.scatter(X_tsne[:,0], X_tsne[:,1], c=y)
ax.set_title('TSNE')


# %% K-means clustering

from sklearn.cluster import KMeans

# Define the number of clusters
num_clusters = 2

# Initialize KMeans model
kmeans = KMeans(n_clusters=num_clusters)

# Fit the KMeans model to the data
# kmeans.fit(data_all_1_array)
kmeans.fit(X_sPCA)

# Get cluster centers and labels
cluster_centers = kmeans.cluster_centers_
labels = kmeans.labels_


# Print cluster centers and labels
print("Cluster Centers:\n", cluster_centers)
print("\nLabels:\n", labels)

import matplotlib.pyplot as plt

# Plotting cluster centers
plt.scatter(cluster_centers[:, 0], cluster_centers[:, 1], c='red', marker='x', label='Cluster Centers')

# Plotting data points with labels
for i in range(num_clusters):
    cluster_points = data_all_1_array[labels == i]
    plt.scatter(cluster_points[:, 0], cluster_points[:, 1], label=f'Cluster {i}')

plt.xlabel('Feature 1')
plt.ylabel('Feature 2')
plt.title('KMeans Clustering Result - sPCA')
plt.legend()
plt.show()

# %% 그나마 제일 그럴듯한 plot (but still don't know the real label - whether it is clustered well or not)

import matplotlib.pyplot as plt

# Create a list of feature pairs
feature_pairs = [(i, j) for i in range(5) for j in range(i+1, 5)]

# Calculate the number of rows and columns for subplots
num_rows = 2
num_cols = 5

# Create subplots
fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, 10))

# Iterate over each pair of features
for idx, pair in enumerate(feature_pairs):
    feature1, feature2 = pair
    row = idx // num_cols
    col = idx % num_cols
    
    # Iterate over each cluster
    for i in range(num_clusters):
        cluster_data = data_all_1_array[labels == i]
        axes[row, col].scatter(cluster_data[:, feature1], cluster_data[:, feature2], label=f'Cluster {i}', alpha=0.5)
        axes[row, col].set_xlabel(f'Feature {feature1+1}')
        axes[row, col].set_ylabel(f'Feature {feature2+1}')
        axes[row, col].set_title(f'Feature {feature1+1} vs Feature {feature2+1}')
        axes[row, col].legend()

# Adjust layout
plt.tight_layout()
plt.show()

# %% GMM - 이 코드 이상함 - 다시
from sklearn.mixture import GaussianMixture

def show_data(X, label=None, centroids=None, cmap='tab10', colorbar=False):
    plt.figure(figsize=(5+int(colorbar),5))
    plt.scatter(*X.T, c=label, cmap=cmap)
    if colorbar:
        plt.colorbar()
    if centroids is not None:
        plt.scatter(*centroids.T, marker='*', c='green', s=100)
    plt.show()
    
gmm = GaussianMixture(n_components=5)
gmm.fit(X)

from sklearn.cluster import MeanShift
clustering = MeanShift(bandwidth=2).fit(X)
clustering.labels_

show_data(X, clustering.labels_, centroids=gmm.means_)


# %% K-means clustering - example

import numpy as np
import matplotlib.pyplot as plt

from sklearn import datasets

np.random.seed(0)

def show_data(X, label=None, centroids=None, cmap='tab10', colorbar=False):
    plt.figure(figsize=(5+int(colorbar),5))
    plt.scatter(*X.T, c=label, cmap=cmap)
    if colorbar:
        plt.colorbar()
    if centroids is not None:
        plt.scatter(*centroids.T, marker='*', c='green', s=100)
    plt.show()
    
X, y = datasets.make_blobs(n_samples=500, n_features=2, centers=4, cluster_std=[2.0, 0.9, 1.2, 1.1], center_box=(-10.0, 10.0), random_state=1)
show_data(X)

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

k = 2
km = KMeans(n_clusters=k)
cluster_labels = km.fit_predict(X)

show_data(X, cluster_labels, centroids=km.cluster_centers_)



# %% (Not Included)
# %% PCA - gpt
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

# Define the number of clusters
num_clusters = 2

# Initialize KMeans model
kmeans = KMeans(n_clusters=num_clusters)

# Fit the KMeans model to the original data
kmeans.fit(data_all_1_array)

# Get cluster centers
cluster_centers = kmeans.cluster_centers_

# Perform PCA to reduce dimensionality from 5 to 2
pca = PCA(n_components=5)  # Corrected to match the desired reduced dimensionality
data_all_1_array_pca = pca.fit_transform(data_all_1_array)

# Fit KMeans on the PCA-transformed data
kmeans.fit(data_all_1_array_pca)
labels = kmeans.labels_
cluster_centers_pca = pca.transform(kmeans.cluster_centers_)

# Plotting cluster centers
plt.scatter(cluster_centers_pca[:, 0], cluster_centers_pca[:, 1], c='red', marker='x', label='Cluster Centers')

# Plotting data points with labels
for i in range(num_clusters):
    cluster_points = data_all_1_array_pca[labels == i]
    plt.scatter(cluster_points[:, 0], cluster_points[:, 1], label=f'Cluster {i}')

# Plotting real labels
for i in range(num_clusters):  # Should iterate over unique labels, not num_clusters
    true_cluster_points = data_all_1_array_pca[data_y_all_1_array == i]
    plt.scatter(true_cluster_points[:, 0], true_cluster_points[:, 1], label=f'Real Label {i}', marker='o', alpha=0.5)

plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('KMeans Clustering Result with PCA')
plt.legend()
plt.show()
