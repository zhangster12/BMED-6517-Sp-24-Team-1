from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import numpy as np

exec(open("data_processing.py").read())

all_features_MI = []
for subj_data in dataframes_MI:
    for stim_index in [0, 2]:  # Stimuli 0 and 2 only
        # Access time and feature values using access_df function
        time, feat_values = access_df(dataframes_MI, dataframes_MI.index(subj_data), 
                                       subj_data.index(stim_index), dataframes_MI[0].index(stim_index))
        # Stack time and feature values into a single array
        features = np.stack((time, feat_values), axis=-1)
        all_features_MI.append(features)

all_features_Ht = []
for subj_data in dataframes_Ht:
    for stim_index in [0, 2]:  # Stimuli 0 and 2 only
        # Access time and feature values using access_df function
        time, feat_values = access_df(dataframes_Ht, dataframes_Ht.index(subj_data), 
                                       subj_data.index(stim_index), dataframes_Ht[0].index(stim_index))
        # Stack time and feature values into a single array
        features = np.stack((time, feat_values), axis=-1)
        all_features_Ht.append(features)

# Stack all the features into a single matrix for MI and HT
feature_matrix_MI = np.vstack(all_features_MI)
feature_matrix_Ht = np.vstack(all_features_Ht)

# Normalize the feature matrices
scaler_MI = StandardScaler()
normalized_features_MI = scaler_MI.fit_transform(feature_matrix_MI)

scaler_Ht = StandardScaler()
normalized_features_Ht = scaler_Ht.fit_transform(feature_matrix_Ht)

# Apply k-means clustering to MI and HT separately
k = 2  # Number of clusters
kmeans_MI = KMeans(n_clusters=k, random_state=42)
cluster_labels_MI = kmeans_MI.fit_predict(normalized_features_MI)

kmeans_Ht = KMeans(n_clusters=k, random_state=42)
cluster_labels_Ht = kmeans_Ht.fit_predict(normalized_features_Ht)

# Assuming true labels for MI and HT subjects
true_labels_MI = np.zeros(len(dataframes_MI))
true_labels_Ht = np.ones(len(dataframes_Ht))

# Calculate ARI and NMI for MI
ari_MI = adjusted_rand_score(true_labels_MI, cluster_labels_MI)
nmi_MI = normalized_mutual_info_score(true_labels_MI, cluster_labels_MI)

# Calculate ARI and NMI for HT
ari_Ht = adjusted_rand_score(true_labels_Ht, cluster_labels_Ht)
nmi_Ht = normalized_mutual_info_score(true_labels_Ht, cluster_labels_Ht)

print("Clustering Accuracy for MI:")
print("Adjusted Rand Index (ARI):", ari_MI)
print("Normalized Mutual Information (NMI):", nmi_MI)

print("\nClustering Accuracy for HT:")
print("Adjusted Rand Index (ARI):", ari_Ht)
print("Normalized Mutual Information (NMI):", nmi_Ht)