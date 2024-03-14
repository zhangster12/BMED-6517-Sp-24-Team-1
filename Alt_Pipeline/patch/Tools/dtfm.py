import numpy as np
# from fastdtw import fastdtw
from dtw import *
# from scipy.spatial.distance import euclidean

def distance(signal, template, maxDist=50):

    # SUMMARY
    # This function implements a modified dynamic time warping (DTW) algorithm
    # optimized for peak-matching in addition to waveform morphology matching.
    # This function first extracts the signal and template features and uses it
    # to construct a warping grid that determines anchor points in the warping
    # path. The distance is then computed between signals along each warping path,
    # the optimal path being the one that minimizes the final distance.
    #
    # Arguments (required):
    # - signal      [Nx1]   Signal vector   (Numpy)
    # - template    [Nx1]   Template vector (Numpy)
    #
    # Arguments (optional):
    # - maxDist     Int     Maximum distance between candidate features (default 50)
    #
    # Outputs (dictionary):
    # - distance    Int     Distance between signal and the template (scaled)
    # - length      Int     Length of warped signals
    # - signal      [Mx1]   Warped signal
    # - template    [Mx1]   Warped template
    # - sigPath     [Mx1]   Signal path
    # - tempPath    [Mx1]   Template path

    # EXTRACT FEATURES

    # For the signal segments, get the peaks and valleys
    peaks, valleys = getPeaks(signal)
    fs = np.array(peaks + valleys)
    # Concatenate and sort the features
    idx, fs = np.argsort(fs), np.sort(fs)
    # Determine lables for combined feature vector (peak = 1; valley = 0)
    ts = np.array([1]*len(peaks) + [0]*len(valleys))
    ts = ts[idx]

    # For the template, get the peaks and valleys
    peaks, valleys = getPeaks(template)
    ft = np.array(peaks + valleys)
    # Concatenate and sort the features
    idx, ft = np.argsort(ft), np.sort(ft)
    # Determine labels for combined feature vector (peak = 1; valley = 0)
    tt = np.array([1]*len(peaks) + [0]*len(valleys))
    tt = tt[idx]

    # Set parameters for algorithm
    M, N = len(fs), len(ft)

    # CONSTRUCT GRID

    # Initialize grid to 0
    A = np.zeros((M, N))

    # Where the labels match, set A to 1
    # Remove outliers (where the features are too far apart, set A to 0)
    for n in range(N):
        for m in range(M):
            if ts[m] == tt[n] and np.abs(fs[m] - ft[n]) < maxDist:
                A[m, n] = 1

    # OBTAIN LOSSES FOR EACH PATH

    # Set initial values
    idx_t, idx_s, path_t, path_s = 0, 0, [], []

    # For each template feature...
    for n in range(N):

        # Get the number of possible paths
        numPaths = np.sum(A[:, n])

        # If there is a possible path...
        if numPaths > 0:

            # Set placeholder for possible distances and warp paths
            dist, wp_t, wp_s, coord = [], [], [], []

            # Get indices
            indices = np.where(A[:, n] == 1)
            indices = indices[0]

            # For each possible path...
            for p in range(len(indices)):

                # Set signal feature
                m = indices[p]

                # Segment signal and template
                tmp, sig = template[idx_t:ft[n]+1], signal[idx_s:fs[m]+1]

                # Warp signals
                # distance, path = fastdtw(tmp, sig, dist=euclidean)
                alignment = dtw(sig, tmp, keep_internals=True)
                distance, path_sig, path_tmp = alignment.distance, alignment.index1.tolist(), alignment.index2.tolist()

                # Correct distance by path length
                distance = distance/len(path_tmp)

                # Save values
                dist.append(distance)
                coord.append(list([m, n]))
                # tempPath_t, tempPath_s = [], []
                # for i in range(len(path)):
                #     tempPath_t.append(path[i][0])
                #     tempPath_s.append(path[i][1])
                # wp_t.append(tempPath_t)
                # wp_s.append(tempPath_s)
                wp_t.append(path_tmp)
                wp_s.append(path_sig)

            # Select path with shortest distance
            minDist = min(dist)
            minIdx = dist.index(minDist)

            # Remove invalid paths from grid
            A[:, 0:n+1], A[0:coord[minIdx][0]+1] = 0, 0

            # Save path
            wp_t, wp_s, coord = wp_t[minIdx], wp_s[minIdx], coord[minIdx]
            for i in range(len(wp_t)):
                wp_t[i] += idx_t
            for i in range(len(wp_s)):
                wp_s[i] += idx_s
            # path_t += wp_t
            path_t += wp_t
            path_s += wp_s

            # Update start indices
            idx_t, idx_s = ft[n] + 1, fs[coord[0]] + 1
    
    # Get the warped signals
    warped_s, warped_t = signal[path_s], template[path_t]

    # Get the euclidian distances of the warped signals
    DTFM_distance = np.linalg.norm(warped_s - warped_t)/len(path_t)

    return {"distance":DTFM_distance, "length":len(path_t), "signal":warped_s, "template":warped_t, "sigPath":path_s, "tempPath":path_t}


def getPeaks(signal):

    # SUMMARY
    # This function returns the indices of signal peaks and valleys.
    #
    # Arguments (required):
    # - signal      [Nx1]   Signal vector
    #
    # Outputs:
    # - peaks           [Kx1]   List of peak indices
    # - valleys         [Lx1]   List of valley indices

    # EXTRACT PEAKS AND VALLEYS

    # Differntiate signal
    signalArray = np.array(signal)
    diffSig = np.diff(signalArray)

    # Find peaks and valleys
    peaks, valleys = [], []
    for i in range(1, len(diffSig) - 1):
        if diffSig[i] > 0 and diffSig[i+1] < 0:
            peaks.append(i)
        elif diffSig[i] < 0 and diffSig[i+1] > 0:
            valleys.append(i)

    # Return values
    return peaks, valleys


def score(beats, template, maxDist=50, ell=25):

    # SUMMARY
    # This function calculates the siganl quality index (SQI) of a signal
    # compared with a template using the DTFM distance.
    #
    # Arguments (required):
    # - beats       [NxM]   Signal beats array   (Numpy)
    # - template    [Nx1]   Template vector (Numpy)
    #
    # Arguments (optional):
    # - maxDist     Int     Maximum distance between candidate features (default 50)
    # - lambda (ell)Int     Scaling factor for signal quality index (SQI)
    #
    # Outputs (dictionary):
    # - distance    [Mx1]     Distance between beat and the template
    dist = np.zeros(np.shape(beats)[1])
    for i in range(np.shape(beats)[1]):
        
        # Calculate DTFM distance
        DTFM = distance(beats[:,i], template, maxDist=maxDist)
        dist[i] = DTFM['distance']

    # Calculate SQI
    SQI = np.exp(-ell*dist)

    return SQI


def score(beats, template, maxDist=50, ell=25):

    # SUMMARY
    # This function calculates the siganl quality index (SQI) of a signal
    # compared with a template using the DTFM distance.
    #
    # Arguments (required):
    # - beats       [NxM]   Signal beats array   (Numpy)
    # - template    [Nx1]   Template vector (Numpy)
    #
    # Arguments (optional):
    # - maxDist     Int     Maximum distance between candidate features (default 50)
    # - lambda (ell)Int     Scaling factor for signal quality index (SQI)
    #
    # Outputs (dictionary):
    # - distance    [Mx1]     Distance between beat and the template
    dist = np.zeros(np.shape(beats)[1])
    for i in range(np.shape(beats)[1]):
        
        # Calculate DTFM distance
        DTFM = distance(beats[:,i], template, maxDist=maxDist)
        dist[i] = DTFM['distance']

    # Calculate SQI
    SQI = np.exp(-ell*dist)

    return SQI