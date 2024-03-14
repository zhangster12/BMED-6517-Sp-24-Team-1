# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 13:59:09 2021

@author: jaber
"""

from scipy.io import loadmat
import scipy
from scipy import signal
from scipy.fftpack import fft, ifft
from scipy.signal import hilbert, chirp, find_peaks
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from .pipelineConfig import *

import numpy as np


def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
    return b, a


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_filter(btype, data, cutoff, fs, order=5):

    if btype == 'high':
        b, a = butter_highpass(cutoff, fs, order=order)
        y = signal.filtfilt(b, a, data)
        return y
    elif btype == 'low':
        b, a = butter_lowpass(cutoff, fs, order=order)
        y = signal.filtfilt(b, a, data)
        return y


def filt_lpf(sig_raw, Fs, cutoff=1):

    sig_filt = butter_filter('low', sig_raw, cutoff, Fs)

    return sig_filt


def filt_hpf(sig_raw, Fs, cutoff=1):

    sig_filt = butter_filter('high', sig_raw, cutoff, Fs)

    return sig_filt


def filt_bpf(sig_raw, Fs, lowcutoff=1, highcutoff=40):

    sig_filt = butter_filter('high', sig_raw, lowcutoff, Fs)
    sig_filt = butter_filter('low', sig_filt, highcutoff, Fs)

    return sig_filt


def get_kaiser_bpf(fstop1, fstop2, beta, order, fs):
    
    b = scipy.signal.firwin(order, np.array([fstop1, fstop2]), window=('kaiser',beta), pass_zero='bandpass', scale=True, fs=Fs)
    a = 1.0
    return b,a


def get_filtered_SCG(sig_raw, Fs, lowcutoff=1, highcutoff=40):
    return filt_bpf(sig_raw, Fs, lowcutoff, highcutoff)


def get_filtered_ECG(sig_raw, Fs, lowcutoff=10, highcutoff=40):
    return filt_bpf(sig_raw, Fs, lowcutoff, highcutoff)


def get_filtered_PPG(sig_raw, Fs, lowcutoff=1, highcutoff=8):
    return filt_bpf(sig_raw, Fs, lowcutoff, highcutoff)


def get_padded_filt(sig_raw, filter_padded=5, lowcutoff=0.1, highcutoff=40, Fs=FS_RESAMPLE):

    # reversely pad 5 sec to beginning and end of the signal
    offset = int(filter_padded*Fs)
    sig_raw = np.r_[sig_raw[0:offset][::-1], sig_raw, sig_raw[-offset:][::-1]]

    if lowcutoff == 0:
        sig_filt = filt_lpf(sig_raw, Fs, cutoff=highcutoff)
    elif highcutoff == 0:
        sig_filt = filt_hpf(sig_raw, Fs, cutoff=lowcutoff)
    else:
        sig_filt = filt_bpf(sig_raw, Fs, lowcutoff=lowcutoff, highcutoff=highcutoff)

    sig_filt = sig_filt[offset:-offset]

    return sig_filt


def get_padded_filt_DSwrapper(ts, sig, filter_padded, lowcutoff, highcutoff, Fs):
    downsample_factor = 5
    # ts = patch_dict['time']
    Fs_DS = Fs // downsample_factor  # 100 Hz
    ts_DS = ts[::downsample_factor]

    sig_DS = np.interp(ts_DS, ts, sig)
    sig_filt = get_padded_filt(sig_DS, filter_padded=filter_padded,
                               lowcutoff=lowcutoff, highcutoff=highcutoff, Fs=Fs_DS)
    sig_filt = np.interp(ts, ts_DS, sig_filt)

    return sig_filt


def arg_lowvar(data, lowcutoff=20, highcutoff=5):
    thre_max = np.median(data.var(axis=0))*highcutoff
    thre_min = np.median(data.var(axis=0))/lowcutoff
    indince_clean = np.where((data.var(axis=0) < thre_max) &
                             (data.var(axis=0) > thre_min))[0][:-1]
    return indince_clean


def arg_CI(data, confidence_interv=80):
    value_up = np.percentile(data, int((100-confidence_interv)/2+confidence_interv))
    value_lw = np.percentile(data, int((100-confidence_interv)/2))
    indices_CI = np.where((data <= value_up) & (data >= value_lw))[0]
    return indices_CI


def arg_sqi(data_sqi, confidence_interv=80):
    value_lw = np.percentile(data_sqi, int(100-confidence_interv))
    indices_good = np.where(data_sqi >= value_lw)[0]

    return indices_good


def get_corr(y1, y2):
    '''
    This function computes the normalized cross-correlation of signals y1 and y2. 
    y1 dim: (N,)
    y2 dim: (N,)
    corr dim: scalar
    '''
    n = len(y1)
    corr = signal.correlate(y2, y1, mode='same') / np.sqrt(signal.correlate(y1,
                                                                            y1, mode='same')[int(n/2)] * signal.correlate(y2, y2, mode='same')[int(n/2)])
    return corr


def sig_smoothing(data, lowcutoff, fs=2000):
    sig_filt = np.zeros(data.shape)
    for i in range(data.shape[1]):
        sig_filt[:, i] = butter_filter('high', data[:, i], lowcutoff, fs)
    return sig_filt


def get_smooth(data, N=5):
    # padding prevent edge effect
    data_padded = np.pad(data, (N//2, N-1-N//2), mode='edge')
    data_smooth = np.convolve(data_padded, np.ones((N,))/N, mode='valid')
    return data_smooth


def get_causal_moving_average(data, window=5):
    len_data = len(data)
    avg_data = data
    for i in range(len_data):
        low_ind = max(i-window, 0)
        upper_ind = min(i+window, len_data)
        avg_data[i] = np.mean(data[low_ind:upper_ind])

    return avg_data


def get_smooth_time_based(data, t_data, t_win=5):
    '''
    Smooth data with using previous 't_win' seconds as window size.

    Parameters
    --------------------------
    t_data (Nx1):
        Time vector corresponding to data vector. MUST BE IN SECONDS.
    t_win (int):
        Size of window in seconds. 
    '''
    smoothed_data = data
    for i,t in enumerate(t_data):
        # Get indices from current time point minus the time window
        # up to the current time
        idx = np.where(np.logical_and(t_data <= t, t_data >= (t-t_win)))[0]

        # Take average of those data points
        if len(idx) != 0:
            smoothed_data[i] = np.mean(data[idx])

    return smoothed_data        



def convert_samples_to_sec(samples, Fs=FS_RESAMPLE):

    seconds = samples/Fs
    return seconds


def convert_samples_to_ms(samples, Fs=FS_RESAMPLE):

    msec = 1000*samples/Fs
    return msec


def convert_samples_to_min(samples, Fs=FS_RESAMPLE):

    mins = np.round(samples/(Fs*60),3)
    return mins


def convert_samples_to_hrs(samples, Fs=FS_RESAMPLE):

    hrs = samples/(Fs*60**2)
    return hrs


def convert_sec_to_samples(seconds, Fs=FS_RESAMPLE):

    samples = seconds*Fs
    return np.asarray(samples, dtype=int)


def medfilt(x, k):
    """Apply a length-k median filter to a 1D array x.
    Boundaries are extended by repeating endpoints.
    """
    assert k % 2 == 1, "Median filter length must be odd."
    assert x.ndim == 1, "Input must be one-dimensional."
    k2 = (k - 1) // 2
    y = np.zeros((len(x), k), dtype=x.dtype)
    y[:, k2] = x
    for i in range(k2):
        j = k2 - i
        y[j:, i] = x[:-j]
        y[:j, i] = x[0]
        y[:-j, -(i+1)] = x[j:]
        y[-j:, -(i+1)] = x[-1]
    return np.median(y, axis=1)


def arg_smooth(data, N=5, thre_scale=5):
    # remove outliers that are too far away from moving average
    # 5 point moving average
    data_smooth = get_smooth(data, N=N)

    # find outlier threshold
    smooth_diff = np.abs(data-data_smooth)
    indices_CI = arg_CI(smooth_diff, confidence_interv=80)
    thre = smooth_diff[indices_CI].mean()*thre_scale

    # exclude outlier indices
    indices_good = np.where(smooth_diff < thre)[0]
    return indices_good


def sliding_window(data, window_length, overlap:float=0):
    '''
    _summary_

    Parameters
    ----------
    data : _type_
        _description_
    window_length : _type_
        _description_
    overlap : float, optional
        0 to 1, by default 0 for no overlap

    Returns
    -------
    _type_
        _description_
    '''
    # Sliding window function adapted from Matlab Cardio repo
    windowed = []
    FLAG = True
    idx = 0
    counter = 0
    indices = []
    # The amount to go forward for each window, ensure it's at least 1 or else
    # you will get stuck here infinitely (don't ask me how I know that...)
    jump = max(int(np.floor((1 - overlap)*window_length)), 1)

    # While end of signal hasn't been reached
    while FLAG:
        # Determine whether the end of the signal has been reached
        if idx + window_length >= len(data):

            # If so, compute the buffer size and append to the final window
            bufferLen = window_length - (len(data) - idx)
            # Pad with nans if needed
            buffer = np.empty((bufferLen))
            buffer[:] = np.nan
            finalWindow = data[idx:]
            if len(buffer) != 0:
                finalWindow = np.append(finalWindow, buffer)
            # else:

            windowed.append(finalWindow)

            # Append window indices to index placeholder
            indices.append([idx, len(data)])

            # Throw the flag
            FLAG = False

        else:

            # Get the end index
            endIdx = idx + window_length

            # Append the window to the return value
            windowed.append(data[idx:endIdx])

            # Append window indices to index placeholder
            indices.append([idx, endIdx])

        # Increment the index to obtain the desired overlap
        idx = idx + jump

    return windowed, indices


def sliding_windows(data, window_length, overlap):
    """
    Create sliding windows of time-series data.

    Parameters:
        data (list or numpy array): Time-series data.
        window_length (int): Length of each window.
        overlap (int): Number of data points to overlap between consecutive windows.

    Returns:
        list of numpy arrays: List of sliding windows.
        list of tuples: List of tuples containing start and end indices of each window.
    """
    windows = []
    window_indices = []
    start = 0
    while start < len(data):
        end = min(start + window_length, len(data))
        if (end - start) < window_length:
            break
        windows.append(data[start:end])
        window_indices.append((start, end))
        start += window_length - overlap

    return windows, window_indices



def sliding_time_window(data_time, time_window, time_shift=1):
    '''
    Generates window indices based on time, not just element-based sliding 
    window.

    Parameters
    ----------
    data_time : TYPE
        DESCRIPTION.
    time_window : TYPE
        Length of window in time (s).
    time_shift : TYPE, optional
        Number of seconds to increment. The default is 1.

    Returns
    -------
    windowed_data : TYPE
        DESCRIPTION.
    windowed_indices : TYPE
        DESCRIPTION.

    '''

    windowed_data = []
    windowed_indices = []

    if len(data_time) < 2:
        return [], []
    # Get first and last time points
    init_time = data_time[0]
    fin_time = data_time[-1]

    # Available time
    time_avail = fin_time-init_time

    # Compute amount to shift the window and the resulting number of windows
    num_wins = int(np.ceil((time_avail-time_window)/time_shift) + 1)

    # Iterate through number of windows finding the indices and data vector
    for i in range(num_wins-1):
        # Get current bounds of window
        start_win = init_time + (i*time_shift)
        end_win = start_win + time_window

        curr_win_times = data_time[np.nonzero(
            (data_time >= start_win) & (data_time < end_win))]

        # If no available time points exist in the current window
        if len(curr_win_times) == 0:
            continue

        # Get indices corresponding to the window times
        # start_idx = np.argmin(np.abs(data_time - start_win))
        # end_idx = np.argmin(np.abs(data_time - end_win))
        start_idx = np.where(data_time == curr_win_times[0])[0][0]
        end_idx = np.where(data_time == curr_win_times[-1])[0][0]

        # Store these indices
        windowed_indices.append([start_idx, end_idx+1])
        # Store the time vector included in this window
        # windowed_data.append(data_time[start_idx:end_idx])
        windowed_data.append(curr_win_times)

    # Get final window, force last point to be included
    curr_win_times = data_time[(data_time >= (fin_time-time_window)) & (data_time <= fin_time)]
    start_idx = np.where(data_time == curr_win_times[0])[0][0]
    end_idx = np.where(data_time == curr_win_times[-1])[0][0]

    windowed_indices.append([start_idx, end_idx+1])
    windowed_data.append(curr_win_times)

    return windowed_data, windowed_indices


def get_time_ensemble_beats(beats, i_beats, t_win_len=10):

    # Get bounds of windows (in seconds)
    bound = (t_win_len-1)/2
    # Convert beat indices to time (seconds)
    time_beats = convert_samples_to_sec(i_beats)

    ensemble_beats = []

    for i, t_beat in enumerate(time_beats):
        # Update search window
        lower_bnd = t_beat - bound
        upper_bnd = t_beat + bound
        # Get indices of beats within this window
        i_beats_win = np.nonzero((time_beats >= lower_bnd) & (time_beats <= upper_bnd))[0]
        # If there's too few beats in this window, don't average them
        if len(i_beats_win) < 3:
            ensemble_beats.append(beats[:, i])
        # Otherwise average the beats in this window together
        else:
            ensemble_beats.append(beats[:, i_beats_win].mean(axis=1))

    # Stack the beats back into a numpy array
    ensemble_beats = np.stack(ensemble_beats, axis=1)

    return ensemble_beats


def get_ensemble_beats(sig_beats, N_enBeats=5):
    # sig_beats = (N_samples, N_beats)
    # N_enBeats = # of beats to average together

    # this function looks at the past N beats of a signal (including itself) and
    # uses Woody's method to produce an ensemble beat
    ensemble_beats = []

    # Determine upper and lower window bound
    # bound = np.floor((N_enBeats-1)/2)

    for i_beat in range(sig_beats.shape[1]):
        # If looking at one of the first N beats
        if i_beat-N_enBeats < 0:
            if i_beat == 0:
                ensemble_beats.append(sig_beats[:, i_beat])
            else:
                ensemble_beats.append(sig_beats[:, :i_beat].mean(axis=1))
        else:
            tmp_ens_beats = create_template_woodys(sig_beats[:, i_beat-N_enBeats:i_beat])
            ensemble_beats.append(tmp_ens_beats)
            # ensemble_beats.append(sig_beats[:, i_beat-N_enBeats:i_beat].mean(axis=1))

    ensemble_beats = np.stack(ensemble_beats, axis=1)

    return ensemble_beats


def rolling_ensemble(beats, win_len=10, use_woodys=False):
    # beats (M,N): N signals of length M
    # win_len (scalar): length of rolling window

    # Get the number of signals
    num_beats = np.shape(beats)[1]

    # Determine upper and lower window bound
    bound = int(np.ceil(win_len/2))

    # Placeholder for return value
    ensembled_beats = beats

    # Iterate through each signal slice
    for i in range(num_beats):

        # Determine the upper and lower bounds
        low_bound = max(0, i-bound)
        up_bound = min(num_beats, i+bound)

        if use_woodys:
            # Use Woody's method to average the signals
            ensembled_beats[:, i] = create_template_woodys(beats[:, low_bound:up_bound])
        else:
            ensembled_beats[:, i] = np.mean(beats[:, low_bound:up_bound], axis=1)

    return np.squeeze(ensembled_beats)

def create_template_woodys(beats):
    '''
    Create template via Woody's method of ensemble averaging

    Arguments:
        beats
    '''
    # Normalize signal (align_signals function handles the normalization)
    # sig = (signals - np.mean(signals))/np.std(signals)

    # Initialize template as first beat
    template = beats[:, 0]

    beat_aligned = np.zeros(np.shape(beats))
    beat_aligned[:, 0] = template

    for i in range(np.shape(beats)[1]-1):

        # Align beat to template
        beat_aligned[:, i+1], _ = align_signals(beats[:, i+1], template)

        # Update template by averaging the aligned beats
        template = np.mean(beat_aligned[:, 0:i+1], axis=1)

    return template


def align_signals(sig1, sig2, fs=FS_RESAMPLE, max_lag=0.1, plotFlag=False):
    '''
    Aligns two signals using max cross correlation to compute the lag. Note, only 
    the first signal is modified for alignment.

    Parameters:
        sig1 (Nx1): signal to be modified for alignment
        sig2 (Nx1): reference signal to be aligned to
        fs (scalar): the sampling frequency of the signal
        max_lag (scalar): the maximum number of lags allowed (in percent of signal length)
        plotFlag (bool): indicate whether to plot the alignment

    Returns:
        sig1c (Nx1): aligned signal 1
        lag_samples (scalar): number of samples that sig1 lags or leads

    '''

    # Signal 1 is always the modified one
    sig1n = (sig1 - np.mean(sig1))/np.std(sig1)
    sig2n = (sig2 - np.mean(sig2))/np.std(sig2)
    # If max_lag is to be used
    if max_lag != 0:
        max_lag = int(max_lag*len(sig1n))
        # Note the weird indexing happens because scipy has no max lag parameter
        corr = signal.correlate(sig1n, sig2n, mode='full')[
            len(sig1n)-max_lag-1:len(sig1n)+max_lag]
        lags = signal.correlation_lags(sig1n.size, sig2n.size, mode='full')[
            len(sig1n)-max_lag-1:len(sig1n)+max_lag]
    else:
        corr = signal.correlate(sig1n, sig2n, mode='full')
        lags = signal.correlation_lags(sig1n.size, sig2n.size, mode='full')

    lag_samples = lags[np.argmax(corr)]
    lag_s = (lag_samples / fs)

    # sig1 lags sig2
    if lag_samples > 0:
        # Truncate beginning and pad at end
        sig1c = np.concatenate(
            [sig1[np.abs(lag_samples):], np.zeros((np.abs(lag_samples)))])
        sig1c_plot = np.concatenate(
            [sig1n[np.abs(lag_samples):], np.zeros((np.abs(lag_samples)))])

    # sig2 lags sig1
    elif lag_samples < 0:
        # Truncate end and pad at beginning
        sig1c = np.concatenate([np.zeros((np.abs(lag_samples))),
                               sig1[:len(sig1)-np.abs(lag_samples)]])
        sig1c_plot = np.concatenate(
            [np.zeros((np.abs(lag_samples))), sig1n[:len(sig1n)-np.abs(lag_samples)]])
    else:
        sig1c = sig1
        sig1c_plot = sig1n

    if plotFlag:
        fig = plt.figure(figsize=(16, 10), dpi=80)
        ax = fig.add_subplot(2, 1, 1)
        ax.plot(sig1n, color='red', alpha=0.5)
        ax.plot(sig2n, color='blue', alpha=0.5)
        ax.legend(['Signal 1', 'Signal 2'])
        ax.set_ylim(-15, 15)
        ax.title.set_text("Original")

        ax = fig.add_subplot(2, 1, 2)
        ax.plot(sig1c_plot, color='red', alpha=0.5)
        ax.plot(sig2n, color='blue', alpha=0.5)
        ax.legend(['Signal 1', 'Signal 2'])
        ax.set_ylim(-15, 15)
        ax.title.set_text("Aligned, Number of seconds lag:" + str(lag_s)+' '+'seconds')

    return sig1c, lag_samples


def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'same') / w


def my_ceil(arr, decimal=0):
    return np.ceil(arr*(10**-decimal))/(10**-decimal)


def my_floor(arr, decimal=0):
    return np.floor(arr*(10**-decimal))/(10**-decimal)


def dict_interpolation(raw_dict, fs_resample=FS_RESAMPLE):

    if 'time' in raw_dict.keys():
        t_start = raw_dict['time'][0]
        t_end = raw_dict['time'][-1]
        ecg_time = raw_dict['time']
        ppg_time = raw_dict['time']
        accel_time = raw_dict['time']
        temp_time = raw_dict['time']
    else:
        ecg_time = raw_dict['ecg_time']
        ppg_time = raw_dict['ppg_time']
        accel_time = raw_dict['accel_time']
        
        if 'env_time' in raw_dict.keys():
            temp_time = raw_dict['env_time']
    
        t_start = np.max([ecg_time[0], ppg_time[0], accel_time[0]])
        t_end = np.min([ecg_time[-1], ppg_time[-1], accel_time[-1]])

    time_interp = np.arange(my_ceil(t_start, decimal=-3)*fs_resample,
                            my_floor(t_end, decimal=-3)*fs_resample+1)/fs_resample

    patch_dict = {}

    # ECG
    patch_dict['ecg'] = np.interp(time_interp, ecg_time, raw_dict['ecg'])

    # ACC
    patch_dict['scg_LT'] = np.interp(
        time_interp, accel_time, raw_dict['accel_x'])
    patch_dict['scg_HF'] = np.interp(
        time_interp, accel_time, raw_dict['accel_y'])
    patch_dict['scg_DV'] = np.interp(
        time_interp, accel_time, raw_dict['accel_z'])

    # TODO: find out which PPG arr this belongs to (in terms of its physical location)
    # PPG array 1 (which one is it? left or right?)
    patch_dict['ppg_g_1'] = np.interp(
        time_interp, ppg_time, raw_dict['ppg_g_1'])
    patch_dict['ppg_r_1'] = np.interp(
        time_interp, ppg_time, raw_dict['ppg_r_1'])
    patch_dict['ppg_ir_1'] = np.interp(
        time_interp, ppg_time, raw_dict['ppg_ir_1'])

    # PPG array 2 (which one is it? left or right?)
    patch_dict['ppg_g_2'] = np.interp(
        time_interp, ppg_time, raw_dict['ppg_g_2'])
    patch_dict['ppg_r_2'] = np.interp(
        time_interp, ppg_time, raw_dict['ppg_r_2'])
    patch_dict['ppg_ir_2'] = np.interp(
        time_interp, ppg_time, raw_dict['ppg_ir_2'])

    # Skin temp sensor
    if 'temp_skin' in raw_dict.keys():
        patch_dict['temp_skin'] = np.interp(
            time_interp, temp_time, raw_dict['temp_skin'])

    time_interp = time_interp-time_interp[0]
    patch_dict['time'] = time_interp

    return patch_dict


def get_filt_dict(patch_dict, filt_ecg, filt_scg, filt_ppg, Fs=FS_RESAMPLE, verbose=False):

    patch_filt_dict = {}

    if verbose:
        print('Filtering the raw patch signals...')
        print(F"ECG: {filt_ecg}")
        print(F"SCG: {filt_scg}")
        print(F"PPG: {filt_ppg}")

    patch_filt_dict['time'] = patch_dict['time']

    patch_filt_dict['ecg'] = get_padded_filt(
        patch_dict['ecg'], filter_padded=5, lowcutoff=filt_ecg[0], highcutoff=filt_ecg[1], Fs=Fs)

    patch_filt_dict['ppg_ir_1'] = -get_padded_filt(
        patch_dict['ppg_ir_1'], filter_padded=5, lowcutoff=filt_ppg[0], highcutoff=filt_ppg[1], Fs=Fs)
    patch_filt_dict['ppg_r_1'] = -get_padded_filt(
        patch_dict['ppg_r_1'], filter_padded=5, lowcutoff=filt_ppg[0], highcutoff=filt_ppg[1], Fs=Fs)
    patch_filt_dict['ppg_g_1'] = -get_padded_filt(
        patch_dict['ppg_g_1'], filter_padded=5, lowcutoff=filt_ppg[0], highcutoff=filt_ppg[1], Fs=Fs)
    patch_filt_dict['ppg_ir_2'] = -get_padded_filt(
        patch_dict['ppg_ir_2'], filter_padded=5, lowcutoff=filt_ppg[0], highcutoff=filt_ppg[1], Fs=Fs)
    patch_filt_dict['ppg_r_2'] = -get_padded_filt(
        patch_dict['ppg_r_2'], filter_padded=5, lowcutoff=filt_ppg[0], highcutoff=filt_ppg[1], Fs=Fs)
    patch_filt_dict['ppg_g_2'] = -get_padded_filt(
        patch_dict['ppg_g_2'], filter_padded=5, lowcutoff=filt_ppg[0], highcutoff=filt_ppg[1], Fs=Fs)

    patch_filt_dict['scg_LT'] = get_padded_filt(
        patch_dict['scg_LT'], filter_padded=5, lowcutoff=filt_scg[0], highcutoff=filt_scg[1], Fs=Fs)
    patch_filt_dict['scg_HF'] = get_padded_filt(
        patch_dict['scg_HF'], filter_padded=5, lowcutoff=filt_scg[0], highcutoff=filt_scg[1], Fs=Fs)
    patch_filt_dict['scg_DV'] = get_padded_filt(
        patch_dict['scg_DV'], filter_padded=5, lowcutoff=filt_scg[0], highcutoff=filt_scg[1], Fs=Fs)
    patch_filt_dict['pcg_DV'] = get_padded_filt(
        patch_dict['scg_DV'], filter_padded=5, lowcutoff=FILT_PCG[0], highcutoff=FILT_PCG[1], Fs=Fs)

    return patch_filt_dict



def interpolate_signal(raw_sig, time, resample_freq=FS_RESAMPLE):

    time_interp = np.arange(my_ceil(
        time[0], decimal=-3)*resample_freq, my_floor(time[-1], decimal=-3)*resample_freq+1)/resample_freq
    interp_sig = np.interp(time_interp, time, raw_sig)

    return interp_sig, time_interp


def normalize_sig(signal):

    return (signal - np.mean(signal))/np.std(signal)
