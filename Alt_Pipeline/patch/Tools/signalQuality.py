# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 15:57:21 2022

@author: jaber
"""

from .pipelineConfig import*
from .preprocessing import *
import matplotlib.pyplot as plt
from dtw import *

# # ECG SQI
def get_ecg_R_peak_sqi(ecg_filt, i_PT_peaks, i_SE_peaks, min_f_score=0.7, Fs=FS_RESAMPLE, usePT_Flag=False):
    # 1. 10-second windowing
    window = int(10*Fs)
    windowed_sig, indices = sliding_window(ecg_filt, window, overlap=0)
     
    f_score = np.zeros((np.shape(indices)[0]))
    i_fin_peaks = []
    i_union_peaks = []
    # set maximum agreement between peak locations (50-100 ms is a good #)
    max_pk_dist = int(.1*Fs)
    for i, idx_pair in enumerate(indices):
        # 2. get peaks found by both algorithms for this window
        i_ref_peaks = i_PT_peaks[(i_PT_peaks >= idx_pair[0]) & (i_PT_peaks <= idx_pair[1])]
        i_test_peaks = i_SE_peaks[(i_SE_peaks >= idx_pair[0]) & (i_SE_peaks <= idx_pair[1])]
        
        # Make sure this window even has reference peaks, otherwise F score = 0
        if (len(i_ref_peaks) > 0):
            # 3. match the R peaks: get indices (with reference to i_ref_peaks)
            # of the peaks that agree between the two
            i_matching = match_R_peaks(i_ref_peaks, i_test_peaks, max_pk_dist)
            
            num_ref = len(i_ref_peaks)
            num_test = len(i_test_peaks)
            num_matching = len(i_matching)
            
            # 4. calculate the f-score
            # if there are some similar peaks detected between the two methods
            if num_matching > 0:
                true_pos = num_matching
                false_neg = num_ref - true_pos
                false_pos = num_test - true_pos
                se = true_pos/(true_pos + false_neg)
                ppv = true_pos/(false_pos + true_pos)
                f1 = 2*se*ppv/(se+ppv)
                f_score[i] = f1
                
            # If f_score > predetermined threshold (defined parameter)
            if f1 > min_f_score:
                i_union_peaks.append(i_ref_peaks[i_matching])
                i_fin_peaks.append(i_ref_peaks)
    
    if usePT_Flag:
        # 5. save all PT peaks in the good windows
        i_R_peaks = np.concatenate(i_fin_peaks, axis=0)
        # Final check to make sure we have no duplicates
        i_R_peaks = np.unique(i_R_peaks)
    else:
        # 5. alternatively could just use the union of the two methods for each 
        # window, but this will reduce overall # of peaks, some of which might be
        # viable but are still removed...
        i_R_peaks = np.concatenate(i_union_peaks, axis=0)
        # Final check to make sure we have no duplicates
        i_R_peaks = np.unique(i_R_peaks)
    
    return i_R_peaks, f_score
     

def combine_r_peak_methods(ecg_filt, i_peaks1, i_peaks2, min_f_score=ecg_quality_threshold, Fs=FS_RESAMPLE, use_ref_Flag=False):
    '''
    Compare R peak methods and extract the peaks common to both algorithms

    Parameters
    ----------
    ecg_filt : _type_
        _description_
    i_peaks1 : _type_
        _description_
    i_peaks2 : _type_
        _description_
    min_f_score : float, optional
        _description_, by default 0.7
    Fs : _type_, optional
        _description_, by default FS_RESAMPLE
    use_ref_Flag : bool, optional
        _description_, by default False

    Returns
    -------
    _type_
        _description_
    '''
    # 1. 10-second windowing
    window = int(10*Fs)
    _, indices = sliding_window(ecg_filt, window, overlap=0)
     
    f_score = np.zeros((np.shape(indices)[0]))
    i_fin_ref_peaks = []
    i_union_peaks = []
    # set maximum agreement between peak locations (50-100 ms is a good #)
    max_pk_dist = int(.1*Fs)

    for i, idx_pair in enumerate(indices):
        # 2. get peaks found by both algorithms for this window
        i_ref_peaks = i_peaks1[(i_peaks1 >= idx_pair[0]) & (i_peaks1 <= idx_pair[1])]
        i_test_peaks = i_peaks2[(i_peaks2 >= idx_pair[0]) & (i_peaks2 <= idx_pair[1])]
       
        # Make sure this window even has reference peaks, otherwise F score = 0
        if (len(i_ref_peaks) > 0):
            # 3. match the R peaks: get indices (with reference to i_ref_peaks)
            # of the peaks that agree between the two
            i_matching = match_R_peaks(i_ref_peaks, i_test_peaks, max_pk_dist)
        
            num_ref = len(i_ref_peaks)
            num_test = len(i_test_peaks)
            num_matching = len(i_matching)
            
            # 4. calculate the f-score
            # if there are some similar peaks detected between the two methods
            if num_matching > 0:
                true_pos = num_matching
                false_neg = num_ref - true_pos
                false_pos = num_test - true_pos

                f_score[i] = 2*true_pos/(2*true_pos + false_pos + false_neg)

            # If f_score > predetermined threshold (defined parameter)
            if f_score[i] > min_f_score:
                i_union_peaks.append(i_ref_peaks[i_matching])
            # Store reference peaks 
            i_fin_ref_peaks.append(i_ref_peaks)

    if use_ref_Flag:
        # 5. save all ref peaks
        i_R_peaks = np.concatenate(i_fin_ref_peaks, axis=0)
        # Final check to make sure we have no duplicates
        i_R_peaks = np.unique(i_R_peaks)
    else:
        # 5. alternatively could just use the union of the two methods for each 
        # window, but this will reduce overall # of peaks, some of which might be
        # viable but are still removed...
        # Check if there are any overlapping peaks, may be an error so we'll just use the 
        # reference method in this case
        if len(i_union_peaks) == 0:
            i_R_peaks = i_peaks1
        else:
            i_R_peaks = np.concatenate(i_union_peaks, axis=0)
        # Final check to make sure we have no duplicates
        i_R_peaks = np.unique(i_R_peaks)
    
    return i_R_peaks, f_score



def get_hr_interval_outliers_mask(i_peaks, num_MAD=5, Fs=FS_RESAMPLE):
    
    # Computes the candidate HR interval
    hr_cand = 1/np.diff(i_peaks/Fs)*60
    
    # Remove unphysiological heart rates
    i_hr_physio = arg_goodhr(hr_cand, label_range_dict['heart_rate'][1], label_range_dict['heart_rate'][0])
    
    # Remove outliers +- # MAD away from median hr
    # i_MAD_int_outliers,_ = rm_outliers_mad(hr_cand, num_MAD=num_MAD)
    # Moving window outlier removal
    win_len = 60
    overlap = 59/60
    _,i_window = sliding_window(hr_cand, win_len, overlap)
    i_MAD_moving_outlier,_ = rm_outliers_mad(hr_cand, num_MAD=num_MAD, i_slide_win=i_window)
    
    return i_hr_physio & i_MAD_moving_outlier #i_MAD_int_outliers & 
     


def match_R_peaks(i_ref_peaks, i_test_peaks, max_pk_dist):
    # Returns indices of peak from i_ref that are acceptably close to a peak in
    # i_test
    i_matching = []
    for i_peak in i_test_peaks:
        i_tmp = find_nearest(i_ref_peaks,i_peak)
        # Ensure they are within maximum acceptable distance
        if np.abs(i_ref_peaks[i_tmp] - i_peak) <= max_pk_dist:
            i_matching.append(i_tmp)
          
    if len(i_matching) != 0:  
        # Take only the unique ones to avoid copies
        i_matching = np.unique(i_matching)
    
    return i_matching


def find_nearest(array,value):
    # Find index of closest value in array to a specified value
    idx = (np.abs(array-value)).argmin()
    return idx


def arg_goodhr(data, thre_max=180, thre_min=40):
    # exclude outlier indices
    indices_good = (data <= thre_max) & (data >= thre_min)

    return indices_good


# # PPG SQI
def get_ppg_sqi(ppg_r, ppg_ir, ppg_g, ppg_sqi_thre=ppg_sqi_threshold):
    '''
    Uses the green PPG as a template and computes the cross-correlation of each
    beat from the red and ir PPGs as the sqi. Both the red and ir PPG beats 
    must be above the sqi threshold for it to be used
    '''
    ppg_g_temp = np.mean(ppg_g, axis=1)
    # ppg_g_temp = create_template(ppg_g)

    ppg_r_sqi = get_sqi(ppg_r, get_max_xcorr, ppg_g_temp)
    ppg_ir_sqi = get_sqi(ppg_ir, get_max_xcorr, ppg_g_temp)
    
    mask_sqi = (ppg_r_sqi > ppg_sqi_thre) & (ppg_ir_sqi > ppg_sqi_thre)
        
    return mask_sqi



# # General SQI Methods
def get_norm_corr(y1, y2):
    '''
    This function computes the normalized cross-correlation of signals y1 and y2. 
    y1 dim: (N,)
    y2 dim: (N,)
    corr dim: scalar
    '''
    n = len(y1)
    corr = signal.correlate(y2, y1, mode='same') / np.sqrt(signal.correlate(y1, y1, mode='same')[int(n/2)] * signal.correlate(y2, y2, mode='same')[int(n/2)])
    return corr


def get_max_xcorr(beats, beats_template):
    '''
    This function (distance_method) computes the distance [cross-correlation] between the each beat in beats and beats_template
    beats dim: (N,M)
    beats_template dim: (N,)
    xcorr dim: (M,)
    '''
    corr_matrix = np.zeros(beats.shape)
    for i_col in range(beats.shape[1]):
        corr_matrix[:, i_col] = get_norm_corr(beats_template, beats[:,i_col])

    beats_xcorr = corr_matrix.max(axis=0)
    return beats_xcorr


def get_dtw(beats, beats_template):
    '''
    This function (distance_method) computes the distance [dynamic time warping] between the each beat in beats and beats_template
    beats dim: (N,M)
    beats_template dim: (N,)
    beats_dtw dim: (M,)
    '''
    beats_dtw = np.zeros(beats.shape[1])
    for i in range(beats.shape[1]):
        si = beats[:,i]
        # Old method
        # distance = dtw.distance(si, beats_template)
        
        # Warp signals
        alignment = dtw(si, beats_template, keep_internals=True)
        distance = alignment.distance
        # path_sig, path_tmp = alignment.index1.tolist(), alignment.index2.tolist()
        # # Correct distance by path length
        # distance = distance/len(path_tmp)
        beats_dtw[i] = np.exp(-distance/np.abs(beats_template).sum())
    
    return beats_dtw


def get_dtw_single_beat(beat, beat_template):
    '''
    This function (distance_method) computes the distance [dynamic time warping] between the each beat in beats and beats_template
    beats dim: (N,)
    beats_template dim: (N,)
    beats_dtw dim: (1)
    '''
    # Warp signals
    alignment = dtw(beat, beat_template, keep_internals=True)
    distance = alignment.distance
    beat_dtw = np.exp(-distance/np.abs(beat_template).sum())
    
    return beat_dtw


def get_sqi(beats, distance_method, template=None):
    '''
    Returns the sqi score for every beat.
    Parameters:
            beats (NxM): a beat matrix that store M beats of the signal (each beat signal has a dimension of N)
            distance_method (python function): takes a beat matrix and a template and subsequently compute distance measures between each beat and the template 
                                                output is a N dimensional array
            template (N,):  a beat signal used to compute the distance measure. If None, the ensemble average of the best 10% beats will be used as the template. 
                            Be careful if beats contain mostly noisy signals. An external template would be recommended.
    Returns:
            sqi_distances (M,): 
            
    Usage:
            sqi_xcorr = get_sqi(beats, get_max_xcorr, template)
    '''
    
    best_partition = 10 # use the best 10 % to produce the template
    
    if template is None:
        beats_template0 = beats.mean(axis=1)

        sqi_xcorr0 = distance_method(beats, beats_template0)

        indices_template = np.argsort(-sqi_xcorr0)[:sqi_xcorr0.shape[0]//best_partition] # use the top 10% to comptue a new template
        beats_template = beats[:,indices_template].mean(axis=1)
        # beats_template = create_template(beats[:,indices_template])
        sqi_distances = distance_method(beats, beats_template)

    else:
        beats_template = template
        sqi_distances = distance_method(beats, beats_template)

    return sqi_distances



def rm_outliers_mad(values, num_MAD=5, i_slide_win=None):
        
    # Removes outliers based on global median
    if i_slide_win == None:
        # Remove outliers +- # MAD away from median value
        mad = np.nanmedian(np.abs(values-np.nanmedian(values)))
        i_keep = (values < (np.nanmedian(values) + num_MAD * mad)) & (values > (np.nanmedian(values) - num_MAD * mad))
    # Removes outliers based on sliding window median
    else:
        i_keep = np.zeros(len(values),dtype=bool)
        for start, stop in i_slide_win:
            curr_vals = values[start:stop]
            mad = np.nanmedian(np.abs(curr_vals[~np.isnan(curr_vals)]-np.nanmedian(curr_vals)))            
            # Remove outliers +- # MAD away from median value
            mask = (curr_vals <= (np.nanmedian(curr_vals) + num_MAD * mad)) & \
                (curr_vals >= (np.nanmedian(curr_vals) - num_MAD * mad))
            i_keep[start:stop] = np.squeeze(mask)
      
    return i_keep, values[i_keep]


def rmoutliers(data, data_idx, num_mad=4, win_len=60, overlap=59/60, globalFlag=False):
    '''
    Function takes a feature and performs outlier removal based on global medians
    and sliding window based rejection. 
    '''
    if globalFlag:
        # First do global median outlier flagging
        i_keep_s1, data_s1 = rm_outliers_mad(data, num_MAD=num_mad)
    else:
        i_keep_s1 = np.ones(len(data_idx), dtype=bool)

    data_s1 = data[i_keep_s1]
    data_idx = data_idx[i_keep_s1]

    # Next, find outliers based on moving window median amplitudes
    _, i_window = sliding_window(data_idx, win_len, overlap)
    i_keep_s2, data_s2 = rm_outliers_mad(data_s1, num_MAD=num_mad, i_slide_win=i_window)

    data_s2 = data_s1[i_keep_s2]
    data_idx = data_idx[i_keep_s2]

    return data_s2, data_idx


def flag_outliers(data, data_idx, num_mad=4, win_len=60, overlap=59/60, globalFlag=True):
    '''
    Function takes a feature and finds outlier portions based on global medians
    and sliding window based rejection. 
    '''
    if globalFlag:
        # First do global median outlier flagging
        i_out_s1, data_s1 = rm_outliers_mad(data, num_MAD=num_mad)
    else:
        i_out_s1 = np.zeros_like(data, dtype=bool)

    mask = np.ones_like(data_idx, dtype=bool)
    mask[i_out_s1] = False

    # Remove the data that are not outliers
    data_idx = data_idx[mask]
    data_s1 = data[mask]

    # Next, find outliers based on moving window median amplitudes
    _, i_window = sliding_window(data_idx, win_len, overlap)
    i_out_s2, data_s2 = rm_outliers_mad(data_s1, num_MAD=num_mad, i_slide_win=i_window)

    # Remove the data that are not outliers
    mask = np.ones_like(data_s1, dtype=bool)
    mask[i_out_s2] = False
    data_s2 = data_s1[mask]
    data_idx = data_idx[mask]

    return data_s2, data_idx