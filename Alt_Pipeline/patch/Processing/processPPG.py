'''
IMPORT
'''

import numpy as np
from ..Tools.pipelineConfig import *
from ..Tools.preprocessing import *
from ..Tools.segmentation import *
from ..Tools.signalQuality import *


'''
INTERNAL HELPER FUNCTION IMPLEMENTATIONS
'''
def simpleDiff1Peak(signal, smooth=10):
    
    # Placeholder for return value
    candidates = []

    # Differentiate signal and filter
    diff1Sig = np.diff(signal)
    diff1Sig = get_smooth(diff1Sig, N=smooth)
    # diff1Sig = scsignal.lfilter(filter_coeffs[0], filter_coeffs[1], diff1Sig)
    diff2Sig = np.diff(diff1Sig)
    diff2Sig = get_smooth(diff2Sig, N=smooth)
    # diff2Sig = scsignal.lfilter(filter_coeffs[0], filter_coeffs[1], diff2Sig)

    # Find zero-crossings of second derivative
    for j in range(len(diff2Sig)-1):
        if (diff2Sig[j] > 0) & (diff2Sig[j+1] < 0):
            candidates.append(j)

    # Limit valid candidates
    if len(candidates) > 0:
        candidates = np.asarray(candidates)
        candidates = candidates[np.where(candidates <= foot_range[-1])]
        candidates = candidates[np.where(candidates >= foot_range[0])]

    # Determine optimal candidate, if possible
    if len(candidates) > 0:
        i_max_cand = np.argmax(diff1Sig[candidates])
        max_cand = candidates[i_max_cand]
    else:
        max_cand = []

    i_max = np.argmax(diff1Sig)

    return i_max, max_cand, candidates
    

def tangents(signal, smooth=10):

    # Placeholder for return value
    candidates = []

    # find signal minimum
    sigMin = min(signal)

    # differentiate the signal
    diff1Sig = np.diff(signal)

    # Limit the signal to region before global max if it is after the minimum idx of the search range
    i_max = np.argmax(signal)
    if (i_max > foot_range[0]):
        signal = signal[:i_max]

    # Find the first derivative peak and candidates
    i_max_fd, _, candidates_fd = simpleDiff1Peak(signal, smooth)

    # Find slope of the line at each candidate point and return the 
    # intersection with the minimum
    candidateSlopes = np.zeros((len(candidates_fd,)))
    for j in range(len(candidates_fd)):

        # Find slope of line at candidate point
        candidateSlopes[j] = diff1Sig[candidates_fd[j]]

        # Find intersection point
        intersection = int((sigMin - (signal[candidates_fd[j]] - candidateSlopes[j]*candidates_fd[j]))/candidateSlopes[j])

        candidates.append(intersection)
    
    # Limit valid candidates
    candidates = np.asarray(candidates)
    candidates = candidates[np.where(candidates <= foot_range[-1])]
    candidates = candidates[np.where(candidates >= foot_range[0])]

    # Find the best guess
    bestSlope = diff1Sig[i_max_fd]
    if bestSlope == 0:
        best_index = 0
    else:
        best_index = int((sigMin - (signal[i_max_fd] - bestSlope*i_max_fd))/bestSlope)

    if (best_index > foot_range[-1]) | (best_index < foot_range[0]):
        best_index = 0

    return best_index, candidates


def find_candidate_intersects(ppg_sig, i_max, i_min, N=10):
    # Get first derivative and smooth
    diff1_sig = np.diff(ppg_sig)
    diff1_sig = get_smooth(diff1_sig, N=N)
    # Get second derivative and smooth
    diff2_sig = np.diff(diff1_sig)
    diff2_sig = get_smooth(diff2_sig, N=N)

    # Find zero crossings by multiplying offset versions of signal and checking
    # if result is negative (i.e., one positive and one negative value)
    zero_crossings = pos2_neg_zero_crossing(diff2_sig)
    # Get only viable locations (between max and min)
    i_slopes = zero_crossings[(zero_crossings > i_min) & (zero_crossings < i_max)]

    i_candidates = np.zeros(np.shape(i_slopes), dtype=int)
    for i in range(len(i_slopes)):

        # Find slope of line at candidate point
        cand_slope = diff1_sig[i_slopes[i]]

        # Find the intersection point
        i_candidates[i] = int(i_slopes[i] - (ppg_sig[i_slopes[i]] - ppg_sig[i_min])
                              / diff1_sig[i_slopes[i]])

    # Get best guess of the foot index
    i_best_slope = np.argmax(diff1_sig)
    i_guess_foot = int(i_best_slope - (ppg_sig[i_best_slope]-ppg_sig[i_min])
                       / diff1_sig[i_best_slope])
    # If negative indices are calculated due to some incorrect slope, replace with i_guess_foot
    i_candidates[i_candidates < 0] = i_guess_foot

    return i_guess_foot, i_candidates


def pos2_neg_zero_crossing(signal):
    # Find all zero-crossing points where signal switches from positive to negative
    pos = signal > 0
    zero_crossings = (pos[:-1] & ~pos[1:]).nonzero()[0]
    return zero_crossings


def closest_element(value, array):
    '''
    Returns the closest element from array to the specified value.
    '''
    if not np.isscalar(array):
        i_closest = np.abs(array-value).argmin()
        closest_val = array[i_closest]
    else:
        closest_val = array
    return closest_val


'''
EXTERNALLY CALLED FUNCTION IMPLEMENTATIONS
'''
# Based on Matlab Cardio function pulseTime.m
def extract_PPG_fiducial(ppg_dict, memory=30, smooth=10, max_dist_foot=max_dist_foot):


    reset_number = 10
    # Only use beats that passed the SQI
    ppg_beats = ppg_dict['ppg_beats']
    t_beats = ppg_dict['t_beats']
    beats_used = 'ref_beats'
    # total number of beats
    tot_num_beats = np.shape(ppg_beats)[1]
    pulse_locs = np.zeros((tot_num_beats,))

    # init_template = create_template_woodys(ppg_beats[:,:memory])
    init_template = np.mean(ppg_beats[:,:memory], axis=1)
    init_pulse_loc, _ = tangents(init_template, smooth)

    buffer = []
    invalid_track = 0
    buffer.append(init_pulse_loc)

    for i in range(tot_num_beats):

        beat = ppg_beats[:,i]

        bestGuess, candidates = tangents(beat, smooth)
        
        # Get rid of weird negative pulse locations
        candidates = candidates[candidates > 0] 
        bestGuess = max(bestGuess, 0)

        if (len(candidates) > 0):
            # Select candidate that is closest to prior pulse in memory
            dist = np.zeros((len(candidates),))
            for k in range(len(dist)):
                dist[k] = np.mean(buffer - candidates[k])
            # Take the candidate with the minimum distance from the previous pulses
            i_min_dist = np.argmin(np.abs(dist))
            pulse_locs[i] = candidates[i_min_dist]
            
            if (min(dist) > max_dist_foot) & ((invalid_track <= reset_number) | (min(dist) > (2*max_dist_foot))):
                pulse_locs[i] = INVALID_IDX_CODE
                invalid_track +=1
        # no candidates found
        else:
            # if bestGuess exists (not set to zero as placeholder) and is valid (not negative value)
            if bestGuess > 0:
                pulse_locs[i] = bestGuess
            else:
                # Just assign it to the median of previous pulses
                # pulse_locs[i] = int(np.median(buffer))
                pulse_locs[i] = INVALID_IDX_CODE
                invalid_track += 1

        if pulse_locs[i] != INVALID_IDX_CODE:
            invalid_track = 0
            if len(buffer) == memory:
                buffer.pop(0)
            buffer.append(pulse_locs[i])
        
    # Flag the invalid pulse locations 
    flagged_pulse = np.flatnonzero(pulse_locs != INVALID_IDX_CODE)

    # Return the times of valid pulse locations
    t_pulse = t_beats[flagged_pulse]

    return {
        'pulse_locs': pulse_locs,
        'flagged_pulse': flagged_pulse,
        't_pulse': t_pulse,
        'beats_used': beats_used
    }


def select_ppg_array(ppg_1, ppg_2, verbose=True):
    '''
    _summary_

    Parameters
    ----------
    ppg_1 : _type_
        _description_
    ppg_2 : _type_
        _description_
    verbose : bool, optional
        _description_, by default True

    Returns
    -------
    _type_
        _description_
    '''
    # Coefficient of variation-based PPG array selection
    arr1_cv = abs(np.std(ppg_1)/np.mean(ppg_1))
    arr2_cv = abs(np.std(ppg_2)/np.mean(ppg_2))
    
    if arr1_cv <= arr2_cv:
        ppg_arr = '_1'
    else:
        ppg_arr = '_2'

    if verbose:
        print("Using PPG array "+ppg_arr[-1])

    return ppg_arr, [arr1_cv, arr2_cv]



def process_PPG(ppg, ppg_dc, ecg_dict, ensembleSize=15, useTemplateSQI=False,\
                       beat_len=int(0.6*FS_RESAMPLE), Fs=FS_RESAMPLE, verbose=True):
    '''
    _summary_

    Parameters
    ----------
    ppg : _type_
        _description_
    ppg_dc : _type_
        _description_
    ecg_dict : _type_
        _description_
    ensembleSize : int, optional
        _description_, by default 15
    useTemplateSQI : bool, optional
        _description_, by default False
    Fs : _type_, optional
        _description_, by default FS_RESAMPLE
    verbose : bool, optional
        _description_, by default True

    Returns
    -------
    _type_
        _description_
    '''
    # Get beat-segmentation points
    i_seg = ecg_dict['i_beats']
    if verbose:
        print("Begin processing PPG")

    # Separate ppg signal into beats using R peaks
    ppg_beats, i_beats = get_beats_fixed_length(ppg, i_seg, ecg_dict['NN_int'], Fs=Fs, beat_len=beat_len, padDC=False)
    ppg_dc_beats,_ = get_beats_fixed_length(ppg_dc, i_seg, ecg_dict['NN_int'], Fs=Fs, beat_len=beat_len, padDC=True)

    # Variable to keep track of which PPG beats (with reference to i_beats) are kept
    ref_beats = np.empty(1)
    sqi_mask = np.zeros(len(i_beats), dtype=bool)
    # # %% Remove PPG beats with outlying amplitudes
    # % I found that motion artifacts play a huge role in corrupting the PPG
    # % signal. Most motion artifacts tend to corrupt the PPG signal similarly -
    # % by causing spikes in the signal. Hence, I feel by excluding PPG beats
    # % with outlier amplitudes, I may be able to filter some of those out; and
    # % possibly even help the SQI template creation
    ppg_amp = np.zeros(np.shape(ppg_beats)[1], dtype=float)
    # Remember to normalize by the PPG AC beats by the DC value for each beat
    for i in range(len(ppg_amp)):
        if np.mean(ppg_dc_beats[:,i]) == 0:
            print('Uh oh, someone is 0') 
        ppg_amp[i] = (max(ppg_beats[:, i]) - min(ppg_beats[:, i]))/np.mean(ppg_dc_beats[:,i])

    # # Designate outliers based on global median amplitudes
    i_keep_s1, ppg_amp_s1 = rm_outliers_mad(ppg_amp, num_MAD=5)

    ppg_beats_s1 = ppg_beats[:, i_keep_s1]
    ref_beats = np.argwhere(i_keep_s1)

    if verbose:
        print('Stage 1: Global MAD removed ' + str(np.shape(ppg_beats)
              [1]-np.shape(ppg_beats_s1)[1]) + ' out of ' + str(np.shape(ppg_beats)[1]))
        
    # Designate outliers based on moving window median amplitudes
    win_len = 30
    overlap = 29/30
    _, i_window = sliding_window(i_beats[ref_beats], win_len, overlap)
    i_keep_s2, ppg_amp_s2 = rm_outliers_mad(ppg_amp_s1, num_MAD=5, i_slide_win=i_window)

    ppg_beats_s2 = ppg_beats_s1[:, i_keep_s2]
    ref_beats = ref_beats[i_keep_s2]
    ref_beats_amp = ref_beats


    if verbose:
        print('Stage 2: Moving window MAD removed ' + str(np.shape(ppg_beats_s1)
              [1]-np.shape(ppg_beats_s2)[1]) + ' out of ' + str(np.shape(ppg_beats_s1)[1]))

    # If we want to use template-based SQI
    if useTemplateSQI:
        # Normalize PPG beats for DTW SQI, but note that we will keep the non-normalized PPG beats to 
        # prevent any changes to pulse arrival times
        for i in range(np.shape(ppg_beats_s2)[1]):
            ppg_beats_s2[:,i] = (ppg_beats_s2[:,i] - np.mean(ppg_beats_s2[:,i]))/np.std(ppg_beats_s2[:,i])

        # Compute SQI with 30 beat windows
        win_len = 30
        overlap = 0
        _, i_window = sliding_window(i_beats[ref_beats], win_len, overlap)

        sqi_xcorr = np.zeros(np.shape(ppg_beats_s2)[1])
        sqi_dtw = np.zeros(np.shape(ppg_beats_s2)[1])
        for i in range(len(i_window)):
            curr_beats = ppg_beats_s2[:, i_window[i][0]:i_window[i][-1]]
            curr_template = np.mean(curr_beats, axis=1)
            sqi_xcorr[i_window[i][0]:i_window[i][-1]] = get_sqi(curr_beats, get_max_xcorr, curr_template)
            sqi_dtw[i_window[i][0]:i_window[i][-1]] = get_sqi(curr_beats, get_dtw, curr_template)

        # Keep this as the stage1 SQI
        sqi_stg1 = (sqi_dtw + sqi_xcorr)/2
        sqi_stg1_mask = (sqi_stg1 >= sqi_PPG_threshold)
        ppg_beats_s3 = ppg_beats_s2[:, sqi_stg1_mask]

        if verbose:
            print('Stage 3: Moving window DTW removed ' + str(np.shape(ppg_beats_s2)
                [1]-np.shape(ppg_beats_s3)[1]) + ' out of ' + str(np.shape(ppg_beats_s2)[1]))


        # SQI Stage 2
        # Fixed template
        template = np.mean(ppg_beats_s3, axis=1)

        sqi_dtw2 = get_sqi(ppg_beats_s2, get_dtw, template)
        sqi_stg2 = sqi_dtw2
        sqi_stg2_mask = (sqi_stg2 >= sqi_PPG_threshold)
        ppg_beats_s4 = ppg_beats_s2[:,sqi_stg2_mask]

        for i in range(len(sqi_stg2)):
            # If the stage 2 SQI is super low and we previously marked this beat
            # as good with the binary label 1
            if (sqi_stg2[i] < (lenient_thresh+0.1)) & (sqi_stg1[i] > 0.5):
                sqi_stg2_mask[i] = False
            
            # On the other hand, if stage 2 SQI is super high and we previously
            # marked this beat as bad with the binary label 0
            elif (sqi_stg2[i] >= strict_thresh) & (sqi_stg1[i] < 0.5):
                sqi_stg2_mask[i] = True
    
        if verbose:
            print('Stage 4: Global DTW fusion removed ' + str(np.shape(ppg_beats_s2)[1]-np.shape(ppg_beats_s4)[1]) + ' out of ' + str(np.shape(ppg_beats_s2)[1]))

        # Determine which method removed less beats, keep the resulting beats of that method
        if np.sum(sqi_stg1_mask) <= np.sum(sqi_stg2_mask):
            ppg_beats_good_mask = sqi_stg2_mask
        else:
            ppg_beats_good_mask = sqi_stg1_mask

        # Note, we return to the beats from s1, as s2 beats have been standardized
        ppg_beats_s5 = ppg_beats_s1[:, i_keep_s2][:, ppg_beats_good_mask]
        ref_beats = ref_beats[ppg_beats_good_mask]
        sqi_mask[ref_beats] = 1 

    # No sqi
    else:
        # Use beats after amplitude-outlier rejection
        ppg_beats_s5 = ppg_beats_s2

    if verbose:
        print('Final: Using ' + str(np.shape(ppg_beats_s5)[1]) + ' beats out of ' + str(np.shape(ppg_beats)[1]))

    if ensembleSize > 1:
        ppg_beats_s5 = rolling_ensemble(ppg_beats[:,ref_beats], win_len=ensembleSize)
        ppg_beats = rolling_ensemble(ppg_beats[:,ref_beats_amp], win_len=ensembleSize)
        ppg_dc_beats = rolling_ensemble(ppg_dc_beats[:, ref_beats_amp], win_len=ensembleSize)

    # Go ahead and get PPG amplitude from beats that have passed through amplitude-based filtering
    # and have been ensemble averaged
    ppg_amp = np.zeros(np.shape(ppg_beats)[1])
    ppg_amp_ac = np.zeros(np.shape(ppg_beats)[1])
    ppg_dc = np.zeros(np.shape(ppg_beats)[1])
    # Remember to normalize by the PPG AC beats by the DC value for each beat
    for i in range(len(ppg_amp)):
        ppg_amp[i] = (max(ppg_beats[:, i]) - min(ppg_beats[:, i]))/np.mean(ppg_dc_beats[:,i])
        ppg_amp_ac[i] = (max(ppg_beats[:, i]) - min(ppg_beats[:, i]))
        ppg_dc[i] = np.mean(ppg_dc_beats[:,i])
    # Note, ref_beats is now in reference to ecg_dict['i_segment'], where we want it in reference
    # to ecg_dict['i_beats'] for creating a clean feature matrix/dataframe later
    _, ref_beats,_ = np.intersect1d(ecg_dict['i_beats'], i_beats[ref_beats], return_indices=True)
    _, ref_beats_amp,_ = np.intersect1d(ecg_dict['i_beats'], i_beats[ref_beats_amp], return_indices=True)

    # Get times of final set of beats
    t_beats = ecg_dict['t_beats'][ref_beats]
    # Get times of PPG amplitude values
    t_ppg_amp = ecg_dict['t_beats'][ref_beats_amp]

    return {'ppg_beats': ppg_beats_s5,          # The beats that passed the SQI
            't_beats': t_beats,                 # Times of each PPG beat that passed preprocessing
            'ref_beats': ref_beats,             # Mask on ECG R peaks locations that correspond to the final set of PPG beats
            'ref_beats_amp': ref_beats_amp,     # Mask on ECG R peaks for the set of PPG beats after amplitude-based outlier removal
            't_ppg_amp': t_ppg_amp,             # Times of PPG amplitudes
            'ppg_amp': ppg_amp,                 # PPG amplitudes after outlier removal
            'ppg_amp_ac': ppg_amp_ac,
            'ppg_dc': ppg_dc,
            'og_beats': ppg_beats,              # The PPG beats before SQI (but after amplitude-based outlier removal), still moving window averaged
            'dc_beats': ppg_dc_beats,
            }
