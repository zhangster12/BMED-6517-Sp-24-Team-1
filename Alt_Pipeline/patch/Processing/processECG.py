'''
IMPORT
'''
from ..Tools.segmentation import get_beats_fixed_length
from ..Tools.preprocessing import sliding_window
from ..Tools.pipelineConfig import *
from ..Tools.signalQuality import arg_goodhr, rm_outliers_mad

import neurokit2 as nk      # Ensure that this is installed---https://github.com/neuropsychology/NeuroKit

# Common python packages to import
import matplotlib.pyplot as plt
import numpy as np

'''
INTERNAL HELPER FUNCTION IMPLEMENTATIONS
'''
def find_nearest_peak(i_peaks, sig, Fs, window):
    # Looks for nearest peak with reference to the indices, helps for getting
    # R-peak out of the QRS complex
    i_calibrate = np.zeros(i_peaks.shape)
    
    window_length = int(Fs*window)
    start_offset = -window_length
    end_offset = window_length

    for i, i_peak in enumerate(i_peaks):
        i_start = i_peak + start_offset
        i_end = i_peak + end_offset

        if (i_start < 0) or (i_end > sig.shape[0]):
            i_calibrate[i] = i_peak
            continue

        i_max = sig[i_start:i_end].argmax()
        i_calibrate[i] = i_start + i_max

    # Make sure we remove duplicates that have arisen from this
    return np.unique(i_calibrate.astype(int))


def find_common_r_peaks_within_tolerance(rpks_ref, rpks2, rpks3, tolerance_ms=50, Fs=FS_RESAMPLE):
    '''
    Takes R-peak indices found from multiple methods, finds the set of indices that match (within a 
    specified tolerance).

    Parameters
    ----------
    rpks_ref : numpy array
        _description_
    rpks2 : numpy array
        _description_
    rpks3 : numpy array
        _description_
    tolerance_ms : int, optional
        _description_, by default 50
    Fs : _type_, optional
        _description_, by default FS_RESAMPLE

    Returns
    -------
    _type_
        _description_
    '''
    tolerance_samples = int((tolerance_ms/1000)*Fs)
    rpks_common = set()
    # First iterate through reference peaks for first pass
    for rpk in rpks_ref:
        if any(abs(rpk - r2) <= tolerance_samples for r2 in rpks2) or \
        any(abs(rpk - r3) <= tolerance_samples for r3 in rpks3):
            rpks_common.add(rpk)

    #### THE FOLLOWING BLOCK WOULD SWITCH TO 2/3 METHODS NEEDING AGREEMENT RATHER THAN 3/3

    # # Now remove the rpks we've already found from consideration
    # rpks2_remaining = set(rpks2)
    # for rpk in rpks2:
    #     if any(abs(rpk - rcommon) <= tolerance_samples for rcommon in rpks_common):
    #         rpks2_remaining.remove(rpk)
    
    # # Now parse through final sets to compare, rpks2 vs. rpks3
    # for rpk in rpks2_remaining:
    #     if any(abs(rpk -rpk3) <= tolerance_samples for rpk3 in rpks3):
    #         rpks_common.add(rpk)

    # Return this in numpy array format, make sure to only return unique values and in ascending order
    return np.unique(np.sort(np.array(list(rpks_common))))



def extract_r_peaks(ecg_raw, force_inverse=False, peak_tol=pk_tolerance_ms, width_QRS=width_QRS, verbose=True, Fs=FS_RESAMPLE):

    # Sanitization, common to all
    ecg_sig = nk.signal.signal_sanitize(ecg_raw)
    # Cleaning works for first method, neurokit 
    ecg_cleaned = nk.ecg.ecg_clean(ecg_sig, sampling_rate=Fs, method='neurokit')
    # METHOD 1 -------------------------------------------------
    # Neurokit method
    # Find peaks on inverted signal
    _, info = nk.ecg.ecg_peaks(-ecg_cleaned, sampling_rate=FS_RESAMPLE, method='neurokit', correct_artifacts=True, show=False)
    rpks_neuro = np.asarray(info['ECG_R_Peaks'])
    # Also get the peaks for the filtered signal (no potential inversion)
    _, info = nk.ecg.ecg_peaks(ecg_cleaned, sampling_rate=FS_RESAMPLE, method='neurokit', correct_artifacts=True, show=False)
    # Do a quick check to see if inversion resulted in more or less peaks being found
    # This is based off of the observation that single-lead ECG tends to be inverted by neurokit, when a non-inverted ECG
    # actually results in better R-peak detection. We want the method with the most candidate beats.
    if len(np.asarray(info['ECG_R_Peaks'])) > len(rpks_neuro):
        rpks_neuro = np.asarray(info['ECG_R_Peaks'])
        isInverse = False
        ecg_inv = ecg_cleaned
    else:
        isInverse = True
        ecg_inv = -ecg_cleaned

    # METHOD 2 -------------------------------------------------
    # Martinez method, no cleaning needed
    # Check if we should feed it the inverse or regular ECG
    if isInverse:
        ecg2,_ = nk.ecg.ecg_invert(ecg_sig, sampling_rate=Fs, force=isInverse, show=False)
    else:
        ecg2 = ecg_sig

    _, info = nk.ecg.ecg_peaks(ecg2, sampling_rate=FS_RESAMPLE, method='martinez2004', correct_artifacts=True, show=False)
    rpks_martinez = np.asarray(info['ECG_R_Peaks'])
    # METHOD 3 -------------------------------------------------
    # Cleaning for kalidas2017
    ecg_cleaned3 = nk.ecg.ecg_clean(ecg_sig, sampling_rate=FS_RESAMPLE, method='kalidas2017')
    # Check if we need to feed method 3 the inverse ECG
    if isInverse:
        ecg3, _ = nk.ecg.ecg_invert(ecg_cleaned3, sampling_rate=FS_RESAMPLE, force=isInverse, show=False)
    else:
        ecg3 = ecg_cleaned3

    # Kalidas method
    _, info = nk.ecg.ecg_peaks(ecg3, sampling_rate=FS_RESAMPLE, method='kalidas2017', correct_artifacts=True, show=False)
    rpks_kalidas = np.asarray(info['ECG_R_Peaks'])

    # Combine the 3 R-peak methods to get a common set, using neuro as reference
    rpks_common = find_common_r_peaks_within_tolerance(rpks_neuro, rpks_martinez, rpks_kalidas, tolerance_ms=peak_tol, Fs=FS_RESAMPLE)

    if verbose:
        plt.figure(figsize=(16,6))
        plt.plot(ecg_cleaned, color='gray', alpha=0.6)
        plt.scatter(rpks_neuro, ecg_cleaned[rpks_neuro], color='green', label='neurokit', alpha=0.6)
        plt.scatter(rpks_martinez, ecg_cleaned[rpks_martinez], color='red', label='martinez', alpha=0.6)
        plt.scatter(rpks_kalidas, ecg_cleaned[rpks_kalidas], color='blue', label='kalidas', alpha=0.6)
        plt.legend()

    # Calibrate peaks to nearest max peak
    if isInverse:
        rpks_fin = find_nearest_peak(rpks_common, ecg_inv, FS_RESAMPLE, window=width_QRS)
    else:
        rpks_fin = find_nearest_peak(rpks_common, ecg_cleaned, FS_RESAMPLE, window=width_QRS)
    
    # Print out details regarding R-peak process
    if verbose:
        if isInverse:
            print("ECG signal was inverted.")
        print(f"Number of peaks found:\n   Neurokit ({len(rpks_neuro)})\n   Martinez ({len(rpks_martinez)})\n   Kalidas ({len(rpks_kalidas)})\n   Common ({len(rpks_fin)})")

    return rpks_fin, ecg_cleaned, ecg_inv, isInverse


def get_nni_from_rpeaks(i_R_peaks, r_peak_vals, rmOutliers=True, Fs=FS_RESAMPLE):
    
    # Computes the candidate nn intervals (in s)
    nni_cands = np.ediff1d(i_R_peaks)/Fs
    i_2nd_cands = i_R_peaks[:-1] + np.asarray(nni_cands*Fs, dtype=int)
    # Get corresponding R-peak values
    r_peak_vals = r_peak_vals[:-1]

    # Remove unphysiological NNi's corresponding to heart rate range in settings
    i_physio_mask = arg_goodhr(60/nni_cands, label_range_dict['heart_rate'][1], label_range_dict['heart_rate'][0])

    # Moving window outlier removal
    if rmOutliers:
        i_global_MAD_mask,_ = rm_outliers_mad(60/nni_cands, num_MAD=4, i_slide_win=None)
        # Perform removal if R peak value (on filtered ECG signal) is an outlier, indicative of noise
        # i_global_amp_MAD_mask,_ = rm_outliers_mad(r_peak_vals, num_MAD=5, i_slide_win=None)
        i_tot_mask = i_global_MAD_mask & i_physio_mask # & i_global_amp_MAD_mask
    else:
        i_tot_mask = i_physio_mask

    # Use the masks to extract the good NNi's and the starting points of these intervals
    nni_good = nni_cands[i_tot_mask]
    i_nni = i_R_peaks[:-1][i_tot_mask]
    # Get the second peak of each interval
    i_2nd_peaks = i_2nd_cands[i_tot_mask]

    return i_nni, nni_good, i_2nd_peaks


'''
EXTERNALLY CALLED FUNCTION IMPLEMENTATIONS
'''
def process_ECG(ecg_raw, time, force_inverse=False, width_QRS=width_QRS, peak_tol=pk_tolerance_ms, beat_len=beat_length, Fs=FS_RESAMPLE, verbose=True, plotFlag=True):

    i_R_peaks, ecg_filt, ecg_inv, isInverse = extract_r_peaks(ecg_raw, force_inverse=force_inverse, peak_tol=peak_tol, width_QRS=width_QRS, verbose=verbose, Fs=Fs)
    # Find NN intervals from the R peaks
    i_NN, nn_intervals,_ = get_nni_from_rpeaks(i_R_peaks, ecg_inv[i_R_peaks], rmOutliers=True, Fs=Fs)

    # Get ECG beats
    ecg_beats, _ = get_beats_fixed_length(ecg_filt, i_NN, nn_intervals, Fs=Fs, beat_len=beat_len, padDC=False)

    if plotFlag:
        _,axs = plt.subplots(2,1, figsize=(20, 10), sharex=True)
        axs[0].plot(np.arange(len(ecg_filt))/(Fs*60), ecg_filt, color='red', alpha=0.5)
        axs[0].scatter(i_R_peaks/(Fs*60), ecg_filt[i_R_peaks], color='blue')
        i_not_nn = np.setdiff1d(i_R_peaks, i_NN)
        [axs[0].axvline(i_nn/(Fs*60), color='teal', lw=0.4, alpha=0.4, ls=':') for i_nn in i_not_nn]
        axs[0].set_title("Final Peaks")

        # plt.figure(figsize=(20,6))
        axs[1].set_title('HR')
        axs[1].plot(i_NN/(60*Fs), 60/nn_intervals, color='gray', label='NN-HR', alpha=0.3)
        axs[1].legend(loc='upper right')

    if verbose:
        print(f"{len(i_NN)} NN intervals detected from {len(i_R_peaks)} peaks.")

    ecg_dict = {
        'ecg_filt': ecg_filt,
        'ecg_inv': ecg_inv,             # Returns ECG signal that may be inverted
        'is_ecg_inv': isInverse,        # Boolean indicating if ECG was inverted in the final processing
        'time': time,
        'ecg_beats': ecg_beats,
        'i_beats': i_NN,         # Common R peaks found by the three methods
        't_beats': time[i_NN],
        'NN_int': nn_intervals,
        'hr': 60/nn_intervals,
        'i_segment': i_NN,
    }
    
    return ecg_dict