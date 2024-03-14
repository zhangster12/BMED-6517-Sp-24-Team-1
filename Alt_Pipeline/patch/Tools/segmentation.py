# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 10:52:59 2022

@author: jaber
"""

import numpy as np
import matplotlib.pyplot as plt
from .pipelineConfig import *
# from PyTools.plotting import *
from .preprocessing import *

def segment_beats(i_beats, sig_filt_dict):
    """
    Segments all signals based off of R-peaks detected from ecg signal.
    Parameters
    ----------
    i_beats : integer array
        Array of indices for the start of beats
    sig_filt_dict : dictionary
        Dictionary of all filtered signals (ecg, ppg, and scg).
    Returns
    -------
    beats_dict : dictionary
        Dictionary of signals segmented into appropriate beats.
    """
    
    

    # 3.1 - Get beat segmentation for ECG
    ecg_beats, _ = beat_segmentation(sig_filt_dict['ecg'], i_beats, start_offset=0, end_offset=end_offset)
    
    # 3.2 - Get beat segmentation for PPG
    ppg_g_1_beats, _ = beat_segmentation(sig_filt_dict['ppg_g_1'], i_beats, start_offset=0, end_offset=end_offset)
    ppg_r_1_beats, _ = beat_segmentation(sig_filt_dict['ppg_r_1'], i_beats, start_offset=0, end_offset=end_offset)
    ppg_ir_1_beats, _ = beat_segmentation(sig_filt_dict['ppg_ir_1'], i_beats, start_offset=0, end_offset=end_offset)
    ppg_g_2_beats, _ = beat_segmentation(sig_filt_dict['ppg_g_2'], i_beats, start_offset=0, end_offset=end_offset)
    ppg_r_2_beats, _ = beat_segmentation(sig_filt_dict['ppg_r_2'], i_beats, start_offset=0, end_offset=end_offset)
    ppg_ir_2_beats, _ = beat_segmentation(sig_filt_dict['ppg_ir_2'], i_beats, start_offset=0, end_offset=end_offset)
    
    
    # 3.3 - Get beat segmentation for SCG (dorso-ventral axis only)
    scg_dv_beats, _ = beat_segmentation(sig_filt_dict['scg_DV'], i_beats, start_offset=0, end_offset=end_offset)
    t_beats, _ = beat_segmentation(sig_filt_dict['time'], i_beats, start_offset=0, end_offset=end_offset)
    # 4 - Aggregate all beats
    b_time = sig_filt_dict['time']
    beats_dict = {
        'ecg': ecg_beats,
        'ppg_g_1': ppg_g_1_beats,
        'ppg_r_1': ppg_r_1_beats,
        'ppg_ir_1': ppg_ir_1_beats,
        'scg_dv': scg_dv_beats,        
        'ppg_g_2': ppg_g_2_beats,
        'ppg_r_2': ppg_r_2_beats,
        'ppg_ir_2': ppg_ir_2_beats,
        'i_R_peaks': i_beats,
        't_beats': t_beats, 
    }
    
    return beats_dict



def beat_segmentation(sig, i_peaks, start_offset=0, end_offset=250):
    # segment a signal using the start index (i_peaks+start_offset) and the end index (i_peaks+end_offset) 
    
    sig_beats = []
    i_beat_peaks = []
    for i, i_peak in enumerate(i_peaks):
        i_start = int(i_peak + start_offset)
        i_end = int(i_peak + end_offset)

        if i_start < 0 or i_end > sig.shape[0]:
            continue

        sig_seg = sig[i_start:i_end]

        sig_beats.append(sig_seg)
        i_beat_peaks.append(i_peak)

    sig_beats = np.vstack(sig_beats).T
    i_beat_peaks = np.vstack(i_beat_peaks).T.squeeze()

    return sig_beats, i_beat_peaks


def segment_into_beats(signal, i_seg, Fs=FS_RESAMPLE, beat_len = int(0.5*FS_RESAMPLE)):
    
    beats = []
    i_beat_start = []
    # t_peaks = i_NN_peaks/Fs
    
    for idx in i_seg:
        # if (t_peaks[i] + nn_intervals[i]) >= 0.5:
        #     beats.append(signal[i_NN_peaks[i]:i_NN_peaks[i]+nn_intervals[i]/Fs]) 
        # else:
        beats.append(signal[idx:idx+beat_len])
        i_beat_start.append(idx)
            
    return np.asarray(beats).T, np.asarray(i_beat_start)



def segment_into_beats_w_padding(signal, i_beats, intervals, Fs=FS_RESAMPLE, padDC=False):
    '''
    Segment a signal using provided indices and lengths of beats, with zero-padding to a fixed length
    '''
    # Compute length of beats from max interval so no important information is removed
    intervals = np.array(intervals*Fs, dtype=int)
    beat_len = max(intervals)
    # First make sure that all intervals are <= beat_len 
    intervals[intervals > beat_len] = beat_len
    beats = np.zeros(shape=(len(i_beats), beat_len), dtype=float)

    for i, beat_int in enumerate(intervals):
        # Get the portion of the signal related to this beat
        beat = signal[i_beats[i]:i_beats[i]+beat_int]   
        beats[i,:len(beat)] = beat
        if padDC:
            # Set rest of beat to mean of beat (in case signal has DC)
            beats[i,len(beat):] = np.mean(beat)

    # Handle the last beat (which has no associated interval)
    # Just use last beat interval for this one (likely similar)
    last_beat_int = intervals[-1]
    beats[-1, :last_beat_int] = signal[i_beats[-1]:i_beats[-1]+last_beat_int]
    if padDC:
        # Append DC value to end of last beat
        beats[-1, last_beat_int:] = np.mean(signal[i_beats[-1]:i_beats[-1]+last_beat_int])

    return beats.T, i_beats


def get_beats_fixed_length(signal, i_beats, intervals, Fs=FS_RESAMPLE, beat_len=int(0.5*FS_RESAMPLE), padDC=False):
    '''
    Segment a signal using provided indices with zero-padding to a fixed length.

    Parameters
    ----------
    signal : _type_
        _description_
    i_beats : _type_
        _description_
    intervals : _type_
        _description_
    Fs : _type_, optional
        _description_, by default FS_RESAMPLE
    padDC : _type_, optional
        _description_, by default False
        
    Returns
    -------
    _type_
        _description_
    '''
    # Compute beat interval into samples
    intervals = np.array(intervals*Fs, dtype=int)
    # Instantiate beats array (N beats x M samples per beat)
    beats = np.zeros(shape=(len(i_beats), beat_len), dtype=float)

    for i, beat_int in enumerate(intervals):
        # Check if beat interval is shorter than specified beat_len
        if beat_int < beat_len:
            # Only store the duration of the beat interval so as not to include info from the next beat
            beat = signal[i_beats[i]:i_beats[i]+beat_int]
        else:
            # Get the portion of the signal related to this beat
            beat = signal[i_beats[i]:i_beats[i]+beat_len]   
        beats[i,:len(beat)] = beat
        # If we want to pad signals with the mean rather than zeros
        if padDC:
            beats[i,len(beat):] = np.mean(signal[i_beats[i]:i_beats[i]+beat_int])

    # Handle the last beat (which has no associated interval)
    # Just use last beat interval for this one (likely similar)
    last_beat_int = intervals[-1]
    if last_beat_int < beat_len:
        beats[-1, :last_beat_int] = signal[i_beats[-1]:i_beats[-1]+last_beat_int]
    else:
        beats[-1, :last_beat_int] = signal[i_beats[-1]:i_beats[-1]+beat_len]

    # Pad DC on last beat
    if padDC:
        beats[-1, last_beat_int:] = np.mean(signal[i_beats[-1]:i_beats[-1]+last_beat_int])

    return beats.T, i_beats




def plot_ALL_beats(beats_dict, beats_id, beats_colors, subject_id, outputdir=None, show_plot=False):

    t_beat = np.arange(beats_dict['ecg'].shape[0])/Fs
    
    fig = plt.figure(figsize=(16, 10), dpi=80)
    fontsize = 20
    alpha = 0.005
    ymin = []
    ymax = []

    for beat_name, beats in beats_dict.items():
        if beat_name == 'i_R_peaks':
            continue
        if 'DC' in beat_name:
            continue
        if 'ppg' in beat_name:
            ymin.append(np.mean(beats,axis=1).min())
            ymax.append(np.mean(beats,axis=1).max())

    ymin = np.min(np.asarray(ymin))
    ymax = np.max(np.asarray(ymax))

    
    for (beat_name, beat_i, color_name) in zip(beats_dict, beats_id, beats_colors):
        if beat_name == 'i_R_peaks':
            continue
        if 'DC' in beat_name:
            continue
            
        
            
        beats = beats_dict[beat_name]
        if 'scg' in beat_name:
            beats *= 1000
            
            
        ax = fig.add_subplot(2, 4, beat_i)
        ax.set_title(beat_name, fontsize=fontsize)

        ax.plot(t_beat, beats, color='gray', alpha=alpha)
        ax.plot(t_beat, np.mean(beats,axis=1), color=color_dict[color_name], linewidth=3)

        if 'ppg' in beat_name:

            beats_mean = np.mean(beats,axis=1)

            ymin = beats_mean.mean() - beats_mean.std()*4
            ymax = beats_mean.mean() + beats_mean.std()*4
            ax.set_ylim(ymin, ymax)

        if 'scg' in beat_name:
            ymin = -15
            ymax = 15
            ax.set_ylim(ymin, ymax)
            
        if 'ecg' in beat_name:
            beats_mean = np.mean(beats,axis=1)
            
            ymin = beats_mean.min() - beats_mean.std()*2
            ymax = beats_mean.max() + beats_mean.std()*2
            ax.set_ylim(ymin, ymax)

        ax.tick_params(axis='both', which='major', labelsize=13)
        ax.set_xlim(-0.1,end_offset/Fs+0.1)
        
        
        if 'ppg' in beat_name:
            ax.set_ylabel('μA', fontsize=fontsize-3)
        if 'ecg' in beat_name:
            ax.set_ylabel('V', fontsize=fontsize-3)
        if 'scg' in beat_name:
            ax.set_ylabel('mg', fontsize=fontsize-3)

        
    ax.set_xlabel('time (sec)', fontsize=fontsize)
    
    
    fig_name = 'beats_ensemble_sub{}'.format(subject_id)
    
    fig.tight_layout()

    if show_plot == False:
        plt.close(fig)
        plt.close('all')

