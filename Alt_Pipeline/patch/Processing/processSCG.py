# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 10:52:50 2022

@author: jaber
"""

import numpy as np
import scipy as sc
import os
import scipy.signal as scsignal
import matplotlib.pyplot as plt
from ..Tools.pipelineConfig import *
from ..Tools.preprocessing import *
from ..Tools.segmentation import *
from ..Tools.signalQuality import *
from ..Tools.dtfm import *
import matlab.engine
from emd import sift
# from multiprocess import Pool



def get_matlab_engine():
    '''
    Start Matlab engine to enable SCG fiducial detection via previously developed Matlab scripts.

    Returns
    -------
    _type_
        _description_
    '''
    eng = matlab.engine.start_matlab()
    # Take matlab.engine to the correct folder where the Matlab scripts reside
    #eng.cd("../patch/Tools/", nargout=0)
    # for vs code it is in reference to where the folder is located 
    eng.cd("/home/michael/Code/mims-transformer-stress-classification/patch/Tools", nargout=0)
    # Ensure that all of the necessary functions are accessible
    eng.addpath(eng.genpath(eng.pwd()+"/Cardio"))
    eng.addpath(eng.genpath(eng.pwd()+"/Smart"))
    eng.addpath(eng.genpath(eng.pwd()+"/scgFeat_extraction"))
    
    #eng.addpath(eng.genpath(eng.pwd()+"/scgFeat_extraction/"))

    return eng


def extract_SCG_fiducials(scg_dict, num_features_ao=4, num_features_ac=3, ema_factor=1.0, verbose=False, method = 'Jon'):
    '''
    AO and AC fiducials detection utilizing the processSCG function in the Matlab Cardio repo, uses GMM to 
    track fiducials in the presence of transients and noise.

    Parameters
    ----------
    scg_dict : _type_
        _description_
    num_features : int, optional
        _description_, by default 4
    win_len : int, optional
        _description_, by default 0
    ema_factor : float, optional
        _description_, by default 0.5
    verbose : bool, optional
        _description_, by default False

    Returns
    -------
    _type_
        _description_
    '''
    # Only use beats that passed the SQI
    beats = scg_dict['scg_beats']
    t_beats = scg_dict['t_beats']
    pcg_beats = scg_dict['pcg_beats']
    beats_used = 'ref_beats'

    ao_cands = np.zeros((num_features_ao, len(t_beats)))
    ac_cands = np.zeros((num_features_ac, len(t_beats)))
    # Get the matlab engine
    eng = get_matlab_engine()
    if verbose:
        print('Engine started...')

    if method == 'Jon':
        # We're going to do all beats at once, hopefully this isn't prohibitively time consuming to run
        # For AO detection, we're going to limit this to the expected search area
        # Note, the copy at the end is to keep the arrays contiguous 
        # (https://stackoverflow.com/questions/26778079/valueerror-ndarray-is-not-c-contiguous-in-cython)
        ao_beats = beats[ao_range[0]:ao_range[1],:].copy(order='F')
        # Setup logic to catch exception when too many candidates are trying to be found, try again with less
        keep_trying = True
        while keep_trying:
            try:
                ao_vals = eng.cardio.scg.processSCG(ao_beats, num_features_ao, ema_factor, 5, False, nargout=1)
            except eng.MatlabExecutionError as e:
                print("MatlabExecutionError:", e)
                num_features_ao -= 1
                if verbose:
                    print(f"Decreasing number of AO candidates to {num_features_ao}...")
                if num_features_ao == 0:
                    if verbose:
                        print("No AO points detected...")
                    keep_trying = False
            else:
                # AO points found successfully
                keep_trying = False
        if verbose:
            print('AO points found...')
        # Same for AC detection, from starting search range to end of beat
        ac_beats = beats[ac_range[0]:ac_range[1],:].copy(order='F')
        keep_trying = True
        while keep_trying:
            try:
                ac_vals = eng.cardio.scg.processSCG(ac_beats, num_features_ac, ema_factor, 5, False, nargout=1)
            except: 
                num_features_ac -= 1
                if verbose:
                    print(f"Decreasing number of AC candidates to {num_features_ac}...")
                if num_features_ac == 0:
                    if verbose:
                        print("No AC points detected...")
                    keep_trying = False
            else:
                # AC points found successfully
                keep_trying = False
        if verbose:
            print('AC points found...')
        # Get the candidates from the return value
        ao_cands = np.squeeze(ao_vals['features']) + ao_range[0]
        ac_cands = np.squeeze(ac_vals['features']) + ac_range[0]

    elif method == 'Template':
        ao_max_locs, ao_min_locs = eng.general.slideTemplate(beats, 30.0, 30.0, float(FS_RESAMPLE), 
                                        "NumFeatures", num_features_ao, "Axrange", np.array(ao_range).astype(float), nargout = 2)
        ao_cands_template = np.vstack([ao_max_locs, ao_min_locs])
        ao_cands = ao_cands_template[np.argsort(ao_cands_template[:, 0]), :]
        
        ac_max_locs, ac_min_locs = eng.general.slideTemplate(beats, 30.0, 30.0, float(FS_RESAMPLE), 
                                        "NumFeatures", num_features_ac, "Axrange", np.array([200, 500]).astype(float), nargout = 2)
        ac_cands_template = np.vstack([ac_max_locs, ac_min_locs])
        ac_cands = ac_cands_template[np.argsort(ac_cands_template[:, 0]), :]
        
    elif method == 'GMM': 
        ao_vals = eng.sortPeaks(beats, float(FS_RESAMPLE), 'axrange', np.array([1, 200]).astype(float), 
                            'brange',np.array([1, 30]).astype(float), 'promconst', 1.0, 'qconst', 2.0, 'influence', 1/10)
        ao_cands = np.squeeze(ao_vals).astype(int).T
        ac_vals = eng.sortPeaks(beats, float(FS_RESAMPLE), 'axrange',np.array([200, 500]).astype(float), 
                            'brange',np.array([1, 50]).astype(float), 'promconst', 3.0, 'qconst', 2.0, 'influence', 1/10)
        ac_cands = np.squeeze(ac_vals).astype(int).T
        
        
    # Kill engine
    eng.quit()
    # Use PCG to select the closest physiological AO and AC point
    # Get PCG envelope by interpolating from the peaks of the signal
    pcg_beat = np.mean(pcg_beats, axis=1)
    pcg_pks,_ = sc.signal.find_peaks(pcg_beat)
    # This creates an interpolater, just need to pass it a new independent variable 
    cs = sc.interpolate.CubicSpline(pcg_pks, pcg_beat[pcg_pks])
    # This is to create the envelope with same length and sampled intervals as the beat
    i_pcg = np.arange(len(pcg_beat))
    pcg_envelope = cs(i_pcg)
    # pcg_envelope = np.abs(sc.signal.hilbert(np.mean(pcg_beats, axis=1)))
    # Take nearest peak after S1 peak to be the AO point
    i_s1_pk = np.argmax(pcg_envelope[ao_range[0]:ao_range[1]]) + ao_range[0]
    # Take nearest peak before S2 peak to be the AC point
    i_s2_pk = np.argmax(pcg_envelope[ac_range[0]:]) + ac_range[0]

    # Get average AO and AC locations across the candidates to decide which is best
    # nanmean as NaNs haven't yet been removed
    if np.shape(ao_cands)[0] > 1:
        avg_ao_cands = np.nanmean(ao_cands, axis=1)
    else:
        avg_ao_cands = ao_cands
    # Add weight to which candidate we think is closest, average PEP is ~100ms
    ao_weight = np.abs(1/(avg_ao_cands - 50))
    # Use the variance in the locations to decide the best AC candidate
    if len(np.shape(ac_cands)) == 1:
        ac_weight = np.abs(1/(np.nanstd(ac_cands)))
    else:
        ac_weight = np.abs(1/(np.nanstd(ac_cands,axis=1)))
    
    # Get the index of the vector of AO points closest to and after the S1 peak
    ao_idx = np.where(avg_ao_cands > i_s1_pk, avg_ao_cands - i_s1_pk, np.inf).argmin()
    # Get the index of the vector of AC points closest to and before the S2 peak
    # ac_idx = np.where(avg_ac_cands < i_s2_pk, i_s2_pk - avg_ac_cands, -np.inf).argmax()
    # Just use the AC with the least variation, it's location is more difficult anyways
    ac_idx = np.argmax(ac_weight)
    
    if ao_weight[ao_idx] < np.sort(ao_weight)[-2]:
        ao_idx = np.argmax(ao_weight)

    # if ac_weight[ac_idx] < np.sort(ac_weight)[-2]:
    #     ac_idx = np.argmax(ac_weight)

    # Select the chosen ones
    if len(np.shape(ao_cands)) == 1:
        ao_locs = ao_cands
    else:
        ao_locs = np.asarray(ao_cands[ao_idx,:])

    if len(np.shape(ac_cands)) == 1:
        ac_locs = ac_cands
    else:
        ac_locs = np.asarray(ac_cands[ac_idx,:])
    # Get indices of not NaN values
    flagged_ao = np.argwhere(~np.isnan(ao_locs))
    flagged_ac = np.argwhere(~np.isnan(ac_locs))

    t_ao = t_beats[flagged_ao]
    t_ac = t_beats[flagged_ac]

    return {'ao_locs':ao_locs,
            'ac_locs':ac_locs,
            'flagged_ao': flagged_ao,
            'flagged_ac': flagged_ac,
            't_ao': t_ao,
            't_ac': t_ac,
            'beats_used': beats_used,
            'ao_cands':ao_cands,
            'ac_cands':ac_cands,
            'pcg_env': pcg_envelope,
            's1_loc': i_s1_pk,
            's2_loc': i_s2_pk,
    }


def process_SCG(scg, pcg, ecg_dict, ensembleSize=ensemble_length, useTemplateSQI=False,\
                       beat_len=beat_length, Fs=FS_RESAMPLE, verbose=True):
    
    # Get start of NN intervals for segmentation
    i_seg = ecg_dict['i_beats']
    
    if verbose:
        print('Begin processing SCG...')
    
    # Separate scg signal into beats w/ R peaks
    scg_beats, i_beats = get_beats_fixed_length(scg, i_seg, ecg_dict['NN_int'], Fs=Fs, beat_len=beat_len, padDC=False)
    # Separate pcg signal into beats w/ R peaks
    pcg_beats, i_pcgbeats = get_beats_fixed_length(pcg, i_seg, ecg_dict['NN_int'], Fs=Fs, beat_len=beat_len, padDC=False)
    # Variable to keep track of which SCG beats (with reference to i_beats) are kept
    ref_beats = np.empty(1)
    sqi_mask = np.zeros(len(i_beats), dtype=bool)

    # Remove SCG beats with outlying amplitudes
    scg_amp = np.zeros(np.shape(scg_beats)[1])
    
    for i in range(len(scg_amp)):
        scg_amp[i] = max(scg_beats[:,i]) - min(scg_beats[:,i])
    
    # Designate outliers based on global median amplitudes
    i_keep_s1, scg_amp_s1 = rm_outliers_mad(scg_amp, num_MAD=5)
    
    scg_beats_s1 = scg_beats[:,i_keep_s1]
    ref_beats = np.argwhere(i_keep_s1)
    if verbose:
        print('Stage 1: Global MAD removed ' + str(np.shape(scg_beats)[1]-np.shape(scg_beats_s1)[1]) + ' out of ' + str(np.shape(scg_beats)[1]))
            
    # Designate outliers based on moving window median amplitudes
    win_len = 30
    overlap = 29/30
    _,i_window = sliding_window(i_beats[ref_beats], win_len, overlap)
    i_keep_s2, scg_amp_s2 = rm_outliers_mad(scg_amp_s1, num_MAD=5, i_slide_win=i_window)
    
    scg_beats_s2 = scg_beats_s1[:,i_keep_s2]
    ref_beats = ref_beats[i_keep_s2]
    ref_beats_amp = ref_beats

    if verbose:
        print('Stage 2: Moving window MAD removed ' + str(np.shape(scg_beats_s1)[1]-np.shape(scg_beats_s2)[1]) + ' out of ' + str(np.shape(scg_beats_s1)[1]))
    
    # If we want to use template-based SQI
    if useTemplateSQI:

        # Normalize SCG beats
        for i in range(np.shape(scg_beats_s2)[1]):
            scg_beats_s2[:,i] = (scg_beats_s2[:,i] - np.mean(scg_beats_s2[:,i]))/np.std(scg_beats_s2[:,i])
        
        # This function uses DTFM (dynamic time feature matching) - a further
        # constrained version of DTW (dynamic time warping) - to compute SQI values
        # for SCG beats (tailored for SCG, as the number of features can vary beat
        # to beat, unlike other signals such as PPG)
        # Compute SQI with 30 second windows
        win_len = 30
        overlap = 0
        dtfm_lambda = 50   # Set as in Asim's work to increase separation between SQI scores
        _,i_window = sliding_window(i_beats[ref_beats], win_len, overlap)
        
        sqi_dtfm = np.zeros(np.shape(scg_beats_s2)[1])
        for i in range(len(i_window)):
            curr_beats = scg_beats_s2[:,i_window[i][0]:i_window[i][-1]]
            # curr_template = create_template_woodys(curr_beats)
            curr_template = np.mean(curr_beats, axis=1)
            sqi_dtfm[i_window[i][0]:i_window[i][-1]] = score(curr_beats, curr_template, maxDist=int(.025*Fs), ell=dtfm_lambda)
            # sqi_dtfm[i_window[i][0]:i_window[i][-1]] = get_sqi(curr_beats, get_dtw, curr_template)

        sqi_dtfm = sqi_dtfm / max(sqi_dtfm)
        sqi_stg1 = sqi_dtfm
        
        sqi_stg1_mask = (sqi_stg1 >= sqi_SCG_threshold)
        scg_beats_s3 = scg_beats_s2[:,sqi_stg1_mask]
        
        if verbose:
            print('Stage 3: Moving window DTFM removed ' + str(np.shape(scg_beats_s2)[1]-np.shape(scg_beats_s3)[1]) + ' out of ' + str(np.shape(scg_beats_s2)[1]))
        
        # SQI Stage 2   
        # Fixed template based off of the beats kept after moving window DTFM removal
        # template = create_template_woodys(scg_beats_s3)
        template = np.mean(scg_beats_s3, axis=1)
        sqi_dtfm2 = score(scg_beats_s2, template, maxDist=int(.025*Fs), ell=dtfm_lambda)
        # sqi_dtfm2 = get_sqi(scg_beats_s2, get_dtw, template)

        sqi_dtfm2 = sqi_dtfm2 / max(sqi_dtfm2)
        sqi_stg2 = sqi_dtfm2
        sqi_stg2_mask = (sqi_stg2 >= sqi_SCG_threshold)
        scg_beats_s4 = scg_beats_s2[:,sqi_stg2_mask]
        
        # for i in range(len(sqi_stg2)):
        #     # If the stage 2 SQI is super low and we previously marked this beat
        #     # as good with the binary label 1
        #     if (sqi_stg2[i] < lenient_thresh) & (sqi_stg1[i] > 0.5):
        #         sqi_stg2_mask[i] = False
            
        #     # On the other hand, if stage 2 SQI is super high and we previously
        #     # marked this beat as bad with the binary label 0
        #     elif (sqi_stg2[i] >= strict_thresh) & (sqi_stg1[i] < 0.5):
        #         sqi_stg2_mask[i] = True
    
        if verbose:
            print('Stage 4: Global DTFM removed ' + str(np.shape(scg_beats_s2)[1]-np.shape(scg_beats_s4)[1]) + ' out of ' + str(np.shape(scg_beats_s2)[1]))
        
        # Final check involves making sure we dealt with the more prominent
        # problem of false rejections. If we are now rejecting even more beats
        # than we were previously, then the single template is likely not good
        # enough to represent all good beats - there may have been a shift in
        # morphology midway that will now cause the final Woody's template to be a
        # bit muddled, reducing stage 2 SQI values 
        # In this case, we will resort back to our original single stage SQI
        
        # So, only assign the stg2 result to the final good v bad returned if
        # stage 2 resulted in more accepted beats (again, our more prominent
        # problem was false rejections; granted, we now forego fixing any false
        # acceptances...but life is full of tradeoffs)
        if (np.sum(sqi_stg1_mask) <= np.sum(sqi_stg2_mask)) and (len(sqi_stg1_mask)-np.sum(sqi_stg1_mask) > 0.1*len(sqi_stg1_mask)):
            scg_beats_good_mask = sqi_stg2_mask
        else:
            scg_beats_good_mask = sqi_stg1_mask

        # Use the normalized beats, neglecting some amplitude information
        scg_beats_s5 = scg_beats_s2[:,scg_beats_good_mask]
        # Return to s1 beats as s2 beats have been normalized
        # scg_beats_s5 = scg_beats_s1[:,i_keep_s2][:, scg_beats_good_mask]
        ref_beats = ref_beats[scg_beats_good_mask]
        sqi_mask[ref_beats] = 1  

    # No sqi
    else:
        scg_beats_s5 = scg_beats_s2

    if verbose:
        print(f'Final: Using {np.shape(scg_beats_s5)[1]} beats out of {np.shape(scg_beats)[1]}')

    # Get the PCG beats based on the sqi mask for SCG
    og_pcg_beats = np.squeeze(pcg_beats[:, ref_beats_amp])
    pcg_beats = np.squeeze(pcg_beats[:, ref_beats])
    # Note, ref_beats is now in reference to ecg_dict['i_segment'], where we want it in reference
    # to ecg_dict['i_beats'] for creating a clean feature matrix/dataframe later
    _, ref_beats,_ = np.intersect1d(ecg_dict['i_beats'], i_beats[ref_beats], return_indices=True)
    _, ref_beats_amp,_ = np.intersect1d(ecg_dict['i_beats'], i_beats[ref_beats_amp], return_indices=True)

    if ensembleSize > 1:
        scg_beats_s5 = rolling_ensemble(scg_beats[:,ref_beats], win_len=ensembleSize)
        scg_beats = rolling_ensemble(scg_beats[:,ref_beats_amp], win_len=ensembleSize)

        if len(np.shape(pcg_beats)) == 1:
            pcg_beats = rolling_ensemble(pcg_beats[:,np.newaxis], win_len=ensembleSize)
        else:
            pcg_beats = rolling_ensemble(pcg_beats, win_len=ensembleSize)
    
    # Add EMD step to denoise the SCG beats, taking the first IMF as the denoised SCG
    scg_emd_beats = np.zeros(shape=np.shape(scg_beats))
    for i in range(np.shape(scg_beats)[1]):
        imf = sift.sift(scg_beats[:,i])
        # Use only first IMF as 'denoised' SCG beat
        scg_emd_beats[:,i] = imf[:,0]

    scg_emd_beats_s5 = np.zeros(shape=np.shape(scg_beats_s5))
    for i in range(np.shape(scg_beats_s5)[1]):
        imf = sift.sift(scg_beats_s5[:,i])
        # Use only first IMF as 'denoised' SCG beat
        scg_emd_beats_s5[:,i] = imf[:,0]

    # Get times of final set of beats
    t_beats = ecg_dict['t_beats'][ref_beats]

    return {'scg_beats': scg_beats_s5,          # The beats that passed all stages of preprocessing
            't_beats': t_beats,                 # Times of each SCG beat that passed preprocessing
            'ref_beats': ref_beats,             # Mask on ECG R peaks locations that correspond to the final set of SCG beats
            'ref_beats_amp': ref_beats_amp,     # Mask on ECG R peaks for the set of SCG beats after amplitude-based outlier removal
            'pcg_beats': pcg_beats,             # Beat segmented Phonocardiogram for fiducial detection
            'og_pcg_beats': og_pcg_beats,       # The PCG beats before SQI, but after amplitude-based outlier removal
            'og_beats': scg_beats,              # The SCG beats before SQI, but after amplitude-based outlier removal
            'og_emd_beats': scg_emd_beats,      # Denoised SCG beats, only beats that passed amplitude-based outlier removal
            'scg_emd_beats': scg_emd_beats_s5,  # Denoised SCG beats, beats that passed all stages of preprocessing
            }