'''
DEFINE PIPELINE PARAMETERS
'''


# Filter settings
FS_RESAMPLE = 500  # Hz
filter_padded = 5
downsample_factor = 10

FILT_PPG = [1, 8]
FILT_ECG = [0.5, 50] # [3,45] for pyHRV code, #[10, 40] for lab stuff
FILT_SCG = [1, 40]
FILT_PCG = [30, 125]

# Physio-Params
label_range_dict = {
    'br': [4, 45],  # breaths per minute
    'heart_rate': [40, 180],  # beats per minute
    'width_QRS': 0.15  # sec
}
# ECG-specific params
ecg_quality_threshold = 0.5
width_QRS = 0.15
pk_tolerance_ms = 50
offset_SR = int(0.2*FS_RESAMPLE)
end_offset = int(0.2*FS_RESAMPLE)

# PPG and SCG Parameters
ppg_wavelen = 'ppg_ir'
max_dist_foot = int(0.12*FS_RESAMPLE)
# Foot detection range (50-400 ms)
foot_range = [int(FS_RESAMPLE*.05), int(FS_RESAMPLE*.4)]
# foot_range = np.array([int(FS_RESAMPLE*0.01), int(FS_RESAMPLE*.4)])
# SCG Fiducial Settings
max_dist_ao = int(.008*FS_RESAMPLE)
max_dist_ac = int(.025*FS_RESAMPLE)
# AO detection range from start of SCG (0-200 ms)
ao_range = [0, int(FS_RESAMPLE*.2)]
# AC detection range from start of SCG (250-500 ms)
ac_range = [int(FS_RESAMPLE*.25), int(FS_RESAMPLE*.5)]

# Beat parameters
# Specify consistent beat length
beat_length = int(0.7*FS_RESAMPLE)
# Specify ensemble length to average beats
ensemble_length = 15

# # MODIFY TO CONTROL SQI PARAMETERS FOR SCG AND PPG
sqi_PPG_threshold = 0.5
sqi_SCG_threshold = 0.3
strict_thresh = 0.6      
lenient_thresh = 0.3

sqi_flag = False
if sqi_flag:
    ppg_sqi_threshold = 0.65
    scg_sqi_threshold = 0.55
    strict_thresh = 0.6      
    lenient_thresh = 0.3
else:
    ppg_sqi_threshold = 0
    scg_sqi_threshold = 0
    strict_thresh = 0     
    lenient_thresh = 0
# This will be used in place of NaNs to keep the array as type int
INVALID_IDX_CODE = -2