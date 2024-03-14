function [AO, recon_signal] = Choudhary_2019(raw_scg, fs, std_const, preprocess, Plot)
%%% AO detection algorithm
% Base Alg: T. Choudhary. et al. "Automatic Detection of Aortic Valve Opening Using
% Seismocardiography in Health Individuals". IEEE Journal of Biomeidcal and
% Health Informatics. 2019. 
%
% Requirements: None
% Inputs: 
% raw_sig:          scg signal 
% preprocess:       Apply preprocessing as done in paper. Not always
%                       recommended if recon_signal is to be used later
% std_const         Used to threshold and discriminate between AO regions
%                       and non-AO regions (suggested value: 2 - paper)
%                       lower if you see that AO regions are being skipped
% fs:               sampling frequency (Hz)
% Plot:             flag to plot or not 
%
% Goal: 
%       Find AO in the SCG signal without ECG via DWT and envelope
%       detection
% Framework: 
%       A. Split the signal into 10 second non-overlapping intervals 
%       B. Preprocessing (normalization and median filter)
%       C. DWT Enhancement of SCG signal. Choose subbands based off of
%       kurtosis and frequency of each subband and do a relative weighted
%       reconstruction
%       D. Detect AO regions by calculating shannon energy of signal. Then
%       smoothening the signal via autocorrelation. Threshold the signal
%       and then use first order gaussian derivative (FOGD) filters to find
%       the peaks of the autocorrelation signal (i.e the general location
%       of the AO regions). Finally use a finer window around the peaks to
%       find the actual AO points.
% Usage: 
%       Choudhary(scg_signal, Rpeak_locs, 2000, 3, 1)
%%
%%%% A: Segmentation
%%%%
%%%%
% segementation into 10 second nonoverlapping windows  
segment_wind = 10*fs;   

% calculate window bounds 
num_winds = floor(length(raw_scg)/segment_wind);
wind_bnds = [[0:(num_winds-1)]*segment_wind + 1; [1:(num_winds)]*segment_wind];

% vector for plotting 
Plot_vec = zeros(num_winds, 1);

if Plot 
    rand_beat = randi(num_winds);
    Plot_vec(rand_beat) = 1;
end

% hold the AO estimates
AO_points_all = [];
AO_enhanced_sig = zeros(size(raw_scg));

% placeholder for plotting 
i = 1;
% for each window: 
% 1. Preprocess 
% 2. DWT Enhancement 
% 3. AO Detection
for wind_i = wind_bnds
     
    % window the signal 
    raw_sig_seg = raw_scg(wind_i(1):wind_i(2));   
    
    %%%% B: preprocessing stepd:
    %%%% 1. normalized signal
    %%%% 2. 5th order median filter 
    if preprocess
        raw_sig_seg = raw_sig_seg - mean(raw_sig_seg);
        raw_sig_seg = raw_sig_seg/max(abs(raw_sig_seg));
        raw_sig_seg = medfilt1(raw_sig_seg, round(5*fs/1000));
    end
    
    
    %%%% C: apply DWT to accent AO and remove baseline
    %%%%
    %%%%
    wave = 'bior3.9';
    recon = weighted_dwt(raw_sig_seg, wave, fs, Plot_vec(i));
    
    
    %%%% D: detect AO points from the DWT filtered segment 
    %%%%
    %%%%
    AO_points = detect_AO(recon, fs, std_const, Plot_vec(i));
    
    i = i + 1;
    
    % Save the outputs 
    AO_points_all = [AO_points_all; AO_points + wind_i(1) - 1];
    AO_enhanced_sig(wind_i(1):wind_i(2)) = recon;
end

AO = AO_points_all;
recon_signal = AO_enhanced_sig; 

if Plot
    figure;
    ax(1) = subplot(2, 1, 1);
    plot(raw_scg); title('Original Signal')
    ax(2) =  subplot(2, 1, 2);
    plot(AO_enhanced_sig); hold on; scatter(AO_points_all, AO_enhanced_sig(AO_points_all), 20,  'filled'); 
    title('Reconstructed Signal with AO Points Marked')
    linkaxes(ax, 'x');
end

end
%%


function recon_sig = weighted_dwt(noisy_sig, wave, fs, Plot)
% Function: Apply DWT and select bands using kurtosis and frequency. 
% reconstruct the signal using a weighted sum of the subbands 
% Inputs: 
% noisy_sig        (dbl)   Signal to enhance
% wave          (string)   Wavelet used
% fs                (Hz)   sampling frequency

%%% Wavelet decomp and reconstruction by suband 
% 1. Decompose using the mother wavelet (biorothogonal 3.9) 9 lvls
% 2. reconstruct each subband 

% decompose signal into its subbands                    %1.
decomp_lvls = 9;
[c, l] = wavedec(noisy_sig, decomp_lvls, wave); 

% reconstruct each detail and approx subband            %2.
subbnds = zeros(size(noisy_sig, 1), decomp_lvls+1);
for lvl = 1:decomp_lvls
    
    % detail band reconstruction
    subbnds(:, lvl) = wrcoef('d', c, l, wave, lvl);
    
    if lvl == decomp_lvls
        
        subbnds(:, lvl) = wrcoef('d', c, l, wave, lvl);
        
        % approx band reconstruction
        subbnds(:, lvl+1) = wrcoef('a', c, l, wave, lvl);
    end
end

%%% Calculate decision features (Dominant Multiscale Kurtoisis (DMK), Dominant
%%% Multiscale Central Fequency (DMCF)) and extract the bands that follow a
%%% certain criteria. Reconstruct the qualified bands through a weighted sum 
%%% based on kurtosis
% 1. Calculate DKM and DMCF 
% 2. Get the subbands depending on DMK and DMCF 
%           kurtosis >= 4, 6 < central_freq < 45
% 3. Weight by kurtosis to accent AO based on relative squared kurtosis
% 4. Reconstruct the denoised signal by summing the weighted bands 

% DMK 
dmk = kurtosis(subbnds);

% DMCF using power spectrum 
subfreq = periodogram(subbnds,rectwin(size(subbnds, 1)),size(subbnds, 1),fs);
f_range = 0:fs/size(subbnds, 1):fs/2;
dmcf = sum(f_range'.*subfreq)./sum(subfreq);

% get the subbands that follow a criteria 
f_decision = [6, 45]; 
k_decision = 4;
f_accept = [dmcf > f_decision(1) & dmcf < f_decision(2)];
k_accept = [dmk >= k_decision];
accept_bands = and(f_accept, k_accept);

% weight by relative squared kurtosis and reconstruct 
dmk_accept = dmk(accept_bands);
rsdmk = (dmk_accept.^2) ./ sum(dmk_accept.^2);
bands_accept = subbnds(:, accept_bands);
recon_sig = sum(bands_accept.*rsdmk, 2);

if Plot
    % get the bands accepted 
    bands = 1:length(k_accept);
    bands = bands(accept_bands);
    figure;
    subplot(sum(accept_bands)+3, 1, 1);
    plot(noisy_sig);
    title('Original Sig')
    for i = 2:sum(accept_bands)+3-2
        subplot(sum(accept_bands)+3, 1, i);
        plot(bands_accept(:, i-1));
        title(['Subband: ', num2str(bands(i-1))])
    end
    subplot(sum(accept_bands)+3, 1, sum(accept_bands)+3-1);
    plot(recon_sig);
    title('Recon Sig')
    subplot(sum(accept_bands)+3, 2, 2*(sum(accept_bands)+3)-1);
    plot(1:length(dmcf), dmcf); xlabel('Decomp level'); ylabel('Frequency')
    hold on; scatter(find(f_accept == 1), dmcf(f_accept), 'filled'); 
    legend('All Bands', 'Accepted Bands')
    subplot(sum(accept_bands)+3, 2, 2*(sum(accept_bands)+3));
    plot(1:length(k_accept), dmk);  xlabel('Decomp level'); ylabel('Kurtosis')
    hold on; scatter(find(k_accept == 1), dmk(k_accept), 'filled'); 
    legend('All Bands', 'Accepted Bands')

end

end
function AOpeaks = detect_AO(recon_sig, fs, std_const, Plot)
% Function: Isolate AO regions through envelope formation and thresholding. 
% Then extract more accurate AO point on 
% Inputs: 
% recon_sig        (dbl)   Signal to enhance
% fs                (Hz)   sampling frequency

%%% AO localization 
% 1. Create envelope on reconstructed signal based on shannon energy
% 2. calculate autocorrelation on windowed shannon energy with FIXED lag
% 3. Smooth the autocorr signal with low pass
% 4. Estimate the periodicity 
% 5. Create a pulse train to isolate AO regions 
% Inputs
% recon_sig  [Nx1]  DWT AO enhanced signal
% fs          (Hz)  Sampling frequency
% std_const   (dbl) Threshold to remove AO regions from autocorrelation
%                       signal
acc_wind_rad = 20;          % autocorrelation window radius 
recon_sig = recon_sig - mean(recon_sig);


% envelope function shannon energy                                  %1.
Senergy = @(x) -1*(x.^2 .* log2(x.^2));
env_sig = Senergy(recon_sig);


% nan corrections 
env_sig(isnan(env_sig)) = 0;

% Autocorrelation on windowed SE                                    %2.
% fixed lag on ACC
lag = 10*fs/1000;

% window length and radius 
wind_len = acc_wind_rad;
wind_rad = round(wind_len/2);

% consecutive window shifts in samples 
wind_shift = 1;

% calculate how mnay windows we need to create
acc_len = floor((length(env_sig) - lag - wind_len - 1)/wind_shift+1);

% calculate the acc signal 
acc_env = zeros(acc_len, 1);
for i = 1:length(acc_env)
    
    % window the envelope signal and its lag by P
    bounds = [1+(i-1)*wind_shift , 1+(i-1)*wind_shift+wind_len];
    w_env = env_sig(bounds(1):bounds(2));
    w_env_shift = env_sig(bounds(1)+lag:bounds(2)+lag);
    
    % calculate autocorr at that time instance 
    acc_env(i) = sum(w_env.*w_env_shift);
end

% Make the autocorr signal the same length as the original signal
acc_env = [zeros(wind_rad, 1); acc_env; zeros(length(env_sig)-length(acc_env)-wind_rad, 1)];

% low pass filter autocorr signal                           %3.
cutoff = 65;
[b, a] = butter(6, cutoff/(fs/2));
acc_sig = filtfilt(b, a, acc_env);

% estimate average beat length from smoothened signal and original signal 
acorr_env = xcorr(acc_sig);
trunc_acorr_env = acorr_env(length(acc_sig):end);
acorr_recon = xcorr(recon_sig);
trunc_acorr_recon = abs(acorr_recon(length(acc_sig):end));
trunc_acorr_env = (trunc_acorr_env - min(trunc_acorr_env))/trapz(trunc_acorr_env - min(trunc_acorr_env));
trunc_acorr_recon = (trunc_acorr_recon - min(trunc_acorr_recon))/trapz(trunc_acorr_recon - min(trunc_acorr_recon));

% calculate the range of periodicities to search around
bpm_range = [150, 40];
int_lens = round(1./bpm_range*60*fs);  
int_range = int_lens(1):int_lens(2);



% threshold by 2* the std of the autocorr signal                     %4. 
threshold = std_const * std(acc_sig);
acc_thresh = acc_sig; acc_thresh(acc_thresh <= threshold) = 0;


%%% Peak Detection
%%% 1. Find peaks of thresholded autocorrelation signal using FOGD filters
%%% 2. Refine AO point using autocorrelation peaks as reference and setting
%%% a window around it and finding the max within that window 


% gaussian parameters for FOGD                                          %1.
I = 2.5*fs; 
alpha = 0.05*fs;
gauss = gausswin(I, alpha);

% FOGD window creation 
d = gauss(1:end-1) - gauss(2:end);

% Apply FOGD window on autocorrelation signal 
fogd = conv(d, acc_thresh);

% readjust the length of the fogd such that the center of the window is
% applied at the start of the signal 
diff_len = length(fogd) - length(recon_sig);
fogd = fogd((1 + diff_len/2):(length(fogd)-diff_len/2));

% now find negative zero corssings for autocorr peaks                   %2. 
zerocross =   fogd(1:end-1).*fogd(2:end) < 0 ;
positives = fogd > 0;
peaks = and([0; zerocross], positives);

% finally better locate peak of each scg ao points using recon signal   %3.



% get the list of peaks from the fogd peak detection 
valid_peaks = find(peaks == 1);

% window to search around the autocorr peak 
search_len = 0.3*fs; 
search_rad = round(search_len/2);
    
% begin looking around each point to find a better AO est 
AOpeaks = zeros(size(valid_peaks));
for i = 1:length(valid_peaks)
    
    % set up the window to look around taking into consideration the bounds
    % of the recon signal, let the peak be in the center of the window 
    loc = valid_peaks(i);
    bounds = [max((loc-search_rad), 1), min((loc+search_rad), length(recon_sig))];
    
    % in edge cases if the peak is not in the center of the search window
    % and is asymmetric with fewer point on the left, recalculate the
    % location of the peak in relation to the lower bounds of the search
    % window 
    if loc - search_rad < 0
        % locs is closer to the left of window 
        % acutal center = expected center - (expected window length - actual window length)
        center = (search_rad+1 - ((search_rad*2+1) - (diff(bounds)+1)));
    else
        % locs is in the center even if window is shifted right
        center = 1+search_rad;
    end

    % find a better estimate of AO using the maximum 
    AOpeaks(i) = find(recon_sig(bounds(1):bounds(2)) == max(recon_sig(bounds(1):bounds(2))), 1, 'first');
    
    % realign the located AO peak to the recon signal 
    AOpeaks(i) = (AOpeaks(i) - center)+loc;
end
if Plot
    figure; 
    subplot(6, 1, 1);
    plot(recon_sig); title('Recon Signal')
    subplot(6, 1, 2); 
    plot(env_sig); title('Shannon Entropy Envelope Sig')
    subplot(6, 1, 3);
    plot(acc_sig); title('Autocorrelation Time Series')
    subplot(6, 1, 4);
    plot(acc_thresh); title('Thresholded Autocorrelation Time Series')
    subplot(6, 1, 5);
    plot(fogd); hold on; scatter(valid_peaks, zeros(length(valid_peaks), 1), 'filled');
    title('FOGD Output with Zeros Marked (corresponds to peaks)')
    subplot(6, 1, 6);
    plot(recon_sig); hold on; scatter(AOpeaks, recon_sig(AOpeaks));
    title('SCG with AO points annotated')
    
    
end
end

% Limitations:
% 1. If AO regions in a 10 seconds window are not similar in energy, the
% lower AO regions will dissapear after thresholding causing more false
% negatives 
% 2. Kurtosis does not really do much and in-band noise will not really be
% filtered out if kurtosis is used. 
% 3. Performance dependent on the wavelet used 
