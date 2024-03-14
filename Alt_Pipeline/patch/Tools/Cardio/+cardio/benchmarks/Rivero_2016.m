function [IM,AO,smoothed_SCG] = Rivero_2016(raw_SCG,Fs,plots_on)

% code written by Jacob Kimball on March 1, 2021
% Implementation of Rivero_2016 method of IM and AO extraction from https://ieeexplore.ieee.org/document/7795816

% focus of article: finding AO and IM points in SCG without ECG
% Process: 
%   1. gaussian smoothing filter
%   2. pre-analysis to determine the scale to utilize in the CWT
%   3. rough estimates of points from the CWT and a decision rule
%   4. fine point detection 

if(strcmp(plots_on,'plotsOn')); plotsOn = true; else; plotsOn = false; end

fprintf('Gaussian Smoothing Filter\n'); 
% Gaussian Smoothing Filter: in the paper they used the CEBS dataset with a sampling
% frequency of 5kHz and so used a filter length of 170 samples. 
% calculate the proportional filter length based on Fs: 
filter_length = round((Fs/5000)*170);
smoothed_SCG = smoothdata(raw_SCG,'gaussian',filter_length);
if(plotsOn); figure(); plot(raw_SCG); hold on; plot(smoothed_SCG); end

fprintf('Pre-analysis\n') 
% pre-analysis: analyze coefficients from CWT in the range with scale a =
% [25:25:600], then determine which scale corresponds to the coefficient with highest
% energy a_max and lowest energy a_min. Select the scale between max and
% min: a_sel = (a_max+a_min)/2. They used the db4 wavelet for analysis. 
% The CEBS recordings are only 5 minutes each. Limit the size of the
% coefficients matrix so you don't overload your system 
sig_length = length(smoothed_SCG); 
length_limit = min(sig_length,Fs*60*5); % limit to 5 minutes
a = 25:25:600; % set scale range 
CWTcoeffs = cwt(smoothed_SCG(1:length_limit),a,'db4'); % find coefficients

% calculate energy of coefficients
coeff_energy = sum((CWTcoeffs.^2)')'; 
[~,max_ind] = max(coeff_energy); 
a_max = a(max_ind); % find max energy scale
[~,min_ind] = min(coeff_energy); 
a_min = a(min_ind); % find min energy scale 

a_sel = round((a_max+a_min)/2); % choose scale 

fprintf('Rough Point Detection\n'); 
% rough estimation has 5 steps: 
% 1. calculate CWT coefficients for the first 1.3 second fragment of the signal,
% corresponding to the cardiac period of a healthy subject (0.6-1.2 sec).
% The detected point in the SCG corresponds to the highest coefficient value.
% 2. determine the next fragment to analyze. This fragment has a duration
% of 0.8 seconds after 0.5 seconds from the point previously detected. 
% 3. calculate the CWT of the next fragment 
% 4. find the highest coefficient absolute value
% 5. implement decision rule for correct point detection
%   a. find the average distance between all detected points in the SCG
%   b. for each detected point, if the distance to the prior point differs
%   from the average by more than 20%, run a correction analysis: 
%       i. find the highest coefficient value in a 0.16 second window
%       centered on the detected point
%       ii. compare this value to the value of the original detected point.
%       If it is more than 70% of the value of the detected point, choose
%       the new value to be the detected value. 

% 1. 
CWTcoeffs = cwt(smoothed_SCG(1:round(Fs*1.3)),a_sel,'db4'); % find coefficients for the first 1.3 seconds
[~,max_ind] = max(abs(CWTcoeffs)); % choose detected point as max abs coefficient location
det_pt = [max_ind]; % start the list 

%2. 
half_sec = round(Fs*0.5); % number of samples in 0.5 seconds
zero8_sec = round(Fs*0.8); % number of samples in 0.8 seconds
frag_start = det_pt(end)+half_sec; % next fragment starts 0.5 seconds after prior detected point
frag_end = frag_start+zero8_sec; % fragment has a duration of 0.8 seconds 

while(frag_end <= sig_length) % detect points through the entire signal 
    % 3. 
    CWTcoeffs = cwt(smoothed_SCG(frag_start:frag_end),a_sel,'db4'); % find coefficients for the fragment
    % 4. 
    [~,max_ind] = max(abs(CWTcoeffs)); % choose detected point as max abs coefficient location
    det_pt = [det_pt max_ind+frag_start-1]; % add to the list 
    
    % 2. (repeat until you've gone through the whole signal) 
    frag_start = det_pt(end)+half_sec; % next fragment starts 0.5 seconds after prior detected point
    frag_end = frag_start+zero8_sec; % fragment has a duration of 0.8 seconds 
end

if(plotsOn); figure(); plot(smoothed_SCG); hold on; plot(det_pt,smoothed_SCG(det_pt),'or'); end

% 5. 
zero08 = round(Fs*0.08); % 0.08 seconds 
ave_diff = diff(det_pt)/(length(det_pt)-1); % average difference between detected points
for i = 2:length(det_pt) 
    interval = det_pt(i)-det_pt(i-1); 
    if(abs(interval-ave_diff)>0.2*ave_diff) % 20% threshold 
        CWTcoeffs = cwt(smoothed_SCG(det_pt(i)-zero08:det_pt(i)+zero08),a_sel,'db4'); % find coefficients for the fragment
        [~,max_ind] = max(abs(CWTcoeffs)); % choose detected point as max abs coefficient location
        % skipping the 70% check. just taking the max. 
        det_pt(i) = max_ind+det_pt(i)-zero08-1; % replace with new value 
    end
end

if(plotsOn); plot(det_pt,smoothed_SCG(det_pt),'xk'); end % plot adjusted points
        

% Fine point detection 
fprintf('Fine Point Detection\n'); 
half_win = Fs*0.1; % half window size for window of 0.2 seconds 
AO = zeros(1,length(det_pt)); % placeholder
IM = zeros(1,length(det_pt)); % placeholder 
for i = 1:length(det_pt)
    if(det_pt(i)-half_win<1) % check if it's too low 
        window = 1:det_pt(i)+half_win; % fix the minimum
    elseif(det_pt(i)+half_win>sig_length)
        window = det_pt(i)-half_win:sig_length; % fix the maximum 
    else
        window = det_pt(i)-half_win:det_pt(i)+half_win; 
    end
    [~,max_ind] = max(smoothed_SCG(window)); 
    AO(i) = max_ind+det_pt(i)-half_win-1; 
end

figure(); plot(smoothed_SCG); hold on; plot(AO,smoothed_SCG(AO),'or'); 

IM = 0; % don't care about IM points right now. Future: find min immediately before AO 

% limitations and things to consider: 
%   1. calibration period for CWT. With signals longer than 5 minutes (ie
%   for continuous monitoring, is it sufficient to calibrate this CWT
%   method with the first 5 minutes of data to find the right scale? Or is
%   this something that needs to occur over time as different physiological
%   states may require different scales? 
%   Additional assumption: db4 wavelet is best for all cases 
%   2. double check why the scale range is limited to 600. Also, when will
%   the energy not simply monotonically increase? 
%   3. This skips over the diastolic portion of the signal, and assumes that the
%   AO power is always higher than AC point
%   Also - many assumptions in choosing the next window 
%   4. the correction scheme requires computing the average beat distance
%   over the entire recording. (assuming that the average will be the correct 
%   interval. What if there are a lot of mistakes?) 
%   Also not sure about the 70% cutoff or 0.16
%   second windows.If centered on the prior detected point, how frequently 
%   will this find a different point? For online processing, this could be switched to an
%   exponential moving average. 



end