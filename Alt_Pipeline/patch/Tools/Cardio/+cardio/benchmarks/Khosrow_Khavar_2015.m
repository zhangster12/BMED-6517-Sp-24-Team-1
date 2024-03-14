function [final_ac, final_im] = Khosrow_Khavar_2015(scg_sig, ecgR, fs, envelope_type, Plot)
% Automatic Annotation of Seismocardiogram With High-Frequency Precordial Accelerations
%
% Requirements: Needs cardio library 
% Inputs: 
% scg_sig:          scg signal 
% ecgR:             ecg r peak locations 
% fs:               sampling frequency (Hz)
% envelope_type     uses either (1) absolute value (2) shannon energy, or 
%                       (3) Hilbert transform to generate envelope profile
% Plot:             flag to plot or not 
%
% Goal: 
%       Find IM and AC points in SCG signal using envelope detection and
%       ECG 
% Framework: 
%       A. Split ACC signal in low pass (SCG) and high pass (HFACC) signal
%       B. Calculate the envelope of the high pass signal
%       C. Find the peaks of the envelopes for each beat
%       D. From there test a set of alpha values that dictate the search
%       range from the left of each envelope peak to find the IM and AC
%       points for each beat. Find the first extrema to the left of the
%       boundary dictated by alpha percent of the peak 
%       E. Using the differences in consecutive fiducial points as the cost
%       function determine the most optimal alpha
%       F. For beats that have a high deviation in IM/AC locations, find a
%       better IM/AC around the median of all IM/AC points 
% Usage: 
%       Khorsrow_Khavar_2015(scg_signal, Rpeak_locs, 2000, 3, 1)
% Hyper parameters: 
%       f1, f2: search radius (ms) around median for resampling IM/AC
%       points
%       alpha_range: alpha values to test


% alpha values to test (in the paper it is set from 0.5:0.1:1 but with the
% pig dataset 0.9:0.02:1 showed much better results)
alpha_range = 0.5:0.1:1;


%%%% A: Splitting the signal 
%%%%
%%%%

% now filter the accleration signal to get 2 signals: SCG (< 30Hz) and 
% HFACC ( > 20Hz) with a 5th order Butterworth filter low pass and high 
% pass filter respectively.
hcutoff = 20;
lcutoff = 40;
filt_ord = 5;
[bhigh,ahigh] = butter(filt_ord,hcutoff/(fs/2), 'high'); 
[blow,alow] = butter(filt_ord,lcutoff/(fs/2), 'low'); 

% get a high pass and low pass version of acclerometer signal
scg_hp = (filtfilt(bhigh, ahigh, scg_sig));
scg_lp = (filtfilt(blow, alow, scg_sig));

%%%% B: Calculating the envelope 
%%%%
%%%%

% HFACC envelope calculation: The envelope of the HFACC signal will be used
% to delineate points on the SCG signal.
% 4 envelope calculations are tested in this paper: 
% 1. CSCW (cardiac sound characteristic waveform) (skipped this for now)
% 2. Absolute Envelope: take the average of the absolute value in a 20ms 
%           sliding windows of 20ms 
% 3. Shannon Energy: take the average of the multiplication of the square
%           of signal and the square of the log of signal in a 20ms 
%           sliding window
% 4. Hilbert: take the average of the magnitude of the hilbert transform in 
%           a 20ms sliding window

% absolute value function 
AbsEnv = @(x) abs(x);

% shannon envelope function
ShEnv = @(x) -1* x.^2 .* log2(x.^2);

% hilbert transform function
Hilb = @(x) envelope(x);

% determine which envelope function to use
switch envelope_type
    case 1 
        env_sig = AbsEnv(scg_hp);
    case 2 
        env_sig = ShEnv(scg_hp);
    case 3
        env_sig = Hilb(scg_hp);
    otherwise 
        error('Input 1 (absolute), 2 (shannon envelope), or 3 (hilbert) for the envelope_type arg')
end

% for the shannon envelope method (in case we take log of of 0) 
env_sig(isnan(env_sig)) = 0; 

% now we want to window by 20ms  
wind_len = 20*(fs/1000);
wind_rad = round(wind_len/2);

% how much to shift each window by (right now we are taking 20ms of
% non-overlapping windows)
wind_shift = 2*wind_rad;

%calculate the number of windows to create assuming a centered window: we
%calculate the start of the last window that's guarenteed to fit and
%calculate the number of windows that can fit in this period given the
%number of samples between the start of each window 
mean_len = floor((length(env_sig) - 2*wind_rad - 1)/wind_shift+1);

% allocate the mean envelope signal and the bound
mean_env = zeros(mean_len, 1);
tot_bounds = zeros(mean_len, 2);

for i = 1:length(mean_env)
    
    % window the envelope signal
    bounds = [1+(i-1)*wind_shift , 1+(i-1)*wind_shift+wind_len];
    tot_bounds(i, :) = bounds';
    w_env = env_sig(bounds(1):bounds(2));
    
    % take the mean of the windowed envelope signal 
    mean_env(i) = mean(w_env);

end

% interpolate to make same length as the
mean_env = interp(mean_env, wind_shift);


% pad behind and ahead to not only make it exactly the same length as the
% envelope signal but also make it so the mean of each window corresponds
% to the middle rather than the beginning of the window
mean_env = [zeros(wind_rad, 1); mean_env; zeros(length(env_sig)-length(mean_env)-wind_rad, 1)];

% get the heart rate for each beat 
hr = diff(ecgR);

%%% C: Find the IM and AC point from each envelope peak 
%%%
%%%
% now we have a mean envelope signal now we need to segment by beats and
% find the peaks associated with each envelope
% 1. The IM window is 1 to 200ms 
% 2. The AC window is a regression based on heart rate (read paper)
% AC window equation assume input hr is in samples
% hr is in samples/beat and we want beats/minute so 
% (1/hr)[beat/sample] * (fs) [samples/s] * (60) [s/min]     % now we have bpm
% [ms] * fs [samples/s] * 1/1000 [s/ms]                     % now we have samples
AC_bound = @(hr, fs) round([477-2*(fs*60/hr), 607 - 2*(fs*60/hr)]*(fs/1000));  
IM_bound = round([1, 200*(fs/1000)]);

% now separate by beats 
env_beats = cardio.general.separateBeats([scg_lp, scg_hp, mean_env, scg_sig], 'indices', ...
    round(ecgR), 'samples', min(diff(ecgR)));
env_beats = env_beats(:, :, 1:end-1);

hp_beats = squeeze(env_beats(2, :, :));
lp_beats = squeeze(env_beats(1, :, :));
scg_beats = squeeze(env_beats(4, :, :));
env_beats = squeeze(env_beats(3, : ,:));

% random beat for plotting 
rng_beat = randperm(length(hr), 1);

%%% now for each beat get the peak of each envelope 
env_IM = zeros(length(hr), 2);
env_AC = zeros(length(hr), 2);
for beat = 1:length(hr)
    
    
    IM_wind = IM_bound;
    AC_wind = AC_bound(hr(beat), fs);
    
    % detect the peaks in the 2 windows 
    env_IM_loc = IM_wind(1) + find(env_beats(IM_wind(1):IM_wind(2), beat) ...
        == max(env_beats(IM_wind(1):IM_wind(2), beat)));
    env_AC_loc = AC_wind(1) + find(env_beats(AC_wind(1):AC_wind(2), beat) ...
        == max(env_beats(AC_wind(1):AC_wind(2), beat)));
    
    % save the location and the actual max of the envelope signal 
    env_IM(beat, :) = [env_IM_loc; env_beats(env_IM_loc, beat)];
    env_AC(beat, :) = [env_AC_loc; env_beats(env_AC_loc, beat)];
    

    
    
end

%%%% D: Calculate the IM/AC points given an alpha value and save the
%%%% beat to beat cost 
%%%%

% placeholder for beats
i = 1; 

% This is the limit that identifies which beat has a high cost (i.e it's
% IM/AC location deviates too much from the pervious beat's IM/AC)
limit = 25*(fs/1000);

% cost array for determining optimal alpha  
cost = zeros(length(alpha_range), 2);

% parameters! look for min or max right or left of alpha bounds (for the
% paper this left is always the direction and IM corresponds to min and AC
% corresponds to max)
im_ext = 'min'; im_dir = 'left';
ac_ext = 'max'; ac_dir = 'left';

% for each alpha find the cost using the cycle to cycle differences 
for alpha = alpha_range
    
    % placeholders for the fid points and the cycle to cycle difference 
    im_locs = zeros(length(hr), 1);
    ac_locs = zeros(length(hr), 1);
    im_ccd = zeros(length(hr)-1, 1);
    ac_ccd = zeros(length(hr)-1, 1);
    
    % set the im and ac alpha 
    alpha_AC = alpha;
    alpha_IM = alpha;
    
    % for each beat search for the first extrema outside the bounds of the
    % alpha bounds 
    for beat = 1:length(hr)
        
        [ac_locs(beat), ac_bnd] = alpha_search(lp_beats(:, beat), env_beats(:, beat),...
            env_AC(beat, 1), alpha_AC, fs, ac_ext, ac_dir);
        [im_locs(beat), im_bnd] = alpha_search(lp_beats(:, beat), env_beats(:, beat), ...
            env_IM(beat, 1), alpha_IM, fs, im_ext, im_dir);
        
        % here we calculate cycle to cycle difference 
        if beat > 1
            im_ccd(beat-1) = diff(im_locs(beat-1:beat));
            ac_ccd(beat-1) = diff(ac_locs(beat-1:beat));
        end
        
    end
    
    % NaN condition just force them to have a high cost to penalize
    im_ccd(isnan(im_ccd)) = limit+10;
    ac_ccd(isnan(ac_ccd)) = limit+10;
    
    % calculate ethe cost for each fid point 
    cost(i, :) = [sum(abs(im_ccd)), sum(abs(ac_ccd))];
    
    % incrememnt for the next alpha
    i = i + 1;
end
%%%% E. Find the most optimal alpha that generates the smallest cost for all
%%%% tested alphas 
%%%%


% ok now given an optimal IM/AC alpha once again find the locs associated 
im_locs = zeros(length(hr), 1);
ac_locs = zeros(length(hr), 1);
im_ccd = zeros(length(hr)-1, 1);
ac_ccd = zeros(length(hr)-1, 1);

% best alpha try to get the most strict one because that gives the best
% results 
alpha_AC = alpha_range((cost(:, 2) == min(cost(:, 2)))); alpha_AC = alpha_AC(end);  
alpha_IM = alpha_range((cost(:, 1) == min(cost(:, 1)))); alpha_IM = alpha_IM(end);

% plot the costs for IM and AC 
if Plot 
    figure; 
    subplot(2, 1, 1);
    plot(alpha_range, cost(:, 1)); hold on;
    xlabel('\alpha'); ylabel('Cost')
    title('IM Point Costs')
    scatter(alpha_IM, cost(alpha_range == alpha_IM, 1))
    legend('Costs', 'Optimal IM_{\alpha}')
    subplot(2, 1, 2);
    plot(alpha_range, cost(:, 2)); hold on;
    title('AC Point Costs')
    scatter(alpha_AC, cost(alpha_range == alpha_AC, 2))
    legend('Costs', 'Optimal IM_{\alpha}')
end

% refine the best IM/AC points with the given alpha (doing it very
% unoptimally as you could just save the IM/AC points from all tried alphas
% in the previous loop and then just choose the column that corresponds to
% the best alpha for either IM or AC)
for beat = 1:length(hr)
    
    [ac_locs(beat), ac_bnd] = alpha_search(lp_beats(:, beat), env_beats(:, beat), ...
        env_AC(beat, 1), alpha_AC, fs, ac_ext, ac_dir);
    [im_locs(beat), im_bnd] = alpha_search(lp_beats(:, beat), env_beats(:, beat), ...
        env_IM(beat, 1), alpha_IM, fs, im_ext, im_dir);
    
    % here we calculate cycle to cycle difference 
    if beat > 1
        im_ccd(beat-1) = diff(im_locs(beat-1:beat));
        ac_ccd(beat-1) = diff(ac_locs(beat-1:beat));
    end
    
    if beat == rng_beat && Plot
        f = figure;
        plot(env_beats(:, beat)/max(env_beats(:, beat))*max(hp_beats(:, beat))); hold on;
        plot(hp_beats(:, beat)); plot(lp_beats(:, beat)); 
        
        % plot 
        scatter(im_locs(beat), lp_beats(im_locs(beat), beat), 'filled');
        scatter(ac_locs(beat), lp_beats(ac_locs(beat), beat), 'filled');
        

        xline(im_bnd(1), 'b--'); 
        xline(ac_bnd(1), 'r--');
        scatter(env_IM(beat, 1), env_beats(env_IM(beat, 1), beat)/max(env_beats(:, beat))*max(hp_beats(:, beat)), ...
            'blue', 'filled')
        scatter(env_AC(beat, 1), env_beats(env_AC(beat, 1), beat)/max(env_beats(:, beat))*max(hp_beats(:, beat)), ...
            'red', 'filled')
        
        legend('HFACC Env', 'HFACC', 'SCG',...
            'est IM', 'est AC', 'IM Upper Bound', 'AC Upper Bound', 'IM Env Peak', 'AC Env Peak')
        title("Alpha IM: " + string(alpha_IM) + " Alpha AC: " + string(alpha_AC) + " Beat: " + string(beat));
    end
    
end

%%% F: replace IM/AC points with high costs by resampling those beats and
%%% refining the search using a window around the median of all found IM/AC points 
% first find the median of all the locations 
%           or median of a window of beats
% then find the ones that go above a certain ccd 
%           i guess this eliminates outliers (but only single ones)

poor_IM = find(abs(im_ccd) > limit)+1;
poor_AC = find(abs(ac_ccd) > limit)+1;

% Determine what type of median to do 
%           0. Global median (paper)
%           1. windowed median (works better)
% is used)
windowed = 0; med_rad = 10;

% determine what kind of search to do (in the paper the min is used):
%           0. max/min within a range (paper)
%           1. closest to median (faster)
closest = 0; 

% placeholder for best IM and AC regions 
final_ac = ac_locs;
final_im = im_locs;

% search radius around median 
f1 = 15*(fs/1000); % 15 ms radius search 
f2 = f1; 

% for all beats with bad AC points find a better point 
for beat = 1:length(poor_AC)
    
    % if windowed get the med_rad beats before and after to get a median
    if windowed
        AC_med = round(median(ac_locs(max(1, poor_AC(beat)-med_rad): ...
            min(length(ac_locs), poor_AC(beat)+med_rad))));
    else
        AC_med = round(median(ac_locs));
    end
    
    % now for all beats with high ccd go through their lp beat and find
    % all max/mins
    if strcmp(ac_ext, 'max')
        [~, locs] = findpeaks(lp_beats(:, poor_AC(beat)));
    else
        [~, locs] = findpeaks(-1*lp_beats(:, poor_AC(beat)));
    end
    
    
    if closest 
        
        % find the closest 
        [~, close_med] = min(abs(locs - AC_med));
        final_ac(poor_AC(beat)) = locs(close_med);
        
    else
        
        % find the extrema within a range 
        if strcmp(ac_ext, 'max')
            
            % find all peaks within a range 
            valid_ext = locs(abs(locs-AC_med) < f1);
            
            % if there contains valid peaks
            if ~isempty(valid_ext)
                % find the index assoicated with the largest peak
                final_ac(poor_AC(beat)) = valid_ext(find(lp_beats(valid_ext, poor_AC(beat)) ...
                    == max(lp_beats(valid_ext, poor_AC(beat))), 1, 'first'));
                
            % otherwise just find the closest peak
            else
                [~, close_med] = min(abs(locs - AC_med));
                final_ac(poor_AC(beat)) = locs(close_med);
            end
            
        else
            
            % find all peaks within a range
            valid_ext = locs(abs(locs-AC_med) < f1);
            
            % if there contains valid peaks 
            if ~isempty(valid_ext)
                % find the index associated with the lowest valley      
                final_ac(poor_AC(beat)) = valid_ext(lp_beats(valid_ext, poor_AC(beat)) ...
                    == min(lp_beats(valid_ext, poor_AC(beat))), 1, 'first');
                
            % otherwise just find the closest peak
            else
                [~, close_med] = min(abs(locs - AC_med));
                final_ac(poor_AC(beat)) = locs(close_med);
            end
            
        end
    end
    
end

% for all beats with poor IM points resample to find better ones 
for beat = 1:length(poor_IM)
    
    % if windowed get the med_rad beats before and after to get a median
    if windowed
        IM_med = round(median(im_locs(max(1, poor_IM(beat)-med_rad):...
            min(length(im_locs), poor_IM(beat)+med_rad))));
    else
        IM_med = round(median(im_locs));
    end
    
    % now for all beats with high ccd go through their lp beat and find
    % all max/mins
    if strcmp(im_ext, 'max')
        [~, locs] = findpeaks(lp_beats(:, poor_IM(beat)));
    else
        [~, locs] = findpeaks(-1*lp_beats(:, poor_IM(beat)));
    end
    
    
    if closest 
        
        % find the closest 
        [~, close_med] = min(abs(locs - IM_med));
        final_im(poor_IM(beat)) = locs(close_med);
        
    else
        
        % find the extrema within a range 
        if strcmp(im_ext, 'max')
            
            % find all peaks within a range 
            valid_ext = locs(abs(locs-IM_med) < f1);
            
            if ~isempty(valid_ext)
                % find the index assoicated with the largest peak
                final_im(poor_IM(beat)) = valid_ext(lp_beats(valid_ext, poor_IM(beat)) ...
                    == max(lp_beats(valid_ext, poor_IM(beat))));
                
            % otherwise just find the closest peak/valley
            else
                [~, close_med] = min(abs(locs - IM_med));
                final_im(poor_IM(beat)) = locs(close_med);
            end
            
        else
            
            % find all peaks within a range
            valid_ext = locs(abs(locs-IM_med) < f1);
            
            if ~isempty(valid_ext)
                % find the index associated with the lowest valley      
                final_im(poor_IM(beat)) = valid_ext(lp_beats(valid_ext, poor_IM(beat))...
                    == min(lp_beats(valid_ext, poor_IM(beat))));
                            % otherwise just find the closest peak
            else
                [~, close_med] = min(abs(locs - IM_med));
                final_im(poor_IM(beat)) = locs(close_med);
            end
            
        end
    end
    
end

if Plot
    figure;
    subplot(2, 1, 1);
    plot(final_im, 'LineWidth', 2); hold on; plot(im_locs); 
    title('IM Locations')
    xlabel('Beats'); ylabel('IM (ms)'); legend('Resampled Estimates', 'Initial Estimates')
    subplot(2, 1, 2);
    plot(final_ac, 'LineWidth', 2); hold on; plot(ac_locs);
    title('AC Locations')
    xlabel('Beats'); ylabel('AC (ms)'); legend('Resampled Estimates', 'Initial Estimates')
end

%clearvars -except ac_dir ac_ext  ac_locs ac_ccd alpha_range baseSCG1 cost ...
%    env_IM env_AC env_beats final_im final_ac fs hp_beats hr ...
%    im_dir im_ext im_locs im_ccd lp_beats limit ...
%    totRAC totRAO pig Plot rng_beat 


end
%%% Given an alpha, an enveloped high frequency beat, the low freqeuncy
%%% beat, and the location of the envelope peak (for IM/AC detection), 
%%% find the best AC or IM point to the left or right of alpha percent of this peak 
function [fid_point, bnds] = alpha_search(lp_beat, env_beat, peak_loc, alpha, fs, extrema, dir)

% search for the closest min or max in a smaller window alpha % around the
% peak of the envelope signal 
if strcmp(extrema, 'max'); findmax = 1; [~, locs] = findpeaks(lp_beat);
elseif strcmp(extrema, 'min'); findmax = 0; [~, locs] = findpeaks(-1 * lp_beat); end

if strcmp(dir, 'left'); left = 1;
elseif strcmp(dir, 'right'); left = 0; end

% start at the peak loc and find the the location of 
% just get an abitrary sized window around the peak say 150ms each side
bounds = [max(1, peak_loc-150*fs/1000), min(length(lp_beat), peak_loc+100*fs/1000)];

% now find the alpha% locations in this bound 
[~, low_bound] = min(abs(env_beat(bounds(1):peak_loc) - alpha*env_beat(peak_loc)));
[~, upp_bound] = min(abs(env_beat(peak_loc:bounds(2)) - alpha*env_beat(peak_loc)));

% for the lower bound get the value closest to the peak loc (end) 
low_alph_bound = low_bound(end) + bounds(1) - 1;
upp_alph_bound = upp_bound(1) + peak_loc - 1;

% find the potential locations 
if left
    pot_locs = locs(locs < low_alph_bound);
else
    pot_locs = locs(locs > upp_alph_bound);
end

% if none exist return nan 
if isempty(pot_locs)
    fid_point = NaN;
else
    if left
        [~, best_fid] = min(abs(pot_locs - low_alph_bound));
        fid_point = pot_locs(best_fid);
    else 
        [~, best_fid] = min(abs(pot_locs - upp_alph_bound));
        fid_point = pot_locs(best_fid);
    end
end

% save hte lower ad upper bounds 
bnds = [low_alph_bound, upp_alph_bound];
end

% Limitations: 
% 1. The AC detection range is dependent on HR and in this paper was found by
% solving a regression between HR and AC location. However, the parameters
% set are not typically transferrable 
% 2. The search range assumes the desired point is always to the left of
% the envelope peak.
% 3. The outlier removal step only gets rid of single point outliers. If
% the IM/AC peak suddenly swithces peaks for long periods of time the
% outlier step will not correct this
% 4. The median used for outliers is a global one which wouldn't work very
% well with dynamic hemodynamic changes (window argument) 




