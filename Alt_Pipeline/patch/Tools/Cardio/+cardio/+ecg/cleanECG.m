function [new_HR, clean_beats] = cleanECG(ECG,Fs,HRV_bound,varargin)
%--------------------------------------------------------------------------
% Created by: Jacob Kimball
% Date Created: May 15, 2019
% Inan Research Lab @ GaTech
%--------------------------------------------------------------------------
% This function takes a filtered ECG trace (no significant DC fluctuations)
% and returns vectors containing R peaks and the resultant HR. 
%
% Best use case: run this function, look for major issues in HR or deltaHR,
% fix them (manually in the ECG trace or change hyperparameters) and then
% run it again. I recommend plotting a mechanical/pressure wave on top of
% the ECG to verify productive heartbeats, especially in the case of PVC.
%
% Inputs: 
% ECG:      filtered ECG
% Fs:       sampling frequency
% HRV bound: expected beat to beat variability limit in HR. In a sedated
% animal, a good number is usually 2. 
%
% Optional inputs: 
% num_thresh:       referring to number of thresholds for automatedHR function and findpeaks. default 30
% peak_sep:         referring to peak separation for automatedHR function and findpeaks. Average is about
%                   850, but if your HR is extremely high or if you have an arrhythmia you want to catch
%                   you may need to drop it a little lower. 
% double:           flag indicating how far back you want to average HR to look for the next beat
% plotsoff:         flag indicating you don't want the plots to pop up
%
% NOTE: IN ORDER FOR THIS FUNCTION TO WORK YOU HAVE TO HAVE A
% DECENTLY CLEAN ECG SIGNAL. If you've filtered your ECG 
% and there's still a lot of baseline wander, consider 
% taking the derivative. This will give you the same 
% relative distances for all timing measures and works 
% well for windowing other signals. 
%
% This function assumes that the first beats are clean. (ie findpeaks makes
% the right decision for the first 3 peaks) 
%
%--------------------------------------------------------------------------
% Parse optional arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'num_thresh'); num_thresholds = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'peak_sep'); peak_separation = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'double'); double_back = true;
        elseif strcmp(varargin{arg}, 'plotsoff'); PPlot = false;
        end
    end
end

% Set defaults
if ~exist('num_thresholds', 'var'); num_thresholds = 30; end
if ~exist('peak_separation', 'var'); peak_separation = 850; end
if ~exist('double_back', 'var'); double_back = false; end
if ~exist('PPlot', 'var'); PPlot = true; end

% run automatedHR, get all the regular peak locations.
[tHR, HR, threshold] = ...
    cardio.ecg.automatedHR(ECG, 'numThresholds', num_thresholds, 'Fs', Fs, 'minDist', peak_separation);

mixed_beats = tHR*Fs;%int32(tHR*Fs); % sample locations of peaks

% run findpeaks with a much smaller interval but the same threshold to get
% all the peaks  
[~,Rpeakslocs] = findpeaks(ECG,'MinPeakHeight',threshold,...
    'MinPeakDistance',30);
% the idea is to find all of the major peaks in the signal (including
% the others previously found, with the exception of a few outliers)

% start with the first few good intervals, [manually verify that the first few beats are good] then look at the next.
%   if this next peak is going to change the HR by more than |HRV_bound|, look in
%   the allpeaks vector until you find a peak that fits in the right
%   window. Then continue.

% HR = beats/min = beats/samples * samples/sec * sec/min
%                = 1/(peak_b-peak_a) * 2000 * 60



% WHY DOES IT NEVER FIND THE FIRST PEAK? It has something to do with the
% original findpeaks function. It doesn't find the first one. 

new_HR = [Fs*60/double(mixed_beats(2)-mixed_beats(1)) Fs*60/double(mixed_beats(3)-mixed_beats(2))]; % first 2 HR (verified good)
clean_beats = [mixed_beats(1) mixed_beats(2) mixed_beats(3)]; % first three good beats
double_dip = [];
for i = 4:length(mixed_beats)-1 % test all of the mixed beats
    beat = i;
%     peak_a = clean_beats(i-1);
    peak_a = clean_beats(end); 
    peak_b = mixed_beats(i);
    test_HR = Fs*60*1/(peak_b-peak_a);
    HR_diff = test_HR - new_HR(end);
    change = HRV_bound;
    if(abs(HR_diff) > change)
        % resample within the designated window of opportunity
        [abs_diff,rb_ind] = min(abs(Rpeakslocs - peak_b)); % find the closest peak in the all_peaks vector
        % I might want to check if abs_diff isn't equal to 0.
        % I MIGHT ALSO WANT TO THROW A FLAG IF IT'S TOO BIG
        
        f_ind = rb_ind; % future index
        pb_check = Rpeakslocs(f_ind);
        future_peaks = [pb_check];
        
        
        % THERE ARE 3 OPTIONS HERE:
        %   1. sign = -1: interval is too long
        %   2. sign = 1: interval too short
        %   3. sign = 0: we missed a beat in there somewhere...maybe?
        
        if(sign(HR_diff) == -1) % interval is too long. look backwards
            %             too_long = 1;
            while((pb_check-peak_a)>peak_separation)
                f_ind = f_ind-1;
                pb_check = Rpeakslocs(f_ind);
                future_peaks = [future_peaks pb_check];
            end
        else
            if(sign(HR_diff) == 1)% interval is too short. Look forwards
                % how many potential peaks do I have to choose from?
                % how many samples is the 2bpm limit at this current HR?
                % samples/beat = fs_b*60/HR
                sample_limit = Fs*60/(new_HR(end)-change); % max allowable decrease in HR
                terminate = false;
                while(((pb_check - peak_a) < sample_limit) && ~terminate)
                    if((f_ind+1) > length(Rpeakslocs))
                        terminate = true;
                    end
                    f_ind = min(f_ind+1,length(Rpeakslocs));
                    pb_check = Rpeakslocs(f_ind);
                    future_peaks = [future_peaks pb_check];
                end
            end
        end
        
        if(new_HR(end) > 200) % this removes the possibility of trying to find HR = inf from landing on the same peak over and over
            HR_avg = double((new_HR(end-1))+new_HR(end-2))/2; % find the average HR over the last two good beats
        else
            if(double_back&&(i>10))
                HR_avg = double(new_HR(end)+new_HR(end-1)+new_HR(end-2)+new_HR(end-3))/4;
            else
                HR_avg = double(new_HR(end));
            end
        end
        
        [min_HRV,ind] = min(abs(HR_avg - Fs*60./(double(future_peaks) - double(peak_a)))); % find the min HRV, but take the average of the last two beats.
        peak_b = future_peaks(ind);
        test_HR = double(double(Fs*60)/(double(peak_b) - double(peak_a)));
    end
    
    if((peak_b-clean_beats(end))>500) % positive and HR < 240
        new_HR = [new_HR double(test_HR)]; % add the HR
        clean_beats = [clean_beats peak_b]; % add the beat
        % otherwise, skip and move to the next
    end
    progress = i/length(mixed_beats);
end

if(PPlot)
    figure();
    plot(ECG)
    hold on;
    plot(mixed_beats,ECG(int32(mixed_beats)),'x')
    plot(clean_beats,ECG(int32(clean_beats)),'o')
    title('difference in beat location'); xlabel('sample');
    legend('ECG','mixed beats','clean beats')
    
    figure();
    plot(HR,'r')
    hold on;
    plot(new_HR,'k')
    ylabel('HR'); xlabel('beat');
    legend('old HR','new HR')
    
    figure();
    plot(diff(HR),'r')
    hold on;
    plot(diff(new_HR),'k')
    ylabel('delta HR'); xlabel('beat');
    legend('old HR','new HR')
end
end
