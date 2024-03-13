close all;
clear all;
clc;

% Path to folder storing RR intervals
RRinterval_data_path = 'RRintervals/';
RRinterval_data_folder = dir([RRinterval_data_path 'sub*']);

Fs = 500; % ECG sampling frequency

for i = 1:length(RRinterval_data_folder)

    sub_id = RRinterval_data_folder(i).name(1:7);

    HRV_filename = ['HRV_feat_all/' sub_id '_HRV_feat.csv'];

    % Store each HRV feature in a table
    HRV_feat = table;

    for j = 0:4

        temp = table;

        % Store the stimulus #
        temp.stim = j;

        % Load the RR intervals
        RRinterval_filename = [RRinterval_data_path sub_id '/stim' int2str(j) ...
            '_RRintervals.csv'];
        RRintervals = readmatrix(RRinterval_filename);

        % Isolate the last two minutes of RR interval data for each stimulus
        RRinterval_t_end = RRintervals(end, 1);
        RRinterval_t_start = RRinterval_t_end - 120;
        RRintervals = RRintervals(find( RRinterval_t_start <= RRintervals(:,1) & ... 
            RRintervals(:,1) <= RRinterval_t_end), 2);

        HRV_feat_all = RRtoHRVDemo(RRintervals, Fs);

        % Here, I will consolidate HRV_feat_all into a set of features that are
        % relevant for stress detection
        
        % Time-Domain HRV 
        temp.SDNN = HRV_feat_all(9);
        temp.RMMSD = HRV_feat_all(10);
        temp.pnn50 = HRV_feat_all(11);
        temp.nn50 = HRV_feat_all(12);
        
        % Nonlinear HRV
        temp.SD1 = HRV_feat_all(13);
        temp.SD2 = HRV_feat_all(14);
        temp.SDratio = HRV_feat_all(15);
        
        % Freq-Domain HRV
        temp.VLF = HRV_feat_all(20);
        temp.LF = HRV_feat_all(21);
        temp.HF = HRV_feat_all(22);
        temp.LFHF = HRV_feat_all(23);
        temp.LFn = HRV_feat_all(24);
        temp.HFn = HRV_feat_all(25);
        temp.POW = HRV_feat_all(26);
        
        HRV_feat = [HRV_feat; temp];
    end

    writetable(HRV_feat, HRV_filename);

end