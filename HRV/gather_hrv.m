close all;
clear all;
clc;

timing_data_path = '../Timing/';
feature_data_path = '../Seis_Process_R2/';


timing_data_folder = dir([timing_data_path '*.mat']);

for i = 1:length(timing_data_folder)

    sub_id = timing_data_folder(i).name(1:7);
    timing = load([timing_data_path sub_id '.mat']).timing;

    feature_file = dir([feature_data_path '*' sub_id(4:end) '*.mat']);

    load([feature_data_path feature_file.name]);

    timing = timing - timing(1);

    timing_part1 = timing(1:4);
    timing_part2 = timing(5:8);

    timing = [timing_part1; timing(4); timing(5); timing_part2];

    feature_folder = ['Consolidated_Features/' sub_id];

    if ~exist(feature_folder, 'dir')
        mkdir(feature_folder);
    end

    HR = physFeatures_beat.HR;
    PPGamp = physFeatures_beat.PPGamp;
    PAT = physFeatures_beat.PAT;
    PEP = physFeatures_beat.PEP;
    PTTrecip = physFeatures_beat.PTT_recip;

    for j = 1:length(timing)/2

        start_time = timing(2*j-1);
        end_time = timing(2*j);

        HR_seg = HR(find(start_time <= HR(:,1) & HR(:,1) <= end_time), :);
        PPGamp_seg = PPGamp(find(start_time <= PPGamp(:,1) & PPGamp(:,1) <= end_time), :);
        PAT_seg = PAT(find(start_time <= PAT(:,1) & PAT(:,1) <= end_time), :);
        PEP_seg = PEP(find(start_time <= PEP(:,1) & PEP(:,1) <= end_time), :);
        PTTrecip_seg = PTTrecip(find(start_time <= PTTrecip(:,1) & PTTrecip(:,1) <= end_time), :);

        writematrix(HR_seg, [feature_folder '/stim' int2str(j-1) '_HR.csv']);
        writematrix(PPGamp_seg, [feature_folder '/stim' int2str(j-1) '_PPGamp.csv']);
        writematrix(PAT_seg, [feature_folder '/stim' int2str(j-1) '_PAT.csv']);
        writematrix(PEP_seg, [feature_folder '/stim' int2str(j-1) '_PEP.csv']);
        writematrix(PTTrecip_seg, [feature_folder '/stim' int2str(j-1) '_PTTrecip.csv']);

    end

end