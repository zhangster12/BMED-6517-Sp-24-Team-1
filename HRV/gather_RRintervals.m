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

    load([feature_data_path feature_file.name], 'RRintervals');

    timing = timing - timing(1);

    timing_part1 = timing(1:4);
    timing_part2 = timing(5:8);

    timing = [timing_part1; timing(4); timing(5); timing_part2];

    feature_folder = ['RRintervals/' sub_id];

    if ~exist(feature_folder, 'dir')
        mkdir(feature_folder);
    end

    for j = 1:length(timing)/2

        start_time = timing(2*j-1);
        end_time = timing(2*j);

        RRseg = RRintervals(find(start_time <= RRintervals(:,1) & RRintervals(:,1) <= end_time), :);
        
        writematrix(RRseg, [feature_folder '/stim' int2str(j-1) '_RRintervals.csv']);

    end

end