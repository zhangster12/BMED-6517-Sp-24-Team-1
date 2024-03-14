% Dependencies: 
%   cardio library https://github.gatech.edu/IRL/Cardio.git
%
% Data sets: 
%   hypovolemia dataset (assumes beat segmented SCG beats)
%
% Reproducing:
%   J. Zia, J. Kimball, C. Rozell and O. T. Inan, "Harnessing the Manifold
% Structure of Cardiomechanical Signals for Physiological Monitoring During 
% Hemorrhage," in IEEE Transactions on Biomedical Engineering, vol. 68, no. 
% 6, pp. 1759-1767, June 2021, doi: 10.1109/TBME.2020.3014040.
%
% Notes: 
%   Add path to dataset (processedData_x.m) 
%   10/3/2022 Only the PCA steps of generating the initial manifold is done
%   in this script. Tracking PEP using the manifold will be added later



% holds scg beats, ao values, and sqi values 
pig_scg = cell(6, 1);
pig_ao = cell(6, 1);
pig_sqi = cell(6, 1);

disp('Loading in Pig Data')
for pig_num = 1:6
    
    load(['processedData_' ,num2str(pig_num), '.mat'])
    Fs = oink.Fs;
    pig_scg{pig_num} = oink.dataset{1}.allAbsolute.sternumSCG_z;
    pig_ao{pig_num} = oink.dataset{1}.allAbsolute.scgRAO;
    clear oink
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start SQI and remove beats 

% first normalize beats (assumed beats are in columns)
pig_scg = cellfun(@(x) normalize(x, 1), pig_scg, 'UniformOutput', false);

% generate templates using the first 100 beats of each pig 
pig_templates = cellfun(@(x) mean(x(:, 1:100), 2), pig_scg, 'UniformOutput', false);

% 
disp('Calculating SQI')
for pig_num = 1:6
    % 1. for better performance remove ('type', Index.DTW) to perform DTFM.
    % DTFM is just more computationally heavy so it might take a while and
    % the paper only uses DTW
    % 2. can also consider using an average of the sqi obtained from using
    % the held out pigs' tempaltes
    % 3. can also consider a moving average (30 beats) to continuously generate a
    % template tN for beat N to get an SQI at each time point
    pig_sqi{pig_num} = cardio.sqi.filter(pig_scg{pig_num}, Fs, inf, ...
            'lambda', 25, 'template', pig_templates{pig_num}, 'tolType', 'percentage', 'type', Index.DTW);
end 

% transpose so beats are on rows     
pig_scg = cellfun(@(x) x', pig_scg, 'UniformOutput', false);

% set a sqi percentile
sqi_thresh = 0.8;

% initialize accepted and rejected beats for each
accepted_beats = cell(6, 1);
rejected_beats = cell(6, 1);

for pig_num = 1:6
    
    % get the sqi values for pig i
    sqi_arr = pig_sqi{pig_num};
    
    % get the sqi values in the sqi percentile and reject those not in that
    % percentile
    [sorted_sqi, idx] = sort(sqi_arr, 'descend');
    accepted_beats_i = idx(1:floor(length(sorted_sqi)*sqi_thresh));
    rejected_beats_i = idx(floor(length(sorted_sqi)*sqi_thresh)+1:end);
    accepted_beats{pig_num} = sort(accepted_beats_i, 'ascend');
    rejected_beats{pig_num} = sort(rejected_beats_i, 'ascend');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting 3 things 
% 1. Representation of PCA and SQI put together 
% 2. Representation of PCA applied to invidual pigs colorcoded by PEP
% 3. Representation fo PCA applied to all pigs colorcoded by Pig

% plotting parameters 
marker_size = 3;
colors = [0, 0.4470, 0.7410;  % Pig 1 color 
    0.8500, 0.3250, 0.0980;   % Pig 2 color
    0.9290, 0.6940, 0.1250;   % Pig 3 Color
    0.4940 0.1840 0.5560;     % Pig 4 Color
    0.4660 0.6740 0.1880;     % Pig 5 Color
    0.6350 0.0780 0.1840];    % Pig 6 Color


%%%%%%%%%%%% Plot 1 
pig_num = 6;
scg_arr = pig_scg{pig_num};

% PCA on original beats 
[coeff, scores] = pca(scg_arr);

% project data onto PCs that you found 
proj_accepted = (scg_arr(accepted_beats{pig_num}, :) - mean(scg_arr(accepted_beats{pig_num}, :), 1)) * coeff(:, 1:3);
proj_rejected = (scg_arr(rejected_beats{pig_num}, :) - mean(scg_arr(rejected_beats{pig_num}, :), 1)) * coeff(:, 1:3);

% PCA on cleaner beats
[cln_coeff, cln_scores] = pca(scg_arr(accepted_beats{pig_num}, :));


disp('Plotting PCA and SQI Steps')
figure; 
ax(1) = subplot(1, 3, 1); 
scatter3(scores(:, 1), scores(:, 2), scores(:, 3), marker_size, 'filled', 'MarkerFaceColor', 'black')
xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
title('PCA Plot before SQI');

ax(2) = subplot(1, 3, 2);
s1 = scatter3(proj_accepted(:, 1), proj_accepted(:, 2), proj_accepted(:, 3), ...
    marker_size, 'filled', 'MarkerFaceColor', [0.7, 0.7, 0.7], 'DisplayName', 'High SQI'); hold on;
s2 = scatter3(proj_rejected(:, 1), proj_rejected(:, 2), proj_rejected(:, 3), ...
    marker_size, 'filled', 'MarkerFaceColor', 'red', 'DisplayName', 'Low SQI');
legend([s1, s2], 'Location', 'best')
title('PCA Plot Marking Low SQI Beats')

ax(3) = subplot(1, 3, 3);
scatter3(cln_scores(:, 1), cln_scores(:, 2), cln_scores(:, 3), marker_size, 'filled', 'MarkerFaceColor', 'black')
title('PCA Plot after Low SQI Beat Removal')
linkprop(ax, 'CameraPosition');
campos([-200   80  300])
set(gcf, 'Position', [400         225        1000        700])


%%%%%%%%%%%% Plot 2
figure; 
disp('Plotting Individual PCA Plots')
for pig_num = 1:6
    
    accepted_beats_i = accepted_beats{pig_num};
    scg_beats_i = pig_scg{pig_num};
    ao_values = pig_ao{pig_num}(accepted_beats_i);
    [coeff, scores] = pca(scg_beats_i(accepted_beats_i, :));
    
    % set a gradient for each color (gray - low value to the pig's color - high value)
    gray = 0.9;
    color_grad = [linspace(gray, colors(pig_num, 1), 100)', ...
        linspace(gray,colors(pig_num, 2), 100)', ...
        linspace(gray,colors(pig_num, 3), 100)'];
    
    bx(pig_num) = subplot(2, 3, pig_num);
    scatter3(scores(:, 1), scores(:, 2), scores(:, 3), marker_size, ao_values, 'filled')
    colormap(bx(pig_num), color_grad);
    c = colorbar('northoutside');
    if pig_num == 1
        ylabel(c, 'PEP (ms)')
    end
    xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
    title(['Pig ', num2str(pig_num)])
    sgtitle('Individual PCA')
    
    % set the camera angle so it's easier to see the curvature
    if pig_num == 1
        campos([-110  150 320])
    elseif pig_num == 2
        campos([40  220 300])
    elseif pig_num == 3
        campos([120  300 280])
    elseif pig_num == 4
        campos([140  -240 300])
    elseif pig_num == 5
        campos([80  -400 180])
    elseif pig_num == 6
        campos([-30  200 320])
    end
end
set(gcf, 'Position', [400         225        1000        700])


%%%%%%%%%%%%% Plot 3
disp('Plotting aggregate PCA Plot')
all_beats = [];
all_colors = [];

% concatenate all beats and save which pig each beat came from
for pig_num = 1:6
   all_beats = [all_beats;pig_scg{pig_num}(accepted_beats{pig_num}, :)];
   all_colors = [all_colors;pig_num*ones(length(accepted_beats{pig_num}), 1)];
end

% apply pca on all beats 
[coeff, scores] = pca(all_beats);


figure; 
for pig_num = 1:6
    scatter3(scores(all_colors == pig_num, 1), scores(all_colors == pig_num, 2), ...
        scores(all_colors == pig_num, 3), marker_size, colors(pig_num, :), 'filled', ...
        'DisplayName', ['Pig ', num2str(pig_num)]); hold on;
end
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
legend()
title('PCA all Pigs')
campos([120, -80, 360])
set(gcf, 'Position', [400         225        1000        700])