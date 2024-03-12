function [cleanNN, cleantNN] = RRIntervalPreprocess_skb(rr,fs)
%
    t_rr = cumsum(rr);


idx_remove = find(diff(t_rr) < 1/fs); 
% could make this the refractory period - and have the variable in the settings
% document
rr(idx_remove+1) = [];
t_rr(idx_remove) = [];
% annotations(idx_remove) = [];
clear idx_remove;
    

    % HRV normal to normal (NN) Settings
%     HRVparams.preprocess.per_limit = 0.2;  % Max. percent change in neighboring NN intervals
%     HRVparams.preprocess.gaplimit = 2;   % Longest RR interval allowed
%     HRVparams.preprocess.lowerphysiolim = 60/180;   % minimum RR interval length
%     HRVparams.preprocess.upperphysiolim = 60/30;    % maximum RR interval length


% 4. Remove Large RR intervals Caused by Gaps
% These are not counted towards the total signal removed
idx_remove = find(rr >=2);
rr(idx_remove) = [];
t_rr(idx_remove) = [];
% annotations(idx_remove) = [];
clear idx_remove;


% 6. Find RR Over Given Percentage Change 
% perLimit = 0.2;
% idxRRtoBeRemoved = FindSpikesInRR(rr, perLimit); 

%% call McNames impulse

% call impulse 
    
IHR=1./rr;
C_IHR= mcnames_impulse2(IHR,4,2); %% impulse filtering
rr=1./C_IHR; %% CRR1= corrected RR1


% idx_outliers = find(outliers == 1);
% 
% % Keep count of outliers 
% numOutliers = length(idx_outliers);
% 
% rr_original = rr;
% rr(idx_outliers) = NaN;


NN_Outliers = interp1(t_rr,rr,t_rr,'spline','extrap');
t_Outliers = t_rr;

% 

% 10. Identify Non-physiologic Beats
toohigh = NN_Outliers > 2;    % equivalent to RR = 2
toolow = NN_Outliers < 0.375;     % equivalent to RR = .375

idx_toolow = find(toolow == 1);
NN_NonPhysBeats = NN_Outliers;
NN_NonPhysBeats(idx_toolow) = NaN;
% numOutliers = numOutliers + length(idx_toolow);

% 
% idx_remove = find(NN_Outliers >=2);
% NN_Outliers(idx_remove) = NaN;
% t_rr(idx_remove) = [];
% % annotations(idx_remove) = [];
% clear idx_remove;



% switch HRVparams.preprocess.method_unphysio
%     case 'cub'
        NN_NonPhysBeats = interp1(t_Outliers,NN_NonPhysBeats,t_Outliers,'spline','extrap');
        t_NonPhysBeats = t_Outliers;
%         flagged_beats = logical(outliers(:) + toohigh(:)+ toolow(:));

% end

% if figures
%     hold on;
%     plot(t_NonPhysBeats,NN_NonPhysBeats+.01);
%     hold on; plot(t_NonPhysBeats,toolow,'o')
%     legend('raw','interp1(after outliers removed)',...
%         'interp2(after too low)','toolow')
% end
% 


% 11. Interpolate Through Beats that are Too Fast
toohigh = NN_NonPhysBeats > 2;    % equivalent to RR = 2

idx_outliers_2ndPass = find(logical(toohigh(:)) ~= 0);
NN_TooFastBeats = NN_NonPhysBeats;
NN_TooFastBeats(idx_outliers_2ndPass) = NaN;
% numOutliers = numOutliers + length(idx_outliers_2ndPass);
% if strcmp(HRVparams.preprocess.method_unphysio,'rem')
%     flagged_beats = numOutliers;
% end

% switch HRVparams.preprocess.method_outliers
%     case 'cub'            
        NN_TooFastBeats = interp1(t_NonPhysBeats,NN_TooFastBeats,t_NonPhysBeats,'spline','extrap');
        t_TooFasyBeats = t_NonPhysBeats;
%    end



% 12. Remove erroneous data at the end of a record 
%       (i.e. a un-physiologic point caused by removing data at the end of
%       a record)

while NN_TooFastBeats(end) > 2	% equivalent to RR = 2
    NN_TooFastBeats(end) = [];
    t_TooFasyBeats(end) = [];
end


cleanNN = NN_TooFastBeats;
cleantNN = t_TooFasyBeats;
end


