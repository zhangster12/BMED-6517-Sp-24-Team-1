function obj = extractData(obj, varargin)

% -------------------------------------------------------------------------
% Extracting data from hypovolemia study dataset.
%
% ARGUMENTS (REQ'D)
% - obj         Oink        Hypovolemia dataset object
%
% ARGUMENTS (OPT'L)
% - filter      Bool        Filter data? (default true)
% - separate    Bool        Beat-separate data? (default true)
% - sigLen      Int         Fixed signal length after R-peak
% - levels      [Level]     Blood volume levels to extract
% - path        String      Path to dataset (if not on IRL server)
% - include     FLAG        Include data during bleeding for absolute hypovolemia
% - notruncation    FLAG    Do not truncate the beats (nan-pad the ends)
% -------------------------------------------------------------------------

% Printing progress
disp("Extracting Data:")

% Set defaults for optional arguments
obj.filtered = true; obj.beatSeparated = true; obj.levels = Level.allAbsolute; sigLen = 1000;
obj.path = "/media/Data/Hypovolemia_Study/"; include = false; truncation = true;

% Parse optional arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'filter'); obj.filtered = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'separate'); obj.beatSeparated = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'levels'); obj.levels = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'path'); obj.path = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'sigLen'); sigLen = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'include'); include = true;
        elseif strcmp(varargin{arg}, 'notruncation'); truncation = false;
        end
    end
end

% Add dataset to path
addpath(genpath(obj.path));

% Initialize counter for dataset index
counter = 1;

% Process the data for each of the selected subjects
for i = 1:length(obj.subjects)
    
    % Print progress
    disp("-> Processing Subject " + string(i) + " of " + string(length(obj.subjects)))
    
    % Switch over the current subect
    switch obj.subjects(i)
        
        case 1
            
            % Create global variables
            ECG_b = []; ECG_t = []; AX_1 = []; AY_1 = []; AZ_1 = []; AX_2 = []; 
            AY_2 = []; AZ_2 = []; SPO2_1 = []; SPO2_2 = []; PPG_1 = []; PPG_2 = [];
            P1 = []; P2 = []; P3 = []; P4 = []; HR_b = []; HR_t = []; PCG = [];
            beats_b = []; beats_t = [];
            
            % Process data for Subject 1
            subject1();
            
        case 2
            
            % Create global variables
            ECG_b = []; ECG_t = []; AX_1 = []; AY_1 = []; AZ_1 = []; AX_2 = []; 
            AY_2 = []; AZ_2 = []; SPO2_1 = []; SPO2_2 = []; PPG_1 = []; PPG_2 = [];
            P1 = []; P2 = []; P3 = []; P4 = []; HR_b = []; HR_t = []; PCG = [];
            beats_b = []; beats_t = [];
            
            % Process data for Subject 2
            subject2();
            
        case 3
            
            % Create global variables
            ECG_b = []; ECG_t = []; AX_1 = []; AY_1 = []; AZ_1 = []; AX_2 = []; 
            AY_2 = []; AZ_2 = []; SPO2_1 = []; SPO2_2 = []; PPG_1 = []; PPG_2 = [];
            P1 = []; P2 = []; P3 = []; P4 = []; HR_b = []; HR_t = []; PCG = [];
            beats_b = []; beats_t = [];
            
            % Process data for Subject 3
            subject3();
            
        case 4
            
            % Create global variables
            ECG_b = []; ECG_t = []; AX_1 = []; AY_1 = []; AZ_1 = []; AX_2 = []; 
            AY_2 = []; AZ_2 = []; SPO2_1 = []; SPO2_2 = []; PPG_1 = []; PPG_2 = [];
            P1 = []; P2 = []; P3 = []; P4 = []; HR_b = []; HR_t = []; PCG = [];
            beats_b = []; beats_t = [];
            
            % Process data for Subject 4
            subject4();
            
        case 5
            
            % Create global variables
            ECG_b = []; ECG_t = []; AX_1 = []; AY_1 = []; AZ_1 = []; AX_2 = []; 
            AY_2 = []; AZ_2 = []; SPO2_1 = []; SPO2_2 = []; PPG_1 = []; PPG_2 = [];
            P1 = []; P2 = []; P3 = []; P4 = []; HR_b = []; HR_t = []; PCG = [];
            beats_b = []; beats_t = [];
            
            % Process data for Subject 5
            subject5();
            
        case 6
            
            % Create global variables
            ECG_b = []; ECG_t = []; AX_1 = []; AY_1 = []; AZ_1 = []; AX_2 = []; 
            AY_2 = []; AZ_2 = []; SPO2_1 = []; SPO2_2 = []; PPG_1 = []; PPG_2 = [];
            P1 = []; P2 = []; P3 = []; P4 = []; HR_b = []; HR_t = []; PCG = [];
            beats_b = []; beats_t = [];
            
            % Process data for Subject 6
            subject6();
            
        otherwise
            disp("Error in Oink.extractData(): Subject" + string(obj.subjects(i)) + " not found");
            return
    end
    
    % Increment counter
    counter = counter + 1;
    
end

%% ------------------------------------------------------------------------
% Sub-Function for Subject 1
% -------------------------------------------------------------------------

    function subject1()
        
        % Returns dataset for Subject 1 as a struct:
        % obj.dataset{1}.(interval).ecg_biopac
        % obj.dataset{1}.(interval).ecg_t3
        % obj.dataset{1}.(interval).sternumSCG_x
        % obj.dataset{1}.(interval).sternumSCG_y
        % obj.dataset{1}.(interval).sternumSCG_z
        % obj.dataset{1}.(interval).sternumSCG_z_back
        % obj.dataset{1}.(interval).apexSCG_x
        % obj.dataset{1}.(interval).apexSCG_y
        % obj.dataset{1}.(interval).apexSCG_z
        % obj.dataset{1}.(interval).apexSCG_z_back
        % obj.dataset{1}.(interval).apexSPO2
        % obj.dataset{1}.(interval).femoralSPO2
        % obj.dataset{1}.(interval).apexPPG
        % obj.dataset{1}.(interval).femoralPPG
        % obj.dataset{1}.(interval).aorticPressure
        % obj.dataset{1}.(interval).aorticPressure_back
        % obj.dataset{1}.(interval).femoralPressure
        % obj.dataset{1}.(interval).wedgePressure
        % obj.dataset{1}.(interval).rightAtrialPressure
        
        % -----------------------------------------------------------------
        % Load Dataset and Setup Workspace
        % -----------------------------------------------------------------
        
        % Import data from Biopac
        disp("--> Step 1 of 7: Importing Biopac Data")
        biopac = importdata('Subject1/biopac.txt');
        
        ECG_b = biopac(:, 1);    % Biopac ECG
        SPO2_1 = biopac(:, 2);   % Apex SPO2
        SPO2_2 = biopac(:, 3);   % Femoral SPO2
        PPG_1 = biopac(:, 4);    % Apex PPG
        PPG_2 = biopac(:, 5);    % Femoral PPG
        AX_1 = biopac(:, 6);     % Sternum SCG (x-acceleration)
        AY_1 = biopac(:, 7);     % Sternum SCG (y-acceleration)
        AZ_1 = biopac(:, 8);     % Sternum SCG (z-acceleration)
        AX_2 = biopac(:, 9);     % Apex SCG (x-acceleration)
        AY_2 = biopac(:, 10);    % Apex SCG (y-acceleration)
        AZ_2 = biopac(:, 11);    % Apex SCG (z-acceleration)
        
        % Clear the biopac workspace variable to conserve memory
        clear biopac
        
        % Set the sampling frequency, period, and time vector for Biopac
        % fs = 2000; ts = 1/fs; t_b = 0:ts:ts*(length(AX_1)-1);

        % Import data from T3
        disp("--> Step 2 of 7: Importing T3 Data")
        T3 = importdata('Subject1/t3.txt'); T3 = T3.data;

        % t3 = T3(:,1);     % Time vector
        P1 = T3(:, 2);      % Aortic arch pressure
        P2 = T3(:, 3);      % Femoral artery pressure
        P3 = T3(:, 4);      % Wedge pressure
        P4 = T3(:, 5);      % Right atrial pressure
        ECG_t = T3(:, 6);    % T3 ECG (1)
        
        % Clear the T3 workspace variable to conserve memory
        clear T3
        
        % -----------------------------------------------------------------
        % Remove NaNs and Resample T3 Data
        % -----------------------------------------------------------------
        
        disp("--> Step 3 of 7: Cleaning Data")
        
        % The T3 data sometimes has NaNs. These come from various forms of noise.
        % However, they seem to happen in different locations in each signal and
        % not every signal has NaNs, which means we replace the NaNs rather than
        % remove them. 

        % Which signals have NaNs?
        %   1. ECG1 - Likely electrical noise
        %   2. ECG2 - Likely electrical noise or due to contact
        %   3. P2 - Likely due to contact
        
        % Remove NaNs
        for j = 1:length(ECG_t)
            if isnan(ECG_t(j)); ECG_t(j) = ECG_t(j-1); end
            if isnan(P1(j)); P1(j) = P1(j-1); end
            if isnan(P2(j)); P2(j) = P2(j-1); end
            if isnan(P3(j)); P3(j) = P3(j-1); end
            if isnan(P4(j)); P4(j) = P4(j-1); end
        end
        
        % Re-sample signals from 1kHz to 2kHz
        % t3 = resample(t3, 2, 1);      % Time vector
        ECG_t = resample(ECG_t, 2, 1);  % T3 ECG (1)
        P1 = resample(P1, 2, 1);      % Aortic arch pressure
        P2 = resample(P2, 2, 1);      % Femoral artery pressure
        P3 = resample(P3, 2, 1);      % Wedge pressure
        P4 = resample(P4, 2, 1);      % Right atrial pressure
        
        % -----------------------------------------------------------------
        % Filter ECG Data
        % -----------------------------------------------------------------
        
        disp("--> Step 4 of 7: Filtering Data")
        
        if obj.filtered
            
            % ECG
            % Standard Kaiser window BP filter |5/10-25\40|
            Hd = cardio.general.createBPF(5, 40, obj.Fs, 'fpass1', 10, 'fpass2', 25, 'kaiser', 'order', 20);
            ECG_b = filtfilt(Hd.Numerator, 1, ECG_b); ECG_t = filtfilt(Hd.Numerator, 1, ECG_t);
            
            % SCG
            % Standard Kaiser window BP filter |0.5/1-20\40|
            Hd = cardio.general.createBPF(0.5, 40, obj.Fs, 'fpass1', 1, 'fpass2', 20, 'kaiser', 'order', 20);
            AX_1 = filtfilt(Hd.Numerator, 1, AX_1); AX_2 = filtfilt(Hd.Numerator, 1, AX_2);
            AY_1 = filtfilt(Hd.Numerator, 1, AY_1); AY_2 = filtfilt(Hd.Numerator, 1, AY_2);
            AZ_1 = filtfilt(Hd.Numerator, 1, AZ_1); AZ_2 = filtfilt(Hd.Numerator, 1, AZ_2);
            
            % PPG / Arterial Pressures
            % Standard Kaiser window BP filter |0.5/1-9\10|
            Hd = cardio.general.createBPF(0.5, 10, obj.Fs, 'fpass1', 1, 'fpass2', 9, 'kaiser', 'order', 20);
            PPG_1 = filtfilt(Hd.Numerator, 1, PPG_1); PPG_2 = filtfilt(Hd.Numerator, 1, PPG_2);
            P1 = filtfilt(Hd.Numerator, 1, P1); P2 = filtfilt(Hd.Numerator, 1, P2);
            P3 = filtfilt(Hd.Numerator, 1, P3); P4 = filtfilt(Hd.Numerator, 1, P4);
            
            % PCG
            % Standard Kaiser window HP filter (50Hz)
            Hd = cardio.general.createHPF('Fs', obj.Fs, 'Fc', 50, 'order', 50);
            PCG = filtfilt(Hd.Numerator, 1, AZ_1);
            
        end
        
        % -----------------------------------------------------------------
        % Process ECG (BIOPAC)
        % -----------------------------------------------------------------
        
        disp("--> Step 5 of 7: Processing ECG Data (Biopac)")
        
        % 1. Check polarity. (Are the R-peaks upside down?) 
        % 2. Take the difference, check its shape.
        % 3. Append a 0 as the last element so it is the same length as the
        % other signals. 

        % Differentiate ECG signal
        diff1_ecgb = diff(-ECG_b);

        % Append a 0 so it's the same length as everything else
        ECG_bio = [diff1_ecgb; zeros(1)]; clear diff1_ecgb; 
        
        % Fix outliers so ECG_bio can be processed by cardio.ecg.automatedHR
        ECG_bio_trim = zeros(size(ECG_bio)); limit = 0.015;
        for j = 1:length(ECG_bio)
           if(abs(ECG_bio(j)) > limit)
               % Keep the relative distance so you don't lose the peaks
               ECG_bio_trim(j) = (0.01*(abs(ECG_bio(j))-limit)+limit)*sign(ECG_bio(j)); 
           else; ECG_bio_trim(j) = ECG_bio(j); 
           end
        end

        % Clear extra workspace variable
        clear ECG_bio
        
        % Set parameters for further further cleaning and beat-separation
        % HR is going to be below 140 BPM.
        peak_separation = 845; HRV_bound = 2;
        ECG_bt_2 = ECG_bio_trim; clear ECG_bio_trim;

        % Fix the ~4 peaks that aren't getting picked up properly
        for j = 1:length(ECG_bt_2)
            if((j > 2004829) && (j < 2004836))
                ECG_bt_2(j) = ECG_bt_2(j)*5;
            elseif (j > 9581944) && (j < 9581956)
                ECG_bt_2(j) = ECG_bt_2(j) + 0.021;
            elseif (j > 10904774) && (j < 10904796)
                ECG_bt_2(j) = ECG_bt_2(j) + 0.029;
            elseif (j > 26739383) && (j < 26739399)
                ECG_bt_2(j) = ECG_bt_2(j) + 0.014;
            end
        end

        % Extract heart rate and beat-separated ECG
        [HR_b, beats_b] = cardio.ecg.cleanECG(ECG_bt_2, obj.Fs, HRV_bound, 'peak_sep', peak_separation, 'double', 'plotsoff');
        
        % -----------------------------------------------------------------
        % Process ECG (T3)
        % -----------------------------------------------------------------
        
        disp("--> Step 6 of 7: Processing ECG Data (T3)")
        
        % Cut out the death sequence that messes everything up
        ECG_t3_trim = ECG_t(1:36397700);
        for j = 1:length(ECG_t3_trim)
            if(ECG_t3_trim(j) > 0.3); ECG_t3_trim(j) = 0.01*ECG_t3_trim(j)+0.3; end
        end

        % Rename workspace variable
        ECG_tt_2 = ECG_t3_trim; clear ECG_t3_trim;

        % Fix the ~5 peaks that aren't getting picked up properly
        for j = 1:length(ECG_tt_2)
            if j > 5320023 &&  j < 5320033
                ECG_tt_2(j) = ECG_tt_2(j) + 0.25;	% Large T wave followed by really low QRS
            elseif j > 5374801 && j < 5374813
                ECG_tt_2(j) = ECG_tt_2(j) + 0.18;	% Large T wave followed by really low QRS
            elseif j > 24232261 && j < 24232275
                ECG_tt_2(j) = ECG_tt_2(j) + 0.3;	% Just noise in the ECG
            elseif j > 25633162 && j < 25633170
                ECG_tt_2(j) = ECG_tt_2(j) + 2.25;
            elseif j > 31026851 && j < 31026863
                ECG_tt_2(j) = ECG_tt_2(j) + 0.55;
            end
        end

        % Extract heart rate and beat-separated ECG
        [HR_t, beats_t] = cardio.ecg.cleanECG(ECG_tt_2, obj.Fs, HRV_bound, 'peak_sep', peak_separation, 'double', 'plotsoff');

        % Because you cut out more noise at the end of T3, shorten the bio
        % heartbeats off the end so you're working with the same number of beats: 
        beats_b = double(beats_b(1:length(beats_t))); HR_b = double(HR_b(1:length(HR_t))); 
        beats_t = double(beats_t); HR_t = double(HR_t);
        
        % obj.dataset{counter}.beats_biopac = beats_b; obj.dataset{counter}.beats_t3 = beats_t;
        % obj.dataset{counter}.HR_biopac = HR_b; obj.dataset{counter}.HR_t3 = HR_t;
        
        % -----------------------------------------------------------------
        % Append Dataset by Label
        % -----------------------------------------------------------------
        
        disp("--> Step 7 of 7: Appending Dataset by Labels")
        
        % Call the trimToIdx function
        trimToIdx_s1(obj.levels, include);
        
        % For each level this subject has, find the beat indices
        allLevelsToBeats(1, include);
        
    end

%% ------------------------------------------------------------------------
% Sub-Function for Subject 2
% -------------------------------------------------------------------------

    function subject2()
        
        % NOTE: Subject 2 does not have interval Level.all
        % NOTE: The first element in any cell array is from the RELATIVE
        % trial while the second element is from the ABSOLUTE trial.
    
        % Returns dataset for Subject 2 as a struct:
        % obj.dataset{2}.(interval).ecg_biopac
        % obj.dataset{2}.(interval).ecg_t3
        % obj.dataset{2}.(interval).sternumSCG_x
        % obj.dataset{2}.(interval).sternumSCG_y
        % obj.dataset{2}.(interval).sternumSCG_z
        % obj.dataset{2}.(interval).sternumSCG_z_back
        % obj.dataset{2}.(interval).apexSCG_x
        % obj.dataset{2}.(interval).apexSCG_y
        % obj.dataset{2}.(interval).apexSCG_z
        % obj.dataset{2}.(interval).apexSCG_z_back
        % obj.dataset{2}.(interval).apexSPO2
        % obj.dataset{2}.(interval).femoralSPO2
        % obj.dataset{2}.(interval).apexPPG
        % obj.dataset{2}.(interval).femoralPPG
        % obj.dataset{2}.(interval).aorticPressure
        % obj.dataset{2}.(interval).femoralPressure
        % obj.dataset{2}.(interval).wedgePressure
        % obj.dataset{2}.(interval).rightAtrialPressure
        
        % -----------------------------------------------------------------
        % Load Dataset and Setup Workspace
        % -----------------------------------------------------------------

        % Import data from T3
        disp("--> Step 1 of 7: Importing T3 Data")
        T3 = importdata('Subject2/t3.txt');
        
        % Window the data appropriately (into relative and absolute)
        T3_rel = T3(2090763:25072198, :);
        T3_abs = T3(26365236:51698256, :);
        
        % Clear the T3 workspace to conserve memory
        clear T3
        
        % Import data from Biopac (relative and absolute trial periods)
        disp("--> Step 2 of 7: Importing Biopac Data")
        biopac_rel = importdata('Subject2/biopac_1.txt');
        biopac_abs = importdata('Subject2/biopac_2.txt');
        
        % Correct for errors in alignment between Biopac and T3
        T3_rel = T3_rel(3141886:end, :); temp = length(biopac_rel) - length(T3_rel) + 5000;
        biopac_rel = biopac_rel(temp:end, :); biopac_rel = biopac_rel(231:end,:); 
        
        % Assign data to variables (T3)
        P1{1} = T3_rel(:, 2); P1{2} = T3_abs(:, 2);     % Aortic arch pressure
        P2{1} = T3_rel(:, 3); P2{2} = T3_abs(:, 3);     % Femoral artery pressure
        P3{1} = T3_rel(:, 4); P3{2} = T3_abs(:, 4);     % Wedge pressure
        P4{1} = T3_rel(:, 5); P4{2} = T3_abs(:, 5);     % Right atrial pressure
        ECG_t{1} = T3_rel(:, 7); ECG_t{2} = T3_abs(:, 7);	% T3 ECG (2)
        
        % Clear the T3 workspace to conserve memory
        clear T3_rel T3_abs
        
        % Assign data to variables (Biopac)
        ECG_b{1} = biopac_rel(:, 1); ECG_b{2} = biopac_abs(:, 1);	% Biopac ECG
        SPO2_1{1} = biopac_rel(:, 2); SPO2_1{2} = biopac_abs(:, 2);	% Apex SPO2
        SPO2_2{1} = biopac_rel(:, 3); SPO2_2{2} = biopac_abs(:, 3);	% Femoral SPO2
        PPG_1{1} = biopac_rel(:, 4); PPG_1{2} = biopac_abs(:, 4);	% Apex PPG
        PPG_2{1} = biopac_rel(:, 5); PPG_2{2} = biopac_abs(:, 5);	% Femoral PPG
        AX_1{1} = biopac_rel(:, 6); AX_1{2} = biopac_abs(:, 6);     % Sternum SCG (x-acceleration)
        AY_1{1} = biopac_rel(:, 7); AY_1{2} = biopac_abs(:, 7);     % Sternum SCG (y-acceleration)
        AZ_1{1} = biopac_rel(:, 8); AZ_1{2} = biopac_abs(:, 8);     % Sternum SCG (z-acceleration)
        AX_2{1} = biopac_rel(:, 9); AX_2{2} = biopac_abs(:, 9);     % Apex SCG (x-acceleration)
        AY_2{1} = biopac_rel(:, 10); AY_2{2} = biopac_abs(:, 10);   % Apex SCG (y-acceleration)
        AZ_2{1} = biopac_rel(:, 11); AZ_2{2} = biopac_abs(:, 11);	% Apex SCG (z-acceleration)
        
        % Clear the Biopac workspace to conserve memory
        clear biopac_rel biopac_abs
        
        % -----------------------------------------------------------------
        % Clean T3 Data (NaN Removal)
        % -----------------------------------------------------------------
        
        disp("--> Step 3 of 7: Cleaning Data")
        
        % Remove NaNs where necessary
        for k = 1:2
            for j = 1:length(ECG_t{k})
                if isnan(ECG_t{k}(j)); ECG_t{k}(j) = ECG_t{k}(j-1); end
                if isnan(P1{k}(j)); P1{k}(j) = P1{k}(j-1); end
                if isnan(P2{k}(j)); P2{k}(j) = P2{k}(j-1); end
                if isnan(P3{k}(j)); P3{k}(j) = P3{k}(j-1); end
                if isnan(P4{k}(j)); P4{k}(j) = P4{k}(j-1); end
            end
        end
        
        % -----------------------------------------------------------------
        % Filter Data
        % -----------------------------------------------------------------
        
        disp("--> Step 4 of 7: Filtering Data")
        
        if obj.filtered
            
            % ECG
            % Standard Kaiser window BP filter |5/10-25\40|
            Hd = cardio.general.createBPF(5, 40, obj.Fs,'fpass1', 10, 'fpass2', 25, 'kaiser', 'order', 20);
            for j = 1:2
                ECG_b{j} = filtfilt(Hd.Numerator, 1, ECG_b{j});
                ECG_t{j} = filtfilt(Hd.Numerator, 1, ECG_t{j});
            end
            
            % SCG
            % Standard Kaiser window BP filter |0.5/1-20\40|
            Hd = cardio.general.createBPF(0.5, 40, obj.Fs, 'fpass1', 1,'fpass2', 20, 'kaiser', 'order', 20);
            for j = 1:2
                AX_1{j} = filtfilt(Hd.Numerator, 1, AX_1{j}); AX_2{j} = filtfilt(Hd.Numerator, 1, AX_2{j});
                AY_1{j} = filtfilt(Hd.Numerator, 1, AY_1{j}); AY_2{j} = filtfilt(Hd.Numerator, 1, AY_2{j});
                AZ_1{j} = filtfilt(Hd.Numerator, 1, AZ_1{j}); AZ_2{j} = filtfilt(Hd.Numerator, 1, AZ_2{j});
            end
            
            % PPG / Arterial Pressures
            % Standard Kaiser window BP filter |0.5/1-9\10|
            Hd = cardio.general.createBPF(0.5, 10, obj.Fs, 'fpass1', 1, 'fpass2', 9, 'kaiser', 'order', 20);
            for j = 1:2
                PPG_1{j} = filtfilt(Hd.Numerator, 1, PPG_1{j}); PPG_2{j} = filtfilt(Hd.Numerator, 1, PPG_2{j});
                P1{j} = filtfilt(Hd.Numerator, 1, P1{j}); P2{j} = filtfilt(Hd.Numerator, 1, P2{j});
                P3{j} = filtfilt(Hd.Numerator, 1, P3{j}); P4{j} = filtfilt(Hd.Numerator, 1, P4{j});
            end
            
            % PCG
            % Standard Kaiser window HP filter (50Hz)
            Hd = cardio.general.createHPF('Fs', obj.Fs, 'Fc', 50, 'order', 50);
            for j = 1:2; PCG{j} = filtfilt(Hd.Numerator, 1, AZ_1{j}); end
            
        end
        
        % -----------------------------------------------------------------
        % Process ECG (BIOPAC)
        % -----------------------------------------------------------------
        
        disp("--> Step 5 of 7: Processing ECG Data (Biopac)")
        
        % 1. Check polarity. (Are the R-peaks upside down?) 
        % 2. Take the difference, check its shape.
        % 3. Append a 0 as the last element so it is the same length as the
        % other signals. 

        % Differentiate ECG signal
        diff1_ecgb = cell(2, 1); for j = 1:2; diff1_ecgb{j} = diff(-ECG_b{j}); end

        % Append a 0 so it's the same length as everything else
        ECG_bio = cell(2, 1); for j = 1:2; ECG_bio{j} = [diff1_ecgb{j}; zeros(1)]; end
        clear diff1_ecgb

        % Fix outliers so ECG_bio can be processed by cardio.ecg.automatedHR
        ECG_bio_trim = cell(2, 1); limit = 0.015;
        for j = 1:2; ECG_bio_trim{j} = zeros(size(ECG_bio{j})); end
        for k = 1:2
            for j = 1:length(ECG_bio{k})
               if(abs(ECG_bio{k}(j)) > limit)
                   % Keep the relative distance so you don't lose the peaks
                   ECG_bio_trim{k}(j) = (0.01*(abs(ECG_bio{k}(j))-limit)+limit)*sign(ECG_bio{k}(j)); 
               else; ECG_bio_trim{k}(j) = ECG_bio{k}(j); 
               end
            end
        end

        % Clear extra workspace variables
        clear ECG_bio
        
        % Set parameters for further further cleaning and beat-separation
        peak_separation = 840;
        HRV_bound = cell(2, 1); ECG_bt_2 = cell(2, 1);
        HRV_bound{1} = 50; HRV_bound{2} = 2;  % Account for arrhythmias during relative
        for j = 1:2; ECG_bt_2{j} = ECG_bio_trim{j}; end; clear ECG_bio_trim

        % Correct outliers
        for j = 1:length(ECG_bt_2{1})
            if  j > 9019340 && j < 9019440
                ECG_bt_2{1}(j) = ECG_bt_2{1}(j) - 0.015;	% Cut noise
            elseif j > 11849420 && j < 11849490
                ECG_bt_2{1}(j) = ECG_bt_2{1}(j)-0.01;       % Large P wave
            elseif j > 11850470 && j < 11850495
                ECG_bt_2{1}(j) = ECG_bt_2{1}(j)-0.01;       % Large T wave
            elseif j > 18287000 && j < 18288133
                ECG_bt_2{1}(j) = ECG_bt_2{1}(j)*0.1;        % Cut noise
            elseif j > 18288200 && j < 18288320
                ECG_bt_2{1}(j) = ECG_bt_2{1}(j)*0.1;        % Cut noise
            end
        end

        % Extract heart rate and beat-separated ECG
        beats_b = cell(2, 1);
        for j = 1:2
            [HR_b{j}, beats_b{j}] = cardio.ecg.cleanECG(ECG_bt_2{j}, obj.Fs, ...
                HRV_bound{j}, 'peak_sep', peak_separation, 'double', 'plotsoff');
        end
        
        % -----------------------------------------------------------------
        % Process ECG (T3)
        % -----------------------------------------------------------------
        
        disp("--> Step 6 of 7: Processing ECG Data (T3)")
        
        % Set bounds and cutoffs
        bound = 0.5; cutoff = 23800000; HRV_bound{1} = 2;
        
        % Bound the ECG signal
        ECG_t3_trim = cell(2, 1); ECG_t3_trim{1} = -ECG_t{1};
        ECG_t3_trim{1}(ECG_t3_trim{1} > bound) = ...
            ECG_t3_trim{1}(ECG_t3_trim{1} > bound) + 0.01*bound;
        
        % Remove the death sequence
        ECG_t3_trim{2} = -ECG_t{2}(1:cutoff);

        % Rename workspace variable
        ECG_tt_2 = ECG_t3_trim; clear ECG_t3_trim
        
        % Correct the errors in the second trial
        for j = 1:length(ECG_tt_2{2})
            if j > 3345650 && j < 3345680
                ECG_tt_2{2}(j) = ECG_tt_2{2}(j) + 0.45;	% Low QRS
            elseif j > 4386580 && j < 4386620
                ECG_tt_2{2}(j) = ECG_tt_2{2}(j) + 0.18;	% Low QRS
            elseif j > 4496665 && j < 4496695
                ECG_tt_2{2}(j) = ECG_tt_2{2}(j) + 0.5;     % Low QRS
            end
        end

        % Extract heart rate and beat-separated ECG
        beats_t = cell(2, 1);
        for j = 1:2
            [HR_t{j}, beats_t{j}] = cardio.ecg.cleanECG(ECG_tt_2{j}, obj.Fs, ...
                HRV_bound{j}, 'peak_sep', peak_separation, 'double', 'plotsoff');
        end
        
        % Because you cut out more noise at the end of T3, shorten the bio
        % heartbeats off the end so you're working with the same number of beats: 
        
        for j = 1:2
            beats_t{j} = double(beats_t{j}); HR_t{j} = double(HR_t{j});
            beats_b{j} = double(beats_b{j}(1:length(beats_t{j})));
            HR_b{j} = double(HR_b{j}(1:length(HR_t{j})));
        end
        % Bug fix
        for j = 1:2
            if length(HR_t{j}) < length(beats_t{j}); HR_t{j}(end + 1) = HR_t{j}(end); end
            if length(HR_b{j}) < length(beats_b{j}); HR_b{j}(end + 1) = HR_b{j}(end); end
        end
        
        % -----------------------------------------------------------------
        % Append Dataset by Label
        % -----------------------------------------------------------------
        
        disp("--> Step 7 of 7: Appending Dataset by Labels")
        
        % Call trimToIdx function
        trimToIdx_s2(2, obj.levels, include);
        
        % For each level this subject has, find the beat indices
        allLevelsToBeats(2, include);
        
    end

%% ------------------------------------------------------------------------
% Sub-Function for Subject 3
% -------------------------------------------------------------------------

    function subject3()
        
        % NOTE: Subject 3 does not have interval Level.all
        % NOTE: The first element in any cell array is from the RELATIVE
        % trial while the second element is from the ABSOLUTE trial.
    
        % Returns dataset for Subject 3 as a struct:
        % obj.dataset{3}.(interval).ecg_biopac
        % obj.dataset{3}.(interval).ecg_t3
        % obj.dataset{3}.(interval).sternumSCG_x
        % obj.dataset{3}.(interval).sternumSCG_y
        % obj.dataset{3}.(interval).sternumSCG_z
        % obj.dataset{3}.(interval).sternumSCG_z_back
        % obj.dataset{3}.(interval).apexSCG_x
        % obj.dataset{3}.(interval).apexSCG_y
        % obj.dataset{3}.(interval).apexSCG_z
        % obj.dataset{3}.(interval).apexSCG_z_back
        % obj.dataset{3}.(interval).apexSPO2
        % obj.dataset{3}.(interval).femoralSPO2
        % obj.dataset{3}.(interval).apexPPG
        % obj.dataset{3}.(interval).femoralPPG
        % obj.dataset{3}.(interval).aorticPressure
        % obj.dataset{3}.(interval).femoralPressure
        % obj.dataset{3}.(interval).wedgePressure
        % obj.dataset{3}.(interval).rightAtrialPressure
        
        % -----------------------------------------------------------------
        % Load Dataset and Setup Workspace
        % -----------------------------------------------------------------

        % Import data from T3
        disp("--> Step 1 of 7: Importing T3 Data")
        T3 = importdata('Subject3/t3.txt');
        
        % Window the data appropriately (into relative and absolute)
        T3_rel = T3(2284208:19081280, :);
        T3_abs = T3(20006844:40624037, :);
        
        % Clear the T3 workspace to conserve memory
        clear T3
        
        % Cut beginning of T3 recording (wrong sampling frequency)
        T3_rel = T3_rel(283844:end, :);
        
        % Import data from Biopac (relative and absolute trial periods)
        disp("--> Step 2 of 7: Importing Biopac Data")
        biopac_rel = importdata('Subject3/biopac_1.txt');
        biopac_abs = importdata('Subject3/biopac_2.txt');
        
        % Correct for errors in alignment between Biopac and T3
        biopac_rel = biopac_rel(601050:end, :);
        
        % Assign data to variables (T3)
        P1{1} = T3_rel(:, 2); P1{2} = T3_abs(:, 2);     % Aortic arch pressure
        P2{1} = T3_rel(:, 3); P2{2} = T3_abs(:, 3);     % Femoral artery pressure
        P3{1} = T3_rel(:, 4); P3{2} = T3_abs(:, 4);     % Wedge pressure
        P4{1} = T3_rel(:, 5); P4{2} = T3_abs(:, 5);     % Right atrial pressure
        ECG_t{1} = T3_rel(:, 6); ECG_t{2} = T3_abs(:, 6);	% T3 ECG (1)
        
        % Clear the T3 workspace to conserve memory
        clear T3_rel T3_abs
        
        % Assign data to variables (Biopac)
        ECG_b{1} = biopac_rel(:, 1); ECG_b{2} = biopac_abs(:, 1);	% Biopac ECG
        SPO2_1{1} = biopac_rel(:, 2); SPO2_1{2} = biopac_abs(:, 2);	% Apex SPO2
        SPO2_2{1} = biopac_rel(:, 3); SPO2_2{2} = biopac_abs(:, 3);	% Femoral SPO2
        PPG_1{1} = biopac_rel(:, 4); PPG_1{2} = biopac_abs(:, 4);	% Apex PPG
        PPG_2{1} = biopac_rel(:, 5); PPG_2{2} = biopac_abs(:, 5);	% Femoral PPG
        AX_1{1} = biopac_rel(:, 6); AX_1{2} = biopac_abs(:, 6);     % Sternum SCG (x-acceleration)
        AY_1{1} = biopac_rel(:, 7); AY_1{2} = biopac_abs(:, 7);     % Sternum SCG (y-acceleration)
        AZ_1{1} = biopac_rel(:, 8); AZ_1{2} = biopac_abs(:, 8);     % Sternum SCG (z-acceleration)
        AX_2{1} = biopac_rel(:, 9); AX_2{2} = biopac_abs(:, 9);     % Apex SCG (x-acceleration)
        AY_2{1} = biopac_rel(:, 10); AY_2{2} = biopac_abs(:, 10);   % Apex SCG (y-acceleration)
        AZ_2{1} = biopac_rel(:, 11); AZ_2{2} = biopac_abs(:, 11);	% Apex SCG (z-acceleration)
        
        % Clear the Biopac workspace to conserve memory
        clear biopac_rel biopac_abs
        
        % -----------------------------------------------------------------
        % Clean T3 Data (NaN Removal)
        % -----------------------------------------------------------------
        
        disp("--> Step 3 of 7: Cleaning Data")
        
        % Remove NaNs where necessary
        for k = 1:2
            for j = 1:length(ECG_t{k})
                if isnan(ECG_t{k}(j)); ECG_t{k}(j) = ECG_t{k}(j-1); end
                if isnan(P1{k}(j)); P1{k}(j) = P1{k}(j-1); end
                if isnan(P2{k}(j)); P2{k}(j) = P2{k}(j-1); end
                if isnan(P3{k}(j)); P3{k}(j) = P3{k}(j-1); end
                if isnan(P4{k}(j)); P4{k}(j) = P4{k}(j-1); end
            end
        end
        
        % -----------------------------------------------------------------
        % Filter Data
        % -----------------------------------------------------------------
        
        disp("--> Step 4 of 7: Filtering Data")
        
        if obj.filtered
            
            % ECG
            % Standard Kaiser window BP filter |5/10-25\40|
            Hd = cardio.general.createBPF(5, 40, obj.Fs,'fpass1', 10, 'fpass2', 25, 'kaiser', 'order', 20);
            for j = 1:2
                ECG_b{j} = filtfilt(Hd.Numerator, 1, ECG_b{j});
                for k = 1:4; ECG_t{j} = filtfilt(Hd.Numerator, 1, ECG_t{j}); end
            end
            
            % SCG
            % Standard Kaiser window BP filter |0.5/1-20\40|
            Hd = cardio.general.createBPF(0.5, 40, obj.Fs, 'fpass1', 1,'fpass2', 20, 'kaiser', 'order', 20);
            for j = 1:2
                AX_1{j} = filtfilt(Hd.Numerator, 1, AX_1{j}); AX_2{j} = filtfilt(Hd.Numerator, 1, AX_2{j});
                AY_1{j} = filtfilt(Hd.Numerator, 1, AY_1{j}); AY_2{j} = filtfilt(Hd.Numerator, 1, AY_2{j});
                AZ_1{j} = filtfilt(Hd.Numerator, 1, AZ_1{j}); AZ_2{j} = filtfilt(Hd.Numerator, 1, AZ_2{j});
            end
            
            % PPG / Arterial Pressures
            % Standard Kaiser window BP filter |0.5/1-9\10|
            Hd = cardio.general.createBPF(0.5, 10, obj.Fs, 'fpass1', 1, 'fpass2', 9, 'kaiser', 'order', 20);
            for j = 1:2
                PPG_1{j} = filtfilt(Hd.Numerator, 1, PPG_1{j}); PPG_2{j} = filtfilt(Hd.Numerator, 1, PPG_2{j});
                P1{j} = filtfilt(Hd.Numerator, 1, P1{j}); P2{j} = filtfilt(Hd.Numerator, 1, P2{j});
                P3{j} = filtfilt(Hd.Numerator, 1, P3{j}); P4{j} = filtfilt(Hd.Numerator, 1, P4{j});
            end
            
            % PCG
            % Standard Kaiser window HP filter (50Hz)
            Hd = cardio.general.createHPF('Fs', obj.Fs, 'Fc', 50, 'order', 50);
            for j = 1:2; PCG{j} = filtfilt(Hd.Numerator, 1, AZ_1{j}); end
            
        end
        
        % -----------------------------------------------------------------
        % Process ECG (BIOPAC)
        % -----------------------------------------------------------------
        
        disp("--> Step 5 of 7: Processing ECG Data (Biopac)")
        
        % 1. Check polarity. (Are the R-peaks upside down?) 
        % 2. Take the difference, check its shape.
        % 3. Append a 0 as the last element so it is the same length as the
        % other signals. 

        % Differentiate ECG signal
        diff1_ecgb = cell(2, 1); for j = 1:2; diff1_ecgb{j} = diff(ECG_b{j}); end

        % Append a 0 so it's the same length as everything else
        ECG_bio = cell(2, 1); for j = 1:2; ECG_bio{j} = [diff1_ecgb{j}; zeros(1)]; end
        clear diff1_ecgb

        % Fix outliers so ECG_bio can be processed by cardio.ecg.automatedHR
        ECG_bio_trim = cell(2, 1); limit = 0.015;
        for j = 1:2; ECG_bio_trim{j} = zeros(size(ECG_bio{j})); end
        for k = 1:2
            for j = 1:length(ECG_bio{k})
               if(abs(ECG_bio{k}(j)) > limit)
                   % Keep the relative distance so you don't lose the peaks
                   ECG_bio_trim{k}(j) = (0.01*(abs(ECG_bio{k}(j))-limit)+limit)*sign(ECG_bio{k}(j)); 
               else; ECG_bio_trim{k}(j) = ECG_bio{k}(j); 
               end
            end
        end

        % Clear extra workspace variables
        clear ECG_bio
        
        % Set parameters for further further cleaning and beat-separation
        peak_separation = 750;
        ECG_bt_2 = cell(2, 1); HRV_bound = 1;  % Account for arrhythmias during relative
        for j = 1:2; ECG_bt_2{j} = ECG_bio_trim{j}; end; clear ECG_bio_trim

        % Correct outliers
        ECG_bt_2{1}(1775551:1775801) = ECG_bt_2{1}(1775551:1775801)*0.1;        % Large T-wave
        ECG_bt_2{1}(7759651:7759801) = ECG_bt_2{1}(7759651:7759801)*0.1;        % Large T-wave
        ECG_bt_2{2}(20200520:20200560) = ECG_bt_2{2}(20200520:20200560)*0.1;    % Large S overcompensation
        ECG_bt_2{2}(17851440:17851480) = ECG_bt_2{2}(17851440:17851480)*0.1;    % Large S overcompensation
        ECG_bt_2{2}(12349900:12350020) = ECG_bt_2{2}(12349900:12350020)*0.1;    % PVC removed
        ECG_bt_2{2}(17847050:17847100) = ECG_bt_2{2}(17847050:17847100)*0.1;    % PCC removed
        ECG_bt_2{2}(18094750:18094900) = ECG_bt_2{2}(18094750:18094900)*0.1;    % PVC removed
        ECG_bt_2{2}(19252920:19253050) = ECG_bt_2{2}(19252920:19253050)*0.1;    % PVC removed

        % Extract heart rate and beat-separated ECG
        beats_b = cell(2, 1);
        for j = 1:2
            [HR_b{j}, beats_b{j}] = cardio.ecg.cleanECG(ECG_bt_2{j}, obj.Fs, ...
                HRV_bound, 'peak_sep', peak_separation, 'double', 'plotsoff');
        end
        
        % -----------------------------------------------------------------
        % Process ECG (T3)
        % -----------------------------------------------------------------
        
        disp("--> Step 6 of 7: Processing ECG Data (T3)")
        
        % Set bounds and variables
        HRV_bound = 75; ECG_t3_trim = ECG_t;

        % Rename workspace variable
        ECG_tt_2 = cell(2, 1); for j = 1:2; ECG_tt_2{j} = diff(ECG_t3_trim{j}); end
        clear ECG_t3_trim
        
        % Correct for outliers
        ECG_tt_2{1}(7396827:7396847) = ECG_tt_2{1}(7396827:7396847) + 0.0015;   % Low QRS
        ECG_tt_2{1}(143077:143117) = ECG_tt_2{1}(143077:143117) + 0.0015;       % Low QRS
        ECG_tt_2{2}(12349400:12349520) = ECG_tt_2{2}(12349400:12349520)*0.1;    % Removed PVC
        ECG_tt_2{2}(17846300:17846330) = ECG_tt_2{2}(17846300:17846330)*0.1;    % Removed PVC
        ECG_tt_2{2}(18094030:18094080) = ECG_tt_2{2}(18094030:18094080)*0.1;    % Removed PVC
        ECG_tt_2{2}(19252100:19252450) = ECG_tt_2{2}(19252100:19252450)*0.1;    % Removed PVC

        % Extract heart rate and beat-separated ECG
        beats_t = cell(2, 1);
        for j = 1:2
            [HR_t{j}, beats_t{j}] = cardio.ecg.cleanECG(ECG_tt_2{j}, obj.Fs, ...
                HRV_bound, 'peak_sep', peak_separation, 'double', 'plotsoff');
        end
        
        % Because you cut out more noise at the end of T3, shorten the bio
        % heartbeats off the end so you're working with the same number of beats: 
        
        % Relative hypovolemia
        beats_t{1} = double(beats_t{1}); HR_t{1} = double(HR_t{1});
        HR_b{1} = double(HR_b{1});  beats_b{1} = double(beats_b{1});
        
        % Absolute hypovolemia
        beats_b{2} = double(beats_b{2}(1:length(beats_t{2})));
        HR_b{2} = double(HR_b{2}(1:length(HR_t{2})));
        beats_t{2} = double(beats_t{2}); HR_t{2} = double(HR_t{2});
        % Bug fix
        for j = 1:2
            if length(HR_t{j}) < length(beats_t{j}); HR_t{j}(end + 1) = HR_t{j}(end); end
            if length(HR_b{j}) < length(beats_b{j}); HR_b{j}(end + 1) = HR_b{j}(end); end
        end
        
        % -----------------------------------------------------------------
        % Append Dataset by Label
        % -----------------------------------------------------------------
        
        disp("--> Step 7 of 7: Appending Dataset by Labels")
        
        % Call trimToIdx function
        trimToIdx_s2(3, obj.levels, include);
        
        % For each level this subject has, find the beat indices
        allLevelsToBeats(3, include);
        
    end

%% ------------------------------------------------------------------------
% Sub-Function for Subject 4
% -------------------------------------------------------------------------

    function subject4()

        % NOTE: Subject 4 does not have interval Level.all
        % NOTE: The first element in any cell array is from the RELATIVE
        % trial while the second element is from the ABSOLUTE trial.

        % Returns dataset for Subject 4 as a struct:
        % obj.dataset{4}.(interval).ecg_biopac
        % obj.dataset{4}.(interval).ecg_t3
        % obj.dataset{4}.(interval).sternumSCG_x
        % obj.dataset{4}.(interval).sternumSCG_y
        % obj.dataset{4}.(interval).sternumSCG_z
        % obj.dataset{4}.(interval).sternumSCG_z_back
        % obj.dataset{4}.(interval).apexSCG_x
        % obj.dataset{4}.(interval).apexSCG_y
        % obj.dataset{4}.(interval).apexSCG_z
        % obj.dataset{4}.(interval).apexSCG_z_back
        % obj.dataset{4}.(interval).apexSPO2
        % obj.dataset{4}.(interval).femoralSPO2
        % obj.dataset{4}.(interval).apexPPG
        % obj.dataset{4}.(interval).femoralPPG
        % obj.dataset{4}.(interval).aorticPressure
        % obj.dataset{4}.(interval).femoralPressure
        % obj.dataset{4}.(interval).wedgePressure
        % obj.dataset{4}.(interval).rightAtrialPressure

        % -----------------------------------------------------------------
        % Load Dataset and Setup Workspace
        % -----------------------------------------------------------------

        % Import data from T3
        disp("--> Step 1 of 7: Importing T3 Data")
        T3 = importdata('Subject4/t3.txt');

        % Window the data appropriately (into relative and absolute)
        T3_rel = T3(5230292:13834583, :);
        T3_abs = T3(14239365:34048624, :);

        % Clear the T3 workspace to conserve memory
        clear T3

        % Import data from Biopac (relative and absolute trial periods)
        disp("--> Step 2 of 7: Importing Biopac Data")
        biopac_rel = importdata('Subject4/biopac_1.txt');
        biopac_abs = importdata('Subject4/biopac_2.txt');

        % Assign data to variables (T3)
        P1{1} = T3_rel(:, 3); P1{2} = T3_abs(:, 3);     % Aortic arch pressure
        P2{1} = T3_rel(:, 2); P2{2} = T3_abs(:, 2);     % Femoral artery pressure
        P3{1} = T3_rel(:, 4); P3{2} = T3_abs(:, 4);     % Wedge pressure
        P4{1} = T3_rel(:, 5); P4{2} = T3_abs(:, 5);     % Right atrial pressure
        ECG_t{1} = T3_rel(:, 6); ECG_t{2} = T3_abs(:, 6);	% T3 ECG (1)

        % Clear the T3 workspace to conserve memory
        clear T3_rel T3_abs

        % Assign data to variables (Biopac)
        ECG_b{1} = biopac_rel(:, 1); ECG_b{2} = biopac_abs(:, 1);	% Biopac ECG
        SPO2_1{1} = biopac_rel(:, 2); SPO2_1{2} = biopac_abs(:, 2);	% Apex SPO2
        SPO2_2{1} = biopac_rel(:, 3); SPO2_2{2} = biopac_abs(:, 3);	% Femoral SPO2
        PPG_1{1} = biopac_rel(:, 4); PPG_1{2} = biopac_abs(:, 4);	% Apex PPG
        PPG_2{1} = biopac_rel(:, 5); PPG_2{2} = biopac_abs(:, 5);	% Femoral PPG
        AX_1{1} = biopac_rel(:, 6); AX_1{2} = biopac_abs(:, 6);     % Sternum SCG (x-acceleration)
        AY_1{1} = biopac_rel(:, 7); AY_1{2} = biopac_abs(:, 7);     % Sternum SCG (y-acceleration)
        AZ_1{1} = biopac_rel(:, 8); AZ_1{2} = biopac_abs(:, 8);     % Sternum SCG (z-acceleration)
        AX_2{1} = biopac_rel(:, 9); AX_2{2} = biopac_abs(:, 9);     % Apex SCG (x-acceleration)
        AY_2{1} = biopac_rel(:, 10); AY_2{2} = biopac_abs(:, 10);   % Apex SCG (y-acceleration)
        AZ_2{1} = biopac_rel(:, 11); AZ_2{2} = biopac_abs(:, 11);	% Apex SCG (z-acceleration)

        % Clear the Biopac workspace to conserve memory
        clear biopac_rel biopac_abs

        % -----------------------------------------------------------------
        % Clean T3 Data (NaN Removal)
        % -----------------------------------------------------------------

        disp("--> Step 3 of 7: Cleaning Data")

        % Remove NaNs where necessary
        for k = 1:2
            for j = 1:length(ECG_t{k})
                if isnan(ECG_t{k}(j)); ECG_t{k}(j) = ECG_t{k}(j-1); end
                if isnan(P1{k}(j)); P1{k}(j) = P1{k}(j-1); end
                if isnan(P2{k}(j)); P2{k}(j) = P2{k}(j-1); end
                if isnan(P3{k}(j)); P3{k}(j) = P3{k}(j-1); end
                if isnan(P4{k}(j)); P4{k}(j) = P4{k}(j-1); end
            end
        end

        % -----------------------------------------------------------------
        % Filter Data
        % -----------------------------------------------------------------

        disp("--> Step 4 of 7: Filtering Data")

        if obj.filtered

            % ECG
            % Standard Kaiser window BP filter |5/10-25\40|
            Hd = cardio.general.createBPF(5, 40, obj.Fs,'fpass1', 10, 'fpass2', 25, 'kaiser', 'order', 20);
            for j = 1:2
                ECG_b{j} = filtfilt(Hd.Numerator, 1, ECG_b{j});
                for k = 1:4; ECG_t{j} = filtfilt(Hd.Numerator, 1, ECG_t{j}); end
            end

            % SCG
            % Standard Kaiser window BP filter |0.5/1-20\40|
            Hd = cardio.general.createBPF(0.5, 40, obj.Fs, 'fpass1', 1,'fpass2', 20, 'kaiser', 'order', 20);
            for j = 1:2
                AX_1{j} = filtfilt(Hd.Numerator, 1, AX_1{j}); AX_2{j} = filtfilt(Hd.Numerator, 1, AX_2{j});
                AY_1{j} = filtfilt(Hd.Numerator, 1, AY_1{j}); AY_2{j} = filtfilt(Hd.Numerator, 1, AY_2{j});
                AZ_1{j} = filtfilt(Hd.Numerator, 1, AZ_1{j}); AZ_2{j} = filtfilt(Hd.Numerator, 1, AZ_2{j});
            end

            % PPG / Arterial Pressures
            % Standard Kaiser window BP filter |0.5/1-9\10|
            Hd = cardio.general.createBPF(0.5, 10, obj.Fs, 'fpass1', 1, 'fpass2', 9, 'kaiser', 'order', 20);
            for j = 1:2
                PPG_1{j} = filtfilt(Hd.Numerator, 1, PPG_1{j}); PPG_2{j} = filtfilt(Hd.Numerator, 1, PPG_2{j});
                P1{j} = filtfilt(Hd.Numerator, 1, P1{j}); P2{j} = filtfilt(Hd.Numerator, 1, P2{j});
                P3{j} = filtfilt(Hd.Numerator, 1, P3{j}); P4{j} = filtfilt(Hd.Numerator, 1, P4{j});
            end

            % PCG
            % Standard Kaiser window HP filter (50Hz)
            Hd = cardio.general.createHPF('Fs', obj.Fs, 'Fc', 50, 'order', 50);
            for j = 1:2; PCG{j} = filtfilt(Hd.Numerator, 1, AZ_1{j}); end

        end

        % -----------------------------------------------------------------
        % Process ECG (BIOPAC)
        % -----------------------------------------------------------------

        disp("--> Step 5 of 7: Processing ECG Data (Biopac)")

        % 1. Check polarity. (Are the R-peaks upside down?) 
        % 2. Take the difference, check its shape.
        % 3. Append a 0 as the last element so it is the same length as the
        % other signals. 

        % Differentiate ECG signal
        diff1_ecgb = cell(2, 1); for j = 1:2; diff1_ecgb{j} = diff(ECG_b{j}); end

        % Append a 0 so it's the same length as everything else
        ECG_bio = cell(2, 1); for j = 1:2; ECG_bio{j} = [diff1_ecgb{j}; zeros(1)]; end
        clear diff1_ecgb

        % Fix outliers so ECG_bio can be processed by cardio.ecg.automatedHR
        ECG_bio_trim = cell(2, 1); limit = 0.015;
        for j = 1:2; ECG_bio_trim{j} = zeros(size(ECG_bio{j})); end
        for k = 1:2
            for j = 1:length(ECG_bio{k})
               if(abs(ECG_bio{k}(j)) > limit)
                   % Keep the relative distance so you don't lose the peaks
                   ECG_bio_trim{k}(j) = (0.01*(abs(ECG_bio{k}(j))-limit)+limit)*sign(ECG_bio{k}(j)); 
               else; ECG_bio_trim{k}(j) = ECG_bio{k}(j); 
               end
            end
        end

        % Clear extra workspace variables
        clear ECG_bio

        % Set parameters for further further cleaning and beat-separation
        peak_separation = 750;
        ECG_bt_2 = cell(2, 1); HRV_bound = cell(2, 1); HRV_bound{1} = 2; HRV_bound{2} = 75;
        for j = 1:2; ECG_bt_2{j} = ECG_bio_trim{j}; end; clear ECG_bio_trim

        % Fix noise at end of signal
        ECG_bt_2{2} = fix_diff(ECG_bt_2{2});
        
        % Correct outliers (Relative)
        ECG_bt_2{1}(740535:740570) = ECG_bt_2{1}(740535:740570)*0.1;
        ECG_bt_2{1}(742035:742080) = ECG_bt_2{1}(742035:742080)*0.1;
        ECG_bt_2{1}(2379640:2379680) = ECG_bt_2{1}(2379640:2379680)*0.1;
        ECG_bt_2{1}(6801510:6801550) = ECG_bt_2{1}(6801510:6801550)*0.1;
        ECG_bt_2{1}(7903410:7903460) = ECG_bt_2{1}(7903410:7903460)*0.1;
        ECG_bt_2{1}(8363650:8363690) = ECG_bt_2{1}(8363650:8363690)*0.1;
%         ECG_bt_2{1}(16207280:16207310) = ECG_bt_2{1}(16207280:16207310)*0.1;
%         ECG_bt_2{1}(16496590:16496630) = ECG_bt_2{1}(16496590:16496630)*0.1;

        % Correct outliers (Absolute)
        ECG_bt_2{2}(8265150:8265320) = ECG_bt_2{2}(8265150:8265320)*0.1;
        ECG_bt_2{2}(8285030:8285180) = ECG_bt_2{2}(8285030:8285180)*0.1;
        ECG_bt_2{2}(8290760:8290900) = ECG_bt_2{2}(8290760:8290900)*0.1;
        ECG_bt_2{2}(8292670:8292800) = ECG_bt_2{2}(8292670:8292800)*0.1;
        ECG_bt_2{2}(8300420:8300570) = ECG_bt_2{2}(8300420:8300570)*0.1;
        ECG_bt_2{2}(8303040:8303120) = ECG_bt_2{2}(8303040:8303120)*0.1;
        ECG_bt_2{2}(8312700:8312820) = ECG_bt_2{2}(8312700:8312820)*0.1;
        ECG_bt_2{2}(8335180:8335300) = ECG_bt_2{2}(8335180:8335300)*0.1;
        ECG_bt_2{2}(8338370:8338500) = ECG_bt_2{2}(8338370:8338500)*0.1;
        ECG_bt_2{2}(8376400:8376520) = ECG_bt_2{2}(8376400:8376520)*0.1;
        ECG_bt_2{2}(8378200:8378300) = ECG_bt_2{2}(8378200:8378300)*0.1;
        ECG_bt_2{2}(8380010:8380120) = ECG_bt_2{2}(8380010:8380120)*0.1;
        ECG_bt_2{2}(8456110:8456250) = ECG_bt_2{2}(8456110:8456250)*0.1;
        ECG_bt_2{2}(8462250:8462430) = ECG_bt_2{2}(8462250:8462430)*0.1;
        ECG_bt_2{2}(8464380:8464500) = ECG_bt_2{2}(8464380:8464500)*0.1;
        ECG_bt_2{2}(8466330:8466400) = ECG_bt_2{2}(8466330:8466400)*0.1;
        ECG_bt_2{2}(8520850:8520950) = ECG_bt_2{2}(8520850:8520950)*0.1;
        ECG_bt_2{2}(8539530:8539650) = ECG_bt_2{2}(8539530:8539650)*0.1;
        ECG_bt_2{2}(8547750:8547840) = ECG_bt_2{2}(8547750:8547840)*0.1;
        ECG_bt_2{2}(8557850:8557950) = ECG_bt_2{2}(8557850:8557950)*0.1;
        ECG_bt_2{2}(8596970:8597060) = ECG_bt_2{2}(8596970:8597060)*0.1;
        ECG_bt_2{2}(8598850:8598960) = ECG_bt_2{2}(8598850:8598960)*0.1;
        ECG_bt_2{2}(8269020:8269140) = ECG_bt_2{2}(8269020:8269140)*0.1;
        ECG_bt_2{2}(8563950:8564100) = ECG_bt_2{2}(8563950:8564100)*0.1;

        % Extract heart rate and beat-separated ECG
        beats_b = cell(2, 1);
        for j = 1:2
            [HR_b{j}, beats_b{j}] = cardio.ecg.cleanECG(ECG_bt_2{j}, obj.Fs, ...
                HRV_bound{j}, 'peak_sep', peak_separation, 'double', 'plotsoff');
        end

        % -----------------------------------------------------------------
        % Process ECG (T3)
        % -----------------------------------------------------------------

        disp("--> Step 6 of 7: Processing ECG Data (T3)")

        % Set bounds and variables
        HRV_bound = 75; ECG_t3_trim = ECG_t;

        % Rename workspace variable
        ECG_tt_2 = cell(2, 1); for j = 1:2; ECG_tt_2{j} = diff(ECG_t3_trim{j}); end
        clear ECG_t3_trim

        % Correct for outliers (Relative)
        ECG_tt_2{1}(7358200:7359970) = ECG_tt_2{1}(7358200:7359970)*0.1;
        ECG_tt_2{1}(7360100:7361800) = ECG_tt_2{1}(7360100:7361800)*0.1;
        ECG_tt_2{1}(7362000:7363740) = ECG_tt_2{1}(7362000:7363740)*0.1;
        ECG_tt_2{1}(7363800:7365600) = ECG_tt_2{1}(7363800:7365600)*0.1;
        ECG_tt_2{1}(7365700:7367400) = ECG_tt_2{1}(7365700:7367400)*0.1;
        ECG_tt_2{1}(7367600:7369200) = ECG_tt_2{1}(7367600:7369200)*0.1;
        ECG_tt_2{1}(6844600:6844850) = ECG_tt_2{1}(6844600:6844850)*0.1;
        ECG_tt_2{1}(5154000:5154150) = ECG_tt_2{1}(5154000:5154150)*0.1;
        ECG_tt_2{1}(3505300:3505400) = ECG_tt_2{1}(3505300:3505400)*0.1;
        ECG_tt_2{1}(3505500:3505550) = ECG_tt_2{1}(3505500:3505550)*4.0;

        % Correct for outliers (Absolute)
        ECG_tt_2{2}(7651580:7651680) = ECG_tt_2{2}(7651580:7651680)*0.1;
        ECG_tt_2{2}(1115000:1115200) = ECG_tt_2{2}(1115000:1115200)*0.1;
        ECG_tt_2{2}(8284700:8284750) = ECG_tt_2{2}(8284700:8284750)*0.1;
        ECG_tt_2{2}(8376050:8376150) = ECG_tt_2{2}(8376050:8376150)*0.1;
        ECG_tt_2{2}(8377820:8377920) = ECG_tt_2{2}(8377820:8377920)*0.1;
        ECG_tt_2{2}(8379650:8379750) = ECG_tt_2{2}(8379650:8379750)*0.1;
        ECG_tt_2{2}(8455750:8455800) = ECG_tt_2{2}(8455750:8455800)*0.1;
        ECG_tt_2{2}(8461920:8461990) = ECG_tt_2{2}(8461920:8461990)*0.1;
        ECG_tt_2{2}(8334780:8334840) = ECG_tt_2{2}(8334780:8334840)*0.1;
        ECG_tt_2{2}(8334780:8334840) = ECG_tt_2{2}(8334780:8334840)*0.1;
        ECG_tt_2{2}(8337980:8338050) = ECG_tt_2{2}(8337980:8338050)*0.1;
        ECG_tt_2{2}(8465900:8466000) = ECG_tt_2{2}(8465900:8466000)*0.1;
        ECG_tt_2{2}(8598430:8598510) = ECG_tt_2{2}(8598430:8598510)*0.1;
        ECG_tt_2{2}(8290340:8290430) = ECG_tt_2{2}(8290340:8290430)*0.1;
        ECG_tt_2{2}(8264875:8264970) = ECG_tt_2{2}(8264875:8264970)*0.1;
        ECG_tt_2{2}(8297210:8297260) = ECG_tt_2{2}(8297210:8297260) + 10;
        ECG_tt_2{2}(8299060:8299090) = ECG_tt_2{2}(8299060:8299090) + 10;
        ECG_tt_2{2}(8300940:8301000) = ECG_tt_2{2}(8300940:8301000) + 10;
        ECG_tt_2{2}(8301870:8301910) = ECG_tt_2{2}(8301870:8301910) + 10;
        ECG_tt_2{2}(8313170:8313210) = ECG_tt_2{2}(8313170:8313210) + 10;
        ECG_tt_2{2}(8336980:8337030) = ECG_tt_2{2}(8336980:8337030) + 10;
        ECG_tt_2{2}(8454830:8454860) = ECG_tt_2{2}(8454830:8454860) + 10;
        ECG_tt_2{2}(8460950:8461030) = ECG_tt_2{2}(8460950:8461030) + 10;
        ECG_tt_2{2}(8463080:8463130) = ECG_tt_2{2}(8463080:8463130) + 10;
        ECG_tt_2{2}(8468520:8468560) = ECG_tt_2{2}(8468520:8468560) + 10;
        ECG_tt_2{2}(8469550:8469590) = ECG_tt_2{2}(8469550:8469590) + 10;
        ECG_tt_2{2}(8465900:8466000) = ECG_tt_2{2}(8465900:8466000) + 10;
        ECG_tt_2{2}(8473330:8473380) = ECG_tt_2{2}(8473330:8473380) + 10;
        ECG_tt_2{2}(8483670:8483720) = ECG_tt_2{2}(8483670:8483720) + 10;
        ECG_tt_2{2}(8487830:8487880) = ECG_tt_2{2}(8487830:8487880) + 10;
        ECG_tt_2{2}(8489165:8489200) = ECG_tt_2{2}(8489165:8489200) + 10;
        ECG_tt_2{2}(8490560:8490600) = ECG_tt_2{2}(8490560:8490600) + 10;
        ECG_tt_2{2}(8491780:8491810) = ECG_tt_2{2}(8491780:8491810) + 10;
        ECG_tt_2{2}(8493650:8493690) = ECG_tt_2{2}(8493650:8493690) + 10;
        ECG_tt_2{2}(8515190:8515225) = ECG_tt_2{2}(8515190:8515225) + 10;
        ECG_tt_2{2}(8526790:8526830) = ECG_tt_2{2}(8526790:8526830) + 10;
        ECG_tt_2{2}(8553540:8553600) = ECG_tt_2{2}(8553540:8553600) + 10;
        ECG_tt_2{2}(9043260:9043300) = ECG_tt_2{2}(9043260:9043300) + 10;

        % Extract heart rate and beat-separated ECG
        beats_t = cell(2, 1);
        for j = 1:2
            [HR_t{j}, beats_t{j}] = cardio.ecg.cleanECG(ECG_tt_2{j}, obj.Fs, ...
                HRV_bound, 'peak_sep', peak_separation, 'double', 'plotsoff');
        end

        % Because you cut out more noise at the end of T3, shorten the bio
        % heartbeats off the end so you're working with the same number of beats: 

        % Relative hypovolemia
        beats_t{1} = double(beats_t{1}); HR_t{1} = double(HR_t{1});
        HR_b{1} = double(HR_b{1});  beats_b{1} = double(beats_b{1});

        % Absolute hypovolemia
        beats_b{2} = double(beats_b{2}(1:length(beats_t{2})));
        HR_b{2} = double(HR_b{2}(1:length(HR_t{2})));
        beats_t{2} = double(beats_t{2}); HR_t{2} = double(HR_t{2});
        % Bug fix
        for j = 1:2
            if length(HR_t{j}) < length(beats_t{j}); HR_t{j}(end + 1) = HR_t{j}(end); end
            if length(HR_b{j}) < length(beats_b{j}); HR_b{j}(end + 1) = HR_b{j}(end); end
        end

        % -----------------------------------------------------------------
        % Append Dataset by Label
        % -----------------------------------------------------------------

        disp("--> Step 7 of 7: Appending Dataset by Labels")

        % Call trimToIdx function
        trimToIdx_s2(4, obj.levels, include);

        % For each level this subject has, find the beat indices
        allLevelsToBeats(4, include);
        
    end

%% ------------------------------------------------------------------------
% Sub-Function for Subject 5
% -------------------------------------------------------------------------

    function subject5()

        % NOTE: Subject 5 does not have interval Level.all
        % NOTE: The first element in any cell array is from the RELATIVE
        % trial while the second element is from the ABSOLUTE trial.

        % Returns dataset for Subject 5 as a struct:
        % obj.dataset{5}.(interval).ecg_biopac
        % obj.dataset{5}.(interval).ecg_t3
        % obj.dataset{5}.(interval).sternumSCG_x
        % obj.dataset{5}.(interval).sternumSCG_y
        % obj.dataset{5}.(interval).sternumSCG_z
        % obj.dataset{5}.(interval).sternumSCG_z_back
        % obj.dataset{5}.(interval).apexSCG_x
        % obj.dataset{5}.(interval).apexSCG_y
        % obj.dataset{5}.(interval).apexSCG_z
        % obj.dataset{5}.(interval).apexSCG_z_back
        % obj.dataset{5}.(interval).apexSPO2
        % obj.dataset{5}.(interval).femoralSPO2
        % obj.dataset{5}.(interval).apexPPG
        % obj.dataset{5}.(interval).femoralPPG
        % obj.dataset{5}.(interval).aorticPressure
        % obj.dataset{5}.(interval).femoralPressure
        % obj.dataset{5}.(interval).wedgePressure
        % obj.dataset{5}.(interval).rightAtrialPressure

        % -----------------------------------------------------------------
        % Load Dataset and Setup Workspace
        % -----------------------------------------------------------------

        % Import data from T3
        disp("--> Step 1 of 7: Importing T3 Data")
        T3 = importdata('Subject5/t3.txt');

        % Window the data appropriately (into relative and absolute)
        T3_rel = T3(20131986:41551538, :);
        T3_abs = T3(42893327:65069845, :);

        % Clear the T3 workspace to conserve memory
        clear T3

        % Import data from Biopac (relative and absolute trial periods)
        disp("--> Step 2 of 7: Importing Biopac Data")
        biopac_rel = importdata('Subject5/biopac_1.txt');
        biopac_abs = importdata('Subject5/biopac_2.txt');

        % Assign data to variables (T3)
        P1{1} = T3_rel(:, 3); P1{2} = T3_abs(:, 3);     % Aortic arch pressure
        P2{1} = T3_rel(:, 2); P2{2} = T3_abs(:, 2);     % Femoral artery pressure
        P3{1} = T3_rel(:, 4); P3{2} = T3_abs(:, 4);     % Wedge pressure
        P4{1} = T3_rel(:, 5); P4{2} = T3_abs(:, 5);     % Right atrial pressure
        ECG_t{1} = T3_rel(:, 6); ECG_t{2} = T3_abs(:, 6);	% T3 ECG (1)

        % Clear the T3 workspace to conserve memory
        clear T3_rel T3_abs

        % Assign data to variables (Biopac)
        % ECG_b{1} = biopac_rel(:, 1); ECG_b{2} = biopac_abs(:, 1);	% Biopac ECG
        SPO2_1{1} = biopac_rel(:, 2); SPO2_1{2} = biopac_abs(:, 2);	% Apex SPO2
        SPO2_2{1} = biopac_rel(:, 3); SPO2_2{2} = biopac_abs(:, 3);	% Femoral SPO2
        PPG_1{1} = biopac_rel(:, 4); PPG_1{2} = biopac_abs(:, 4);	% Apex PPG
        PPG_2{1} = biopac_rel(:, 5); PPG_2{2} = biopac_abs(:, 5);	% Femoral PPG
        AX_1{1} = biopac_rel(:, 6); AX_1{2} = biopac_abs(:, 6);     % Sternum SCG (x-acceleration)
        AY_1{1} = biopac_rel(:, 7); AY_1{2} = biopac_abs(:, 7);     % Sternum SCG (y-acceleration)
        AZ_1{1} = biopac_rel(:, 8); AZ_1{2} = biopac_abs(:, 8);     % Sternum SCG (z-acceleration)
        AX_2{1} = biopac_rel(:, 9); AX_2{2} = biopac_abs(:, 9);     % Apex SCG (x-acceleration)
        AY_2{1} = biopac_rel(:, 10); AY_2{2} = biopac_abs(:, 10);   % Apex SCG (y-acceleration)
        AZ_2{1} = biopac_rel(:, 11); AZ_2{2} = biopac_abs(:, 11);	% Apex SCG (z-acceleration)

        % Clear the Biopac workspace to conserve memory
        clear biopac_rel biopac_abs

        % -----------------------------------------------------------------
        % Clean T3 Data (NaN Removal)
        % -----------------------------------------------------------------

        disp("--> Step 3 of 7: Cleaning Data")

        % Remove NaNs where necessary
        for k = 1:2
            for j = 1:length(P1{k})
                if isnan(ECG_t{k}(j)); ECG_t{k}(j) = ECG_t{k}(j-1); end
                if isnan(P1{k}(j)); P1{k}(j) = P1{k}(j-1); end
                if isnan(P2{k}(j)); P2{k}(j) = P2{k}(j-1); end
                if isnan(P3{k}(j)); P3{k}(j) = P3{k}(j-1); end
                if isnan(P4{k}(j)); P4{k}(j) = P4{k}(j-1); end
            end
        end

        % -----------------------------------------------------------------
        % Filter Data
        % -----------------------------------------------------------------

        disp("--> Step 4 of 7: Filtering Data")

        if obj.filtered

            % ECG
            % Standard Kaiser window BP filter |5/10-25\40|
            Hd = cardio.general.createBPF(5, 40, obj.Fs,'fpass1', 10, 'fpass2', 25, 'kaiser', 'order', 20);
            for j = 1:2
                % ECG_b{j} = filtfilt(Hd.Numerator, 1, ECG_b{j});
                for k = 1:4; ECG_t{j} = filtfilt(Hd.Numerator, 1, ECG_t{j}); end
            end

            % SCG
            % Standard Kaiser window BP filter |0.5/1-20\40|
            Hd = cardio.general.createBPF(0.5, 40, obj.Fs, 'fpass1', 1,'fpass2', 20, 'kaiser', 'order', 20);
            for j = 1:2
                AX_1{j} = filtfilt(Hd.Numerator, 1, AX_1{j}); AX_2{j} = filtfilt(Hd.Numerator, 1, AX_2{j});
                AY_1{j} = filtfilt(Hd.Numerator, 1, AY_1{j}); AY_2{j} = filtfilt(Hd.Numerator, 1, AY_2{j});
                AZ_1{j} = filtfilt(Hd.Numerator, 1, AZ_1{j}); AZ_2{j} = filtfilt(Hd.Numerator, 1, AZ_2{j});
            end

            % PPG / Arterial Pressures
            % Standard Kaiser window BP filter |0.5/1-9\10|
            Hd = cardio.general.createBPF(0.5, 10, obj.Fs, 'fpass1', 1, 'fpass2', 9, 'kaiser', 'order', 20);
            for j = 1:2
                PPG_1{j} = filtfilt(Hd.Numerator, 1, PPG_1{j}); PPG_2{j} = filtfilt(Hd.Numerator, 1, PPG_2{j});
                P1{j} = filtfilt(Hd.Numerator, 1, P1{j}); P2{j} = filtfilt(Hd.Numerator, 1, P2{j});
                P3{j} = filtfilt(Hd.Numerator, 1, P3{j}); P4{j} = filtfilt(Hd.Numerator, 1, P4{j});
            end

            % PCG
            % Standard Kaiser window HP filter (50Hz)
            Hd = cardio.general.createHPF('Fs', obj.Fs, 'Fc', 50, 'order', 50);
            for j = 1:2; PCG{j} = filtfilt(Hd.Numerator, 1, AZ_1{j}); end

        end

        % -----------------------------------------------------------------
        % Process ECG (BIOPAC)
        % -----------------------------------------------------------------

        disp("--> Step 5 of 7: Processing ECG Data (Biopac)")

        % Set bounds and parameters
        HRV_bound = 2; peak_separation = 840;
        beats_b = cell(2, 1); beats_t = cell(2, 1);
        
        % Import the cleaned ECG data
        ECG_b{1} = importdata('Subject5/biopac_1_clean.mat');
        ECG_b{2} = importdata('Subject5/biopac_2_clean.mat');
        
        % Extract heart rate and beat-separated ECG (Relative)
        [HR_b{1}, beats_b{1}] = cardio.ecg.cleanECG(ECG_b{1}, obj.Fs, ...
            HRV_bound, 'peak_sep', peak_separation, 'double', 'plotsoff');

        % -----------------------------------------------------------------
        % Process ECG (T3)
        % -----------------------------------------------------------------

        disp("--> Step 6 of 7: Processing ECG Data (T3)")
        
        % Import the cleaned ECG data
        ECG_t{1} = importdata('Subject5/t3_clean.mat');
        ECG_t{2} = diff(ECG_t{2});
        
        % Extract heart rate and beat-separated ECG
        for j = 1:2
            [HR_t{j}, beats_t{j}] = cardio.ecg.cleanECG(ECG_t{j}, obj.Fs, ...
                HRV_bound, 'peak_sep', peak_separation, 'double', 'plotsoff');
        end
        
        % Use the T3 data to correct the Biopac data (Absolute)
        [~, peakLocs] = findpeaks(ECG_b{2}, 'MinPeakHeight', 0.005);
        beats_b{2} = []; offset = 10; % Define offset and initialize beats_b
        for j = 1:length(beats_t{2})
            [~, idx] = min(abs(peakLocs - offset - beats_t{2}(j))); % T lags B
            beats_b{2} = [beats_b{2} peakLocs(idx)]; offset = beats_b{2}(j) - beats_t{2}(j);
        end; HR_b{2} = obj.Fs*60./diff(beats_b{2});

        % Because you cut out more noise at the end of T3, shorten the bio
        % heartbeats off the end so you're working with the same number of beats: 

        % Relative hypovolemia
        beats_t{1} = double(beats_t{1}); HR_t{1} = double(HR_t{1});
        HR_b{1} = double(HR_b{1});  beats_b{1} = double(beats_b{1});

        % Absolute hypovolemia
        beats_b{2} = double(beats_b{2}(1:length(beats_t{2})));
        HR_b{2} = double(HR_b{2}(1:length(HR_t{2})));
        beats_t{2} = double(beats_t{2}); HR_t{2} = double(HR_t{2});
        % Bug fix
        for j = 1:2
            if length(HR_t{j}) < length(beats_t{j}); HR_t{j}(end + 1) = HR_t{j}(end); end
            if length(HR_b{j}) < length(beats_b{j}); HR_b{j}(end + 1) = HR_b{j}(end); end
        end

        % -----------------------------------------------------------------
        % Append Dataset by Label
        % -----------------------------------------------------------------

        disp("--> Step 7 of 7: Appending Dataset by Labels")

        % Call trimToIdx function
        trimToIdx_s2(5, obj.levels, include);
        
        % For each level this subject has, find the beat indices
        allLevelsToBeats(5, include);

    end

%% ------------------------------------------------------------------------
% Sub-Function for Subject 6
% -------------------------------------------------------------------------

    function subject6()

        % NOTE: Subject 6 does not have interval Level.all
        % NOTE: The first element in any cell array is from the RELATIVE
        % trial while the second element is from the ABSOLUTE trial.

        % Returns dataset for Subject 6 as a struct:
        % obj.dataset{6}.(interval).ecg_biopac
        % obj.dataset{6}.(interval).ecg_t3
        % obj.dataset{6}.(interval).sternumSCG_x
        % obj.dataset{6}.(interval).sternumSCG_y
        % obj.dataset{6}.(interval).sternumSCG_z
        % obj.dataset{6}.(interval).sternumSCG_z_back
        % obj.dataset{6}.(interval).apexSCG_x
        % obj.dataset{6}.(interval).apexSCG_y
        % obj.dataset{6}.(interval).apexSCG_z
        % obj.dataset{6}.(interval).apexSCG_z_back
        % obj.dataset{6}.(interval).apexSPO2
        % obj.dataset{6}.(interval).femoralSPO2
        % obj.dataset{6}.(interval).apexPPG
        % obj.dataset{6}.(interval).femoralPPG
        % obj.dataset{6}.(interval).aorticPressure
        % obj.dataset{6}.(interval).femoralPressure
        % obj.dataset{6}.(interval).wedgePressure
        % obj.dataset{6}.(interval).rightAtrialPressure

        % -----------------------------------------------------------------
        % Load Dataset and Setup Workspace
        % -----------------------------------------------------------------

        % Import data from T3
        disp("--> Step 1 of 7: Importing T3 Data")
        T3 = importdata('Subject6/t3.txt'); T3 = T3.data;

        % Window the data appropriately (into relative and absolute)
        T3_rel = T3(15878950:22543798, :);
        T3_abs = T3(23055559:47036665, :);

        % Clear the T3 workspace to conserve memory
        clear T3

        % Import data from Biopac (relative and absolute trial periods)
        disp("--> Step 2 of 7: Importing Biopac Data")
        biopac_rel = importdata('Subject6/biopac_1.txt');
        biopac_abs = importdata('Subject6/biopac_2.txt');

        % Assign data to variables (T3)
        P1{1} = T3_rel(:, 2); P1{2} = T3_abs(:, 2);     % Aortic arch pressure
        P2{1} = T3_rel(:, 1); P2{2} = T3_abs(:, 1);     % Femoral artery pressure (switched w/CH 1)
        P3{1} = T3_rel(:, 3); P3{2} = T3_abs(:, 3);     % Wedge pressure
        P4{1} = T3_rel(:, 4); P4{2} = T3_abs(:, 4);     % Right atrial pressure
        ECG_t{1} = T3_rel(:, 6); ECG_t{2} = T3_abs(:, 6);	% T3 ECG (1)

        % Clear the T3 workspace to conserve memory
        clear T3_rel T3_abs

        % Assign data to variables (Biopac)
        ECG_b{1} = biopac_rel(:, 1); ECG_b{2} = biopac_abs(:, 1);	% Biopac ECG
        SPO2_1{1} = biopac_rel(:, 2); SPO2_1{2} = biopac_abs(:, 2);	% Apex SPO2
        SPO2_2{1} = biopac_rel(:, 3); SPO2_2{2} = biopac_abs(:, 3);	% Femoral SPO2
        PPG_1{1} = biopac_rel(:, 4); PPG_1{2} = biopac_abs(:, 4);	% Apex PPG
        PPG_2{1} = biopac_rel(:, 5); PPG_2{2} = biopac_abs(:, 5);	% Femoral PPG
        AX_1{1} = biopac_rel(:, 6); AX_1{2} = biopac_abs(:, 6);     % Sternum SCG (x-acceleration)
        AY_1{1} = biopac_rel(:, 7); AY_1{2} = biopac_abs(:, 7);     % Sternum SCG (y-acceleration)
        AZ_1{1} = biopac_rel(:, 8); AZ_1{2} = biopac_abs(:, 8);     % Sternum SCG (z-acceleration)
        AX_2{1} = biopac_rel(:, 9); AX_2{2} = biopac_abs(:, 9);     % Apex SCG (x-acceleration)
        AY_2{1} = biopac_rel(:, 10); AY_2{2} = biopac_abs(:, 10);   % Apex SCG (y-acceleration)
        AZ_2{1} = biopac_rel(:, 11); AZ_2{2} = biopac_abs(:, 11);	% Apex SCG (z-acceleration)

        % Clear the Biopac workspace to conserve memory
        clear biopac_rel biopac_abs

        % -----------------------------------------------------------------
        % Clean T3 Data (NaN Removal)
        % -----------------------------------------------------------------

        disp("--> Step 3 of 7: Cleaning Data")

        % Remove NaNs where necessary
        for k = 1:2
            for j = 1:length(ECG_t{k})
                if isnan(ECG_t{k}(j)); ECG_t{k}(j) = ECG_t{k}(j-1); end
                if isnan(P1{k}(j)); P1{k}(j) = P1{k}(j-1); end
                if isnan(P2{k}(j)); P2{k}(j) = P2{k}(j-1); end
                if isnan(P3{k}(j)); P3{k}(j) = P3{k}(j-1); end
                if isnan(P4{k}(j)); P4{k}(j) = P4{k}(j-1); end
            end
        end

        % -----------------------------------------------------------------
        % Filter Data
        % -----------------------------------------------------------------

        disp("--> Step 4 of 7: Filtering Data")

        if obj.filtered

            % ECG
            % Standard Kaiser window BP filter |5/10-25\40|
            Hd = cardio.general.createBPF(5, 40, obj.Fs,'fpass1', 10, 'fpass2', 25, 'kaiser', 'order', 20);
            for j = 1:2
                ECG_b{j} = filtfilt(Hd.Numerator, 1, ECG_b{j});
                for k = 1:4; ECG_t{j} = filtfilt(Hd.Numerator, 1, ECG_t{j}); end
            end

            % SCG
            % Standard Kaiser window BP filter |0.5/1-20\40|
            Hd = cardio.general.createBPF(0.5, 40, obj.Fs, 'fpass1', 1,'fpass2', 20, 'kaiser', 'order', 20);
            for j = 1:2
                AX_1{j} = filtfilt(Hd.Numerator, 1, AX_1{j}); AX_2{j} = filtfilt(Hd.Numerator, 1, AX_2{j});
                AY_1{j} = filtfilt(Hd.Numerator, 1, AY_1{j}); AY_2{j} = filtfilt(Hd.Numerator, 1, AY_2{j});
                AZ_1{j} = filtfilt(Hd.Numerator, 1, AZ_1{j}); AZ_2{j} = filtfilt(Hd.Numerator, 1, AZ_2{j});
            end

            % PPG / Arterial Pressures
            % Standard Kaiser window BP filter |0.5/1-9\10|
            Hd = cardio.general.createBPF(0.5, 10, obj.Fs, 'fpass1', 1, 'fpass2', 9, 'kaiser', 'order', 20);
            for j = 1:2
                PPG_1{j} = filtfilt(Hd.Numerator, 1, PPG_1{j}); PPG_2{j} = filtfilt(Hd.Numerator, 1, PPG_2{j});
                P1{j} = filtfilt(Hd.Numerator, 1, P1{j}); P2{j} = filtfilt(Hd.Numerator, 1, P2{j});
                P3{j} = filtfilt(Hd.Numerator, 1, P3{j}); P4{j} = filtfilt(Hd.Numerator, 1, P4{j});
            end

            % PCG
            % Standard Kaiser window HP filter (50Hz)
            Hd = cardio.general.createHPF('Fs', obj.Fs, 'Fc', 50, 'order', 50);
            for j = 1:2; PCG{j} = filtfilt(Hd.Numerator, 1, AZ_1{j}); end

        end

        % -----------------------------------------------------------------
        % Process ECG (BIOPAC)
        % -----------------------------------------------------------------

        disp("--> Step 5 of 7: Processing ECG Data (Biopac)")

        % 1. Check polarity. (Are the R-peaks upside down?) 
        % 2. Take the difference, check its shape.
        % 3. Append a 0 as the last element so it is the same length as the
        % other signals. 

        % Differentiate ECG signal
        diff1_ecgb = cell(2, 1); for j = 1:2; diff1_ecgb{j} = diff(ECG_b{j}); end

        % Append a 0 so it's the same length as everything else
        ECG_bio = cell(2, 1); for j = 1:2; ECG_bio{j} = [diff1_ecgb{j}; zeros(1)]; end
        clear diff1_ecgb

        % Fix outliers so ECG_bio can be processed by cardio.ecg.automatedHR
        ECG_bio_trim = cell(2, 1); limit = 0.012;
        for j = 1:2; ECG_bio_trim{j} = zeros(size(ECG_bio{j})); end
        for k = 1:2
            for j = 1:length(ECG_bio{k})
               if(abs(ECG_bio{k}(j)) > limit)
                   % Keep the relative distance so you don't lose the peaks
                   ECG_bio_trim{k}(j) = (0.01*(abs(ECG_bio{k}(j))-limit)+limit)*sign(ECG_bio{k}(j)); 
               else; ECG_bio_trim{k}(j) = ECG_bio{k}(j); 
               end
            end
        end

        % Clear extra workspace variables
        clear ECG_bio

        % Set parameters for further further cleaning and beat-separation
        peak_separation = 840; ECG_bt_2 = cell(2, 1); HRV_bound = 1;
        for j = 1:2; ECG_bt_2{j} = ECG_bio_trim{j}; end; clear ECG_bio_trim

        % Correct outliers (Relative)
        ECG_bt_2{1}(146000:146600) = ECG_bt_2{1}(146000:146600)*0.1;
        
        % Fix errors
        for j = 1:2; ECG_bt_2{j} = fix_diff(ECG_bt_2{j}); end

        % Extract heart rate and beat-separated ECG
        beats_b = cell(2, 1);
        for j = 1:2
            [HR_b{j}, beats_b{j}] = cardio.ecg.cleanECG(ECG_bt_2{j}, obj.Fs, ...
                HRV_bound, 'peak_sep', peak_separation, 'plotsoff');
        end

        % -----------------------------------------------------------------
        % Process ECG (T3)
        % -----------------------------------------------------------------

        disp("--> Step 6 of 7: Processing ECG Data (T3)")

        % Set bounds and variables
        bound_2 = 0.01; ECG_t3_trim = cell(2, 1);
        for j = 1:2; ECG_t3_trim{j} = diff(ECG_t{j}); end
        
        % Correct errors in T3 ECG
        for j = 1:length(ECG_t3_trim{2})
            if ECG_t3_trim{2}(j) > bound_2
                ECG_t3_trim{2}(j) = ECG_t3_trim{2}(j)*0.01 + bound_2;
            end
        end
        
        % Rename workspace variable
        ECG_tt_2 = cell(2, 1); for j = 1:2; ECG_tt_2{j} = fix_diff(ECG_t3_trim{j}); end
        clear ECG_t3_trim

        % Correct for outliers (Relative)
        ECG_tt_2{1}(2937400:2937700) = ECG_tt_2{1}(2937400:2937700)*0.1;
        ECG_tt_2{1}(2937920:2938000) = ECG_tt_2{1}(2937920:2938000) + 0.15;
        ECG_tt_2{1}(6101720:6101750) = ECG_tt_2{1}(6101720:6101750) + 0.15;

        % Extract heart rate and beat-separated ECG (Relative)
        beats_t = cell(2, 1);
        [HR_t{1}, beats_t{1}] = cardio.ecg.cleanECG(ECG_tt_2{1}, obj.Fs, ...
                HRV_bound, 'peak_sep', peak_separation, 'plotsoff');
        
        % Extract heart rate and beat-separated ECG (Absolute)
        [~, peakLocs] = findpeaks(ECG_tt_2{2}, 'MinPeakHeight', 0.0005);
        tPeaks = []; offset = 0;
        for j = 1:length(beats_b{2})
            [~, idx] = min(abs(beats_b{2}(j) - offset - peakLocs));
            tPeaks = [tPeaks peakLocs(idx)]; offset = beats_b{2}(j) - tPeaks(j);
        end; HR_t{2} = obj.Fs*60./diff(tPeaks); beats_t{2} = tPeaks;
            
        % Because you cut out more noise at the end of T3, shorten the bio
        % heartbeats off the end so you're working with the same number of beats: 

        % Relative hypovolemia
        beats_t{1} = double(beats_t{1}); HR_t{1} = double(HR_t{1});
        HR_b{1} = double(HR_b{1});  beats_b{1} = double(beats_b{1});

        % Absolute hypovolemia
        beats_b{2} = double(beats_b{2}(1:length(beats_t{2})));
        HR_b{2} = double(HR_b{2}(1:length(HR_t{2})));
        beats_t{2} = double(beats_t{2}); HR_t{2} = double(HR_t{2});
        % Bug fix
        for j = 1:2
            if length(HR_t{j}) < length(beats_t{j}); HR_t{j}(end + 1) = HR_t{j}(end); end
            if length(HR_b{j}) < length(beats_b{j}); HR_b{j}(end + 1) = HR_b{j}(end); end
        end

        % -----------------------------------------------------------------
        % Append Dataset by Label
        % -----------------------------------------------------------------

        disp("--> Step 7 of 7: Appending Dataset by Labels")

        % Call trimToIdx function
        trimToIdx_s2(6, obj.levels, include);

        % For each level this subject has, find the beat indices
        allLevelsToBeats(6, include);
        
    end

%% ------------------------------------------------------------------------
% Helper Functions
% -------------------------------------------------------------------------

% Trim all data modalities to index
    function trimToIdx_s1(levels, include)
        
        % If the data is to be beat-separated, separate the data beatwise
        % before saving the valid levels
        if obj.beatSeparated
        
            % If the data is not to be truncated, set the signal length to
            % the maximum RR interval, and for each resulting signal
            % nan-pad the end
            if truncation
                len = sigLen;
            else
                len_b = max(diff(beats_b));
                len_t = max(diff(beats_t));
                len = max([len_b, len_t]);
            end
            
            % Beat-separate the data
            ecg_biopac = cardio.general.separateBeats(ECG_b, 'indices', beats_b, 'samples', len, 'nanpad', ~truncation);
            ecg_t3 = cardio.general.separateBeats(ECG_t, 'indices', beats_t, 'samples', len, 'nanpad', ~truncation);
            sternumSCG_x = cardio.general.separateBeats(AX_1, 'indices', beats_b, 'samples', len, 'nanpad', ~truncation);
            sternumSCG_y = cardio.general.separateBeats(AY_1, 'indices', beats_b, 'samples', len, 'nanpad', ~truncation);
            sternumSCG_z = cardio.general.separateBeats(AZ_1, 'indices', beats_b, 'samples', len, 'nanpad', ~truncation);
            sternumSCG_z_back = cardio.general.separateBeats(AZ_1, 'indices', beats_b, 'samples', sigLen, 'backward');
            apexSCG_x = cardio.general.separateBeats(AX_2, 'indices', beats_b, 'samples', len, 'nanpad', ~truncation);
            apexSCG_y = cardio.general.separateBeats(AY_2, 'indices', beats_b, 'samples', len, 'nanpad', ~truncation);
            apexSCG_z = cardio.general.separateBeats(AZ_2, 'indices', beats_b, 'samples', len, 'nanpad', ~truncation);
            apexSCG_z_back = cardio.general.separateBeats(AZ_2, 'indices', beats_b, 'samples', sigLen, 'backward');
            apexSPO2 = cardio.general.separateBeats(SPO2_1, 'indices', beats_b, 'samples', len, 'nanpad', ~truncation);
            femoralSPO2 = cardio.general.separateBeats(SPO2_2, 'indices', beats_b, 'samples', len, 'nanpad', ~truncation);
            apexPPG = cardio.general.separateBeats(PPG_1, 'indices', beats_b, 'samples', len, 'nanpad', ~truncation);
            femoralPPG = cardio.general.separateBeats(PPG_2, 'indices', beats_b, 'samples', len, 'nanpad', ~truncation);
            aorticPressure = cardio.general.separateBeats(P1, 'indices', beats_t, 'samples', len, 'nanpad', ~truncation);
            aorticPressure_back = cardio.general.separateBeats(P1, 'indices', beats_t, 'samples', sigLen, 'backward');
            femoralPressure = cardio.general.separateBeats(P2, 'indices', beats_t, 'samples', len, 'nanpad', ~truncation);
            wedgePressure = cardio.general.separateBeats(P3, 'indices', beats_t, 'samples', len, 'nanpad', ~truncation);
            rightAtrialPressure = cardio.general.separateBeats(P4, 'indices', beats_t, 'samples', len, 'nanpad', ~truncation);
            sternumPCG = cardio.general.separateBeats(PCG, 'indices', beats_b, 'samples', len, 'nanpad', ~truncation);
            sternumPCG_back = cardio.general.separateBeats(PCG, 'indices', beats_b, 'samples', sigLen, 'backward');
            
            % For each level...
            for l = 1:length(levels)
                
                % Extract the current level
                level = levels(l);
                
                % Get the beat indices for the current level
                beatIdx = levelToBeats(level, 1, include); 
                if level == Level.allRelative || level == Level.relBaseline2; beatIdx = beatIdx(1):11251; end
                if isnan(beatIdx(1)); continue; end     % Continue if the level does not exist
                
                % Return beats_b and beats_t
                obj.dataset{counter}.(string(level)).beats_biopac = beats_b(beatIdx);
                obj.dataset{counter}.(string(level)).beats_t3 = beats_t(beatIdx);
                
                % Return the beats in range for each modality
                obj.dataset{counter}.(string(level)).ecg_biopac = ecg_biopac(:, beatIdx);
                obj.dataset{counter}.(string(level)).ecg_t3 = ecg_t3(:, beatIdx);
                obj.dataset{counter}.(string(level)).sternumSCG_x = sternumSCG_x(:, beatIdx);
                obj.dataset{counter}.(string(level)).sternumSCG_y = sternumSCG_y(:, beatIdx);
                obj.dataset{counter}.(string(level)).sternumSCG_z = sternumSCG_z(:, beatIdx);
                obj.dataset{counter}.(string(level)).sternumSCG_z_back = sternumSCG_z_back(:, beatIdx);
                obj.dataset{counter}.(string(level)).apexSCG_x = apexSCG_x(:, beatIdx);
                obj.dataset{counter}.(string(level)).apexSCG_y = apexSCG_y(:, beatIdx);
                obj.dataset{counter}.(string(level)).apexSCG_z = apexSCG_z(:, beatIdx);
                obj.dataset{counter}.(string(level)).apexSCG_z_back = apexSCG_z_back(:, beatIdx);
                obj.dataset{counter}.(string(level)).apexSPO2 = apexSPO2(:, beatIdx);
                obj.dataset{counter}.(string(level)).femoralSPO2 = femoralSPO2(:, beatIdx);
                obj.dataset{counter}.(string(level)).apexPPG = apexPPG(:, beatIdx);
                obj.dataset{counter}.(string(level)).femoralPPG = femoralPPG(:, beatIdx);
                obj.dataset{counter}.(string(level)).aorticPressure = aorticPressure(:, beatIdx);
                obj.dataset{counter}.(string(level)).aorticPressure_back = aorticPressure_back(:, beatIdx);
                obj.dataset{counter}.(string(level)).femoralPressure = femoralPressure(:, beatIdx);
                obj.dataset{counter}.(string(level)).wedgePressure = wedgePressure(:, beatIdx);
                obj.dataset{counter}.(string(level)).rightAtrialPressure = rightAtrialPressure(:, beatIdx);
                obj.dataset{counter}.(string(level)).HR_biopac = HR_b(beatIdx);
                obj.dataset{counter}.(string(level)).HR_t3 = HR_t(beatIdx);
                obj.dataset{counter}.(string(level)).sternumPCG = sternumPCG(:, beatIdx);
                obj.dataset{counter}.(string(level)).sternumPCG_back = sternumPCG_back(:, beatIdx);
                
                % Return the un-segmented ECG signal
                obj.dataset{counter}.(string(level)).ecg_biopac_full = ECG_b(beats_b(beatIdx(1)):beats_b(beatIdx(end)));
                obj.dataset{counter}.(string(level)).ecg_t3_full = ECG_t(beats_b(beatIdx(1)):beats_b(beatIdx(end)));
                
            end
                
        end
        
        % If the data is not beat-separated, convert the beat indices to
        % samples for the current level
        if ~obj.beatSeparated
            
            % For each level...
            for l = 1:length(levels)
                
                % Extract the current level
                level = levels(l);
                
                % Get the beat indices for the current level
                beatIdx = levelToBeats(level, 1, include);
                if level == Level.allRelative || level == Level.relBaseline2; beatIdx(2) = 11251; end
                if isnan(beatIdx(1)); continue; end     % Continue if the level does not exist
                
                % Return the signal vector for each modality
                obj.dataset{counter}.(string(level)).ecg_biopac = ECG_b(beats_b(beatIdx(1)):beats_b(beatIdx(end)));
                obj.dataset{counter}.(string(level)).ecg_t3 = ECG_t(beats_t(beatIdx(1)):beats_t(beatIdx(end)));
                obj.dataset{counter}.(string(level)).sternumSCG_x = AX_1(beats_b(beatIdx(1)):beats_b(beatIdx(end)));
                obj.dataset{counter}.(string(level)).sternumSCG_y = AY_1(beats_b(beatIdx(1)):beats_b(beatIdx(end)));
                obj.dataset{counter}.(string(level)).sternumSCG_z = AZ_1(beats_b(beatIdx(1)):beats_b(beatIdx(end)));
                obj.dataset{counter}.(string(level)).apexSCG_x = AX_2(beats_b(beatIdx(1)):beats_b(beatIdx(end)));
                obj.dataset{counter}.(string(level)).apexSCG_y = AY_2(beats_b(beatIdx(1)):beats_b(beatIdx(end)));
                obj.dataset{counter}.(string(level)).apexSCG_z = AZ_2(beats_b(beatIdx(1)):beats_b(beatIdx(end)));
                obj.dataset{counter}.(string(level)).apexSPO2 = SPO2_1(beats_b(beatIdx(1)):beats_b(beatIdx(end)));
                obj.dataset{counter}.(string(level)).femoralSPO2 = SPO2_2(beats_b(beatIdx(1)):beats_b(beatIdx(end)));
                obj.dataset{counter}.(string(level)).apexPPG = PPG_1(beats_b(beatIdx(1)):beats_b(beatIdx(end)));
                obj.dataset{counter}.(string(level)).femoralPPG = PPG_2(beats_b(beatIdx(1)):beats_b(beatIdx(end)));
                obj.dataset{counter}.(string(level)).aorticPressure = P1(beats_t(beatIdx(1)):beats_t(beatIdx(end)));
                obj.dataset{counter}.(string(level)).femoralPressure = P2(beats_t(beatIdx(1)):beats_t(beatIdx(end)));
                obj.dataset{counter}.(string(level)).wedgePressure = P3(beats_t(beatIdx(1)):beats_t(beatIdx(end)));
                obj.dataset{counter}.(string(level)).rightAtrialPressure = P4(beats_t(beatIdx(1)):beats_t(beatIdx(end)));
                obj.dataset{counter}.(string(level)).sternumPCG = PCG(beats_b(beatIdx(1)):beats_b(beatIdx(end)));
                
            end
            
        end
        
    end

% Trim all data modalities to index (Subject 2+)
    function trimToIdx_s2(subject, levels, include)
        
        % Perform the following for (1) relative and (2) absolute hypovolemia
        for set = 1:2
            
            % If the data is to be beat-separated, separate the data beatwise
            % before saving the valid levels
            if obj.beatSeparated
                
                % If the data is not to be truncated, set the signal length to
                % the maximum RR interval, and for each resulting signal
                % nan-pad the end
                if truncation
                    len = sigLen;
                else
                    len_b = max(diff(beats_b{set}));
                    len_t = max(diff(beats_t{set}));
                    len = max([len_b, len_t]);
                end
                
                % Beat-separate the data
                ecg_biopac = cardio.general.separateBeats(ECG_b{set}, 'indices', beats_b{set}, 'samples', len, 'nanpad', ~truncation);
                ecg_t3 = cardio.general.separateBeats(ECG_t{set}, 'indices', beats_t{set}, 'samples', len, 'nanpad', ~truncation);
                sternumSCG_x = cardio.general.separateBeats(AX_1{set}, 'indices', beats_b{set}, 'samples', len, 'nanpad', ~truncation);
                sternumSCG_y = cardio.general.separateBeats(AY_1{set}, 'indices', beats_b{set}, 'samples', len, 'nanpad', ~truncation);
                sternumSCG_z = cardio.general.separateBeats(AZ_1{set}, 'indices', beats_b{set}, 'samples', len, 'nanpad', ~truncation);
                try sternumSCG_z_back = cardio.general.separateBeats(AZ_1{set}, 'indices', beats_b{set}, 'samples', sigLen, 'backward'); catch; end
                apexSCG_x = cardio.general.separateBeats(AX_2{set}, 'indices', beats_b{set}, 'samples', len, 'nanpad', ~truncation);
                apexSCG_y = cardio.general.separateBeats(AY_2{set}, 'indices', beats_b{set}, 'samples', len, 'nanpad', ~truncation);
                apexSCG_z = cardio.general.separateBeats(AZ_2{set}, 'indices', beats_b{set}, 'samples', len, 'nanpad', ~truncation);
                try apexSCG_z_back = cardio.general.separateBeats(AZ_2{set}, 'indices', beats_b{set}, 'samples', sigLen, 'backward'); catch; end
                apexSPO2 = cardio.general.separateBeats(SPO2_1{set}, 'indices', beats_b{set}, 'samples', len, 'nanpad', ~truncation);
                femoralSPO2 = cardio.general.separateBeats(SPO2_2{set}, 'indices', beats_b{set}, 'samples', len, 'nanpad', ~truncation);
                apexPPG = cardio.general.separateBeats(PPG_1{set}, 'indices', beats_b{set}, 'samples', len, 'nanpad', ~truncation);
                femoralPPG = cardio.general.separateBeats(PPG_2{set}, 'indices', beats_b{set}, 'samples', len, 'nanpad', ~truncation);
                aorticPressure = cardio.general.separateBeats(P1{set}, 'indices', beats_t{set}, 'samples', len, 'nanpad', ~truncation);
                try aorticPressure_back = cardio.general.separateBeats(P1{set}, 'indices', beats_t{set}, 'samples', sigLen, 'backward'); catch; end
                femoralPressure = cardio.general.separateBeats(P2{set}, 'indices', beats_t{set}, 'samples', len, 'nanpad', ~truncation);
                wedgePressure = cardio.general.separateBeats(P3{set}, 'indices', beats_t{set}, 'samples', len, 'nanpad', ~truncation);
                rightAtrialPressure = cardio.general.separateBeats(P4{set}, 'indices', beats_t{set}, 'samples', len, 'nanpad', ~truncation);
                sternumPCG = cardio.general.separateBeats(PCG{set}, 'indices', beats_b{set}, 'samples', len, 'nanpad', ~truncation);
                try sternumPCG_back = cardio.general.separateBeats(PCG{set}, 'indices', beats_b{set}, 'samples', sigLen, 'backward'); catch; end
                
                % For each level...
                for l = 1:length(levels)
                    
                    % Extract the current level
                    level = levels(l);
                    
                    % For the current set, run only the appropriate levels
                    if set == 1 && obj.getParent(level) == Level.allAbsolute; continue; end
                    if set == 2 && obj.getParent(level) == Level.allRelative; continue; end
                    
                    % Get the beat indices for the current level and save
                    beatIdx = levelToBeats(level, subject, include);
                    if isnan(beatIdx(1)); continue; end     % Continue if the level does not exist
                    if level == Level.allRelative || level == Level.relBaseline2; beatIdx = beatIdx(1):length(beats_b{set}); end
                    if subject == 6 && set == 1; beatIdx = beatIdx(1:end-1); end
                    obj.dataset{counter}.(string(level)).indices = beatIdx;
                    
                    % Return beats_b and beats_t
                    obj.dataset{counter}.(string(level)).beats_biopac = beats_b{set}(beatIdx);
                    obj.dataset{counter}.(string(level)).beats_t3 = beats_t{set}(beatIdx);
                    
                    % Return the beats in range for each modality
                    obj.dataset{counter}.(string(level)).ecg_biopac = ecg_biopac(:, beatIdx);
                    obj.dataset{counter}.(string(level)).ecg_t3 = ecg_t3(:, beatIdx);
                    obj.dataset{counter}.(string(level)).sternumSCG_x = sternumSCG_x(:, beatIdx);
                    obj.dataset{counter}.(string(level)).sternumSCG_y = sternumSCG_y(:, beatIdx);
                    obj.dataset{counter}.(string(level)).sternumSCG_z = sternumSCG_z(:, beatIdx);
                    try obj.dataset{counter}.(string(level)).sternumSCG_z_back = sternumSCG_z_back(:, beatIdx); catch; end
                    obj.dataset{counter}.(string(level)).apexSCG_x = apexSCG_x(:, beatIdx);
                    obj.dataset{counter}.(string(level)).apexSCG_y = apexSCG_y(:, beatIdx);
                    obj.dataset{counter}.(string(level)).apexSCG_z = apexSCG_z(:, beatIdx);
                    try obj.dataset{counter}.(string(level)).apexSCG_z_back = apexSCG_z_back(:, beatIdx); catch; end
                    obj.dataset{counter}.(string(level)).apexSPO2 = apexSPO2(:, beatIdx);
                    obj.dataset{counter}.(string(level)).femoralSPO2 = femoralSPO2(:, beatIdx);
                    obj.dataset{counter}.(string(level)).apexPPG = apexPPG(:, beatIdx);
                    obj.dataset{counter}.(string(level)).femoralPPG = femoralPPG(:, beatIdx);
                    obj.dataset{counter}.(string(level)).aorticPressure = aorticPressure(:, beatIdx);
                    try obj.dataset{counter}.(string(level)).aorticPressure_back = aorticPressure_back(:, beatIdx); catch; end
                    obj.dataset{counter}.(string(level)).femoralPressure = femoralPressure(:, beatIdx);
                    obj.dataset{counter}.(string(level)).wedgePressure = wedgePressure(:, beatIdx);
                    obj.dataset{counter}.(string(level)).rightAtrialPressure = rightAtrialPressure(:, beatIdx);
                    obj.dataset{counter}.(string(level)).HR_biopac = HR_b{set}(beatIdx);
                    obj.dataset{counter}.(string(level)).HR_t3 = HR_t{set}(beatIdx);
                    obj.dataset{counter}.(string(level)).sternumPCG = sternumPCG(:, beatIdx);
                    try obj.dataset{counter}.(string(level)).sternumPCG_back = sternumPCG_back(:, beatIdx); catch; end
                    
                    % Return the un-segmented ECG signal
                    obj.dataset{counter}.(string(level)).ecg_biopac_full = ECG_b{set}(beats_b{set}(beatIdx(1)):beats_b{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).ecg_t3_full = ECG_t{set}(beats_t{set}(beatIdx(1)):beats_t{set}(beatIdx(end)));
                    
                end
                
            end
            
            % If the data is not beat-separated, convert the beat indices to
            % samples for the current level
            if ~obj.beatSeparated

                % For each level...
                for l = 1:length(levels)

                    % Extract the current level
                    level = levels(l);
                    
                    % For the current set, run only the appropriate levels
                    if set == 1 && obj.getParent(level) == Level.allAbsolute; continue; end
                    if set == 2 && obj.getParent(level) == Level.allRelative; continue; end

                    % Get the beat indices for the current level
                    beatIdx = levelToBeats(level, subject, include);
                    if level == Level.allRelative || level == Level.relBaseline2; beatIdx = beatIdx(1):length(beats_b{set}); end
                    if isnan(beatIdx(1)); continue; end     % Continue if the level does not exist

                    % Return the signal vector for each modality
                    obj.dataset{counter}.(string(level)).ecg_biopac = ECG_b{set}(beats_b{set}(beatIdx(1)):beats_b{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).ecg_t3 = ECG_t{set}(beats_t{set}(beatIdx(1)):beats_t{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).sternumSCG_x = AX_1{set}(beats_b{set}(beatIdx(1)):beats_b{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).sternumSCG_y = AY_1{set}(beats_b{set}(beatIdx(1)):beats_b{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).sternumSCG_z = AZ_1{set}(beats_b{set}(beatIdx(1)):beats_b{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).apexSCG_x = AX_2{set}(beats_b{set}(beatIdx(1)):beats_b{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).apexSCG_y = AY_2{set}(beats_b{set}(beatIdx(1)):beats_b{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).apexSCG_z = AZ_2{set}(beats_b{set}(beatIdx(1)):beats_b{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).apexSPO2 = SPO2_1{set}(beats_b{set}(beatIdx(1)):beats_b{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).femoralSPO2 = SPO2_2{set}(beats_b{set}(beatIdx(1)):beats_b{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).apexPPG = PPG_1{set}(beats_b{set}(beatIdx(1)):beats_b{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).femoralPPG = PPG_2{set}(beats_b{set}(beatIdx(1)):beats_b{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).aorticPressure = P1{set}(beats_t{set}(beatIdx(1)):beats_t{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).femoralPressure = P2{set}(beats_t{set}(beatIdx(1)):beats_t{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).wedgePressure = P3{set}(beats_t{set}(beatIdx(1)):beats_t{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).rightAtrialPressure = P4{set}(beats_t{set}(beatIdx(1)):beats_t{set}(beatIdx(end)));
                    obj.dataset{counter}.(string(level)).sternumPCG = PCG{set}(beats_b{set}(beatIdx(1)):beats_b{set}(beatIdx(end)));

                end

            end
            
        end
        
    end

% Enable the cleanECG function to choose consistent peaks when the diff
% function has been used to remove noise from the ECG.
    function new_ecg = fix_diff(ecg)
        
        % Only look at peaks that are close to the max - fix the majority
        max_bound = max(ecg); threshold = 0.6*max_bound; new_ecg = ecg;
        [~, peakLocs] = findpeaks(ecg, 'MinPeakHeight', threshold, 'MinPeakDistance', 30);
        
        % For each peak...
        for peak = 2:length(peakLocs)
            if peakLocs(peak) - peakLocs(peak - 1) < 500
                for idx = (peakLocs(peak - 1) + 20):(peakLocs(peak) + 40)
                    new_ecg(idx) = new_ecg(idx)*0.1; % Remove extra peak
                end
            end
        end
        
    end

% Extract beat indices for the specified level for each subject
    function beatIdx = levelToBeats(level, subject, include)
        
        % Load the event indices .mat file
        filename = strcat(obj.path, 'eventIndices.mat'); load(filename, 'eventIndices');
        tab = eventIndices; clear eventIndices;
        if subject == 1
            % For subject 1, add 11250 beats to absBaseline1 and beyond (one
            % file for both levels)
            startIdx = find(tab.event == Event.absBaseline1);
            tab.subject1(startIdx:end) = tab.subject1(startIdx:end) + 11250;
        end
        
        % Convert subject to string
        subject = "subject" + string(subject);
        
        % "include" is a BOOL indicating whether to include bleeding
        % regions during absolute hypovolemia
        
        % Extract subject's data
        idx = tab.(subject);
        
        if ~include
            
            % For the specified level...
            switch level

                case Level.all
                    if strcmp(subject, "subject1")
                        beatIdx = idx(tab.event == Event.relBaseline1):idx(tab.event == Event.startMortality);
                    else; beatIdx = nan;
                    end
                case Level.allRelative
                    beatIdx = idx(tab.event == Event.relBaseline1);
                case Level.allAbsolute
                    beatIdx = idx(tab.event == Event.absBaseline1):idx(tab.event == Event.startMortality);
                case Level.absBaseline1
                    beatIdx = idx(tab.event == Event.absBaseline1):idx(tab.event == Event.absDecrease7);
                case Level.absBaseline2
                    beatIdx = idx(tab.event == Event.absIncrease7):idx(tab.event == Event.startMortality);
                case Level.absDecrease7
                    if ~isnan(idx(tab.event == Event.absDecrease14))
                        beatIdx = idx(tab.event == Event.stopBleed7):idx(tab.event == Event.absDecrease14);
                    else
                        beatIdx = idx(tab.event == Event.stopBleed7):idx(tab.event == Event.absIncrease7);
                    end
                case Level.absDecrease14
                    if ~isnan(idx(tab.event == Event.absDecrease21))
                        beatIdx = idx(tab.event == Event.stopBleed14):idx(tab.event == Event.absDecrease21);
                    else
                        beatIdx = idx(tab.event == Event.stopBleed14):idx(tab.event == Event.absIncrease14);
                    end
                case Level.absDecrease21
                    if ~isnan(idx(tab.event == Event.absDecrease28))
                        beatIdx = idx(tab.event == Event.stopBleed21):idx(tab.event == Event.absDecrease28);
                    else
                        beatIdx = idx(tab.event == Event.stopBleed21):idx(tab.event == Event.absIncrease21);
                    end
                case Level.absDecrease28
                    beatIdx = idx(tab.event == Event.stopBleed28):idx(tab.event == Event.absIncrease28);
                case Level.absIncrease28
                    beatIdx = idx(tab.event == Event.absIncrease28):idx(tab.event == Event.absIncrease21);
                case Level.absIncrease21
                    beatIdx = idx(tab.event == Event.absIncrease21):idx(tab.event == Event.absIncrease14);
                case Level.absIncrease14
                    beatIdx = idx(tab.event == Event.absIncrease14):idx(tab.event == Event.absIncrease7);
                case Level.absIncrease7
                    beatIdx = idx(tab.event == Event.absIncrease7):idx(tab.event == Event.startMortality);
                case Level.relBaseline1
                    beatIdx = idx(tab.event == Event.relBaseline1):idx(tab.event == Event.relative5);
                case Level.relBaseline2
                    beatIdx = idx(tab.event == Event.relBaseline2);
                case Level.relative5
                    if ~isnan(idx(tab.event == Event.relative10))
                        beatIdx = idx(tab.event == Event.relative5):idx(tab.event == Event.relative10);
                    else
                        beatIdx = idx(tab.event == Event.relative5):idx(tab.event == Event.relBaseline2);
                    end
                case Level.relative10
                    if ~isnan(idx(tab.event == Event.relative20))
                        beatIdx = idx(tab.event == Event.relative10):idx(tab.event == Event.relative20);
                    else
                        beatIdx = idx(tab.event == Event.relative10):idx(tab.event == Event.relBaseline2);
                    end
                case Level.relative20
                    beatIdx = idx(tab.event == Event.relative20):idx(tab.event == Event.relBaseline2);
            end
            
        else
            
            % For the specified level...
            switch level

                case Level.all
                    if strcmp(subject, "subject1")
                        beatIdx = idx(tab.event == Event.relBaseline1):idx(tab.event == Event.startMortality);
                    else; beatIdx = nan;
                    end
                case Level.allRelative
                    beatIdx = idx(tab.event == Event.relBaseline1);
                case Level.allAbsolute
                    beatIdx = idx(tab.event == Event.absBaseline1):idx(tab.event == Event.startMortality);
                case Level.absBaseline1
                    beatIdx = idx(tab.event == Event.absBaseline1):idx(tab.event == Event.absDecrease7);
                case Level.absBaseline2
                    beatIdx = idx(tab.event == Event.absIncrease7):idx(tab.event == Event.startMortality);
                case Level.absDecrease7
                    if ~isnan(idx(tab.event == Event.absDecrease14))
                        beatIdx = idx(tab.event == Event.absDecrease7):idx(tab.event == Event.absDecrease14);
                    else
                        beatIdx = idx(tab.event == Event.absDecrease7):idx(tab.event == Event.absIncrease7);
                    end
                case Level.absDecrease14
                    if ~isnan(idx(tab.event == Event.absDecrease21))
                        beatIdx = idx(tab.event == Event.absDecrease14):idx(tab.event == Event.absDecrease21);
                    else
                        beatIdx = idx(tab.event == Event.absDecrease14):idx(tab.event == Event.absIncrease14);
                    end
                case Level.absDecrease21
                    if ~isnan(idx(tab.event == Event.absDecrease28))
                        beatIdx = idx(tab.event == Event.absDecrease21):idx(tab.event == Event.absDecrease28);
                    else
                        beatIdx = idx(tab.event == Event.absDecrease21):idx(tab.event == Event.absIncrease21);
                    end
                case Level.absDecrease28
                    beatIdx = idx(tab.event == Event.absDecrease28):idx(tab.event == Event.absIncrease28);
                case Level.absIncrease28
                    beatIdx = idx(tab.event == Event.absIncrease28):idx(tab.event == Event.absIncrease21);
                case Level.absIncrease21
                    beatIdx = idx(tab.event == Event.absIncrease21):idx(tab.event == Event.absIncrease14);
                case Level.absIncrease14
                    beatIdx = idx(tab.event == Event.absIncrease14):idx(tab.event == Event.absIncrease7);
                case Level.absIncrease7
                    beatIdx = idx(tab.event == Event.absIncrease7):idx(tab.event == Event.startMortality);
                case Level.relBaseline1
                    beatIdx = idx(tab.event == Event.relBaseline1):idx(tab.event == Event.relative5);
                case Level.relBaseline2
                    beatIdx = idx(tab.event == Event.relBaseline2);
                case Level.relative5
                    if ~isnan(idx(tab.event == Event.relative10))
                        beatIdx = idx(tab.event == Event.relative5):idx(tab.event == Event.relative10);
                    else
                        beatIdx = idx(tab.event == Event.relative5):idx(tab.event == Event.relBaseline2);
                    end
                case Level.relative10
                    if ~isnan(idx(tab.event == Event.relative20))
                        beatIdx = idx(tab.event == Event.relative10):idx(tab.event == Event.relative20);
                    else
                        beatIdx = idx(tab.event == Event.relative10):idx(tab.event == Event.relBaseline2);
                    end
                case Level.relative20
                    beatIdx = idx(tab.event == Event.relative20):idx(tab.event == Event.relBaseline2);
            end
            
        end
        
    end

% Extract beat indices for all levels a subject has
    function allLevelsToBeats(subject, include)
        
        % For each level this subject has, find the indices of each level
        levels = enumeration('Level');  % Placeholder for all possible levels
        for l = 1:length(levels)

            % Extract the current level
            level = levels(l);

            % If the subject is not 1, determine the set from the parent
            if subject ~= 1 && obj.getParent(level) == Level.allRelative; set = 1; end
            if subject ~= 1 && obj.getParent(level) == Level.allAbsolute; set = 2; end
            
            % Get the beat indices for the current level and save
            beatIdx = levelToBeats(level, subject, include);
            if subject == 1 && (level == Level.allRelative || level == Level.relBaseline2); ...
                    beatIdx = beatIdx(1):11251; end
            if subject ~= 1 && (level == Level.allRelative || level == Level.relBaseline2); ...
                    beatIdx = beatIdx(1):length(beats_b{set}); end
            if isnan(beatIdx(1)); continue; end     % Continue if the level does not exist
            obj.dataset{counter}.(string(level)).indices = beatIdx;
            
            % For segmented data, return beat indices
            if obj.beatSeparated
                obj.dataset{counter}.(string(level)).indices = beatIdx;
            else
                % Get the sample indices
                if subject ~= 1; samples_b = beats_b{set}; samples_t = beats_t{set}; else; ...
                        samples_b = beats_b; samples_t = beats_t; end
                % For unsegmented data, return sample indices
                obj.dataset{counter}.(string(level)).samples_b = ...
                    [samples_b(beatIdx(1)) samples_b(beatIdx(end))];
                obj.dataset{counter}.(string(level)).samples_t = ...
                    [samples_t(beatIdx(1)) samples_t(beatIdx(end))];
            end

        end
        
    end

end