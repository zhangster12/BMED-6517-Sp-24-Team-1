function slidingMovie(signal1,Fs,start,speed,varargin)
% -------------------------------------------------------------------------
% Created by: Jacob Kimball
% Date Created: July 15, 2019
% Inan Research Lab @ GaTech
% -------------------------------------------------------------------------
% This function streams time series data with the following inputs: 
% signal1:      signal to plot, a vector of size 1xN or Nx1
% Fs:           sampling frequency for signal1
% start:        sample number between 1 and N 
% speed:        options are 'slow','medium','medfast', and 'fast'

% Optional Inputs: 
% 's2'          second signal to plot. must be the same length as signal1. 
% 'f1'          feature to track on signal1. 
% 'f2'          feature to track on signal2. 
%               Note: f1 and f2 must be vectors of size 1xM or Mx1 where
%                     M<N and f(i) is a sample number between 1 and N
% 'f1_style'    plotting style for feature 1
% 'f2_style'    plotting style for feature 2
%               Note: standard plot options such as 'ok' or '^c'apply, or
%                     'vert' which will create a vertical dotted line
% 'grid'        FLAG. asserts 'grid on' in the plot window
% 'window'      size of the plot window in number of samples. default 50000
% 'downsample'  integer value for downsampling to increase processing speed. 
%               Note: suggested options are 2 or 4. All signals, features, 
%                     and the start value are automatically updated with the 
%                     choice of downsample. 
% -------------------------------------------------------------------------
% NOTE: the number at the top of the figure is the true sample number at 
%       the left side of the window


% Example use case: I use this function frequently to compare an ECG signal 
% and its R peaks with a pressure or mechanical waveform to validate the 
% choice of R peaks. (useful with noisy ECG signals) In this case, I would
% run the command below: 

% cardio.general.sliding_movie(ECG,Fs,1,'medfast','f1',Rpeaks,'s2',SCG,'grid')
% -------------------------------------------------------------------------

% Parse optional arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 's2'); signal2 = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'f1'); feature1 = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'f2'); feature2 = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'f1_style'); f1_style = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'f2_style'); f2_style = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'grid'); grid_on = true;
        elseif strcmp(varargin{arg}, 'window'); window = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'downsample'); downsample = varargin{arg + 1};
        end
    end
end

% Set defaults
if ~exist('signal2', 'var'); signal2 = NaN(1,length(signal1)); end
if ~exist('feature1', 'var'); feature1 = 0; end
if ~exist('feature2', 'var'); feature2 = 0; end
if ~exist('f1_style', 'var'); f1_style = 'ob'; end
if ~exist('f2_style', 'var'); f2_style = 'oc'; end
if ~exist('grid_on', 'var'); grid_on = false; end
if ~exist('window', 'var'); window = Fs*5; end
if ~exist('downsample', 'var'); downsample = 1; end
pauselength = 0.0001*2000/Fs; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if(strcmp(speed,'slow')) % ~1x
    tau = 50/downsample;
else
    if(strcmp(speed,'medium')) % ~1.5x
        tau = 120/downsample;
    else
        if(strcmp(speed,'medfast')) % ~2x
            tau = 170/downsample;
        else
            if(strcmp(speed,'fast')) % ~2.5x
                tau = 300/downsample;
            else
                disp('ERROR: speed options are "slow", "medium", "medfast", and "fast"')
            end
        end
    end
end

% downsample everything
signal1 = resample(signal1,1,downsample); 
signal2 = resample(signal2,1,downsample); 
feature1 = feature1/downsample; 
feature2 = feature2/downsample; 
start = start/downsample; 
window = window/downsample; 

sample = start;

% normalize signal 2
range1 = max(signal1(1:window))-min(signal1(1:window));
range2 = max(signal2(1:window))-min(signal2(1:window));
divider = range2/range1; 
sig2 = signal2/divider; % make them have the same 'variance'
mean1 = mean(signal1);
mean2 = mean(sig2);
diffmean = mean1-mean2;
sig2 = sig2+diffmean+max(signal1)+2*mean1; % offset signal 1 and signal 2


% create feature vectors in advance
feat1 = NaN(1,length(signal1))';
if(length(feature1)>1)
    % create a vector of nans and signal1(feature1)
    for i = 1:length(feature1)
        feat1(int32(feature1(i))) = signal1(int32(feature1(i)));
    end
end
feat2 = NaN(1,length(signal2))';
if(length(feature2)>1)
    % create a vector of nans and signal1(feature1)
    for i = 1:length(feature2)
        feat2(int32(feature2(i))) = sig2(int32(feature2(i)));
    end
end

% adjust for vertical features
if(strcmp(f1_style,'vert'))
    f1_style = '.b'; 
    test = [feat1' feat1']'; 
    if(length(test) > length(signal1))
        test = [feat1]; 
        for i = 1:6
            test = [test (feat1+i*2*range1/6)];
        end
    else
        test = [feat1']; 
        for i = 1:6
            test = [test (feat1+i*2*range1/6)'];
        end
    end
    feat1 = test; 
end

if(strcmp(f2_style,'vert'))
    f2_style = '.c'; 
    test = [feat2' feat2']'; 
    if(length(test) > length(signal2))
        test = [feat2]; 
        for i = 1:6
            test = [test (feat2-i*2*range1/6)];
        end
    else
        test = [feat2']; 
        for i = 1:6
            test = [test (feat2-i*2*range1/6)'];
        end
    end
    feat2 = test; 
end
    
figure(); % watch what happens
for i = start:length(floor(signal1-window)/tau)
    plot(signal1(sample:sample+window),'k')
    if(grid_on); grid on; end
    hold on;
    plot(feat1(sample:sample+window,:),f1_style)
    plot(sig2(sample:sample+window),'k')
    plot(feat2(sample:sample+window,:),f2_style)
    title(sample*downsample)
    sample = sample+tau;
    pause(pauselength)
    clf('reset');
end
end
