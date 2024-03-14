function sample_number = timeToSample(~, start_time, Fs, time)

% -------------------------------------------------------------------------
% This function takes a time vector and returns the appropriate sample
% number vector based on absolute start time and sampling frequency. 
%
% ARGUMENTS
% start_time    Military time start time (ex: 0805, 1315, etc) 
% Fs            Sampling freqency of desired series to match with time
% time          Vector of military time values of events
% -------------------------------------------------------------------------

% Create vector of all time-points
all_time = [start_time time]; 

% Perform transformations
extra_minutes = mod(all_time,100);              % Number of minutes past the hour
min_in_hours = 60*(all_time-extra_minutes)/100; % Convert hours to minutes
total_minutes = min_in_hours + extra_minutes;   % Total number of minutes since 0:00
total_seconds = total_minutes*60;               % Total seconds since 0:00
all_samples = total_seconds*Fs;                 % Sample number from seconds

% Return result
sample_number = all_samples(2:end)-all_samples(1); % Subtract the samples from the start time

end