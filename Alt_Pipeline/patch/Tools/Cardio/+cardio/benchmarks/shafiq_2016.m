function [ao, ac] = shafiq_2016(scg, nst, ao_tol, ac_tol, ac_temp_tol, fs, varargin)
% ------------------------------------------------------------------------
% Function description:
%
% This function uses a sliding template methodology to find AO and AC from
% heart-beat segmented SCG signals. At each beat, a template is formed
% using the preceeding beats. AO/AC in the current template are then 
% identified in the template using the AO/AC from the previous template.
% AO/AC in the current beat is then identified using AO/AC from the current
% template. In the original algorithm the AC template is made through an
% ensemble average of AO align
%
% Notes: 
%
% Optional arguments: easy fixes to make the algorithm a bit more 
% robust to actually be used. Original algorithm has a few limitations: 
% 1. With noise, using just previous beat for AO detection is unreliable. 
% 2. Using R-peak aligned SCG beats for AC detection usually works better
% than using AO-peak aligned SCG beats (especially if AO detection is
% wrong)
%
% Input arguments:
% - scg             NxM             M signals of length N 
% - nst             dbl             number of templates to include  
% - ao_tol          (ms)            AO template-to-template/beat tolerance 
% - ac_tol          (ms)            AC template-to-template tolerance 
% - ac_temp_tol     (ms)            AC template-to-beat tolerance 
% - fs              (Hz)            sampling frequency
% Output arguments:
% - ao              Mx1 (samp)      AO points       
% - ac              Mx1 (samp)      AC points
% Optional arguments: 
% - align           flg             determines whether to align by r peaks 
% - independent     flg             determine mo 
% - ao_slide        flg             determine if ao sliding template made
%
% Reference: 
%
% Shafiq, G., Tatinati, S., Ang, W. et al. Automatic Identification of 
% Systolic Time Intervals in Seismocardiogram. Sci Rep 6, 37524 (2016).
% 
% Usage: 
%
% nst = 30; ao_tol = 10, ac_temp_tol = 24, ac_tol = 20,
% [ao,ac] = shafiq(scg, nst, ao_tol, ac_temp_tol, ac_tol);
% [ao,ac] = shafiq(scg, nst, ao_tol, ac_temp_tol, ac_tol, 'align', 'slide')
% ------------------------------------------------------------------------
    
    % parse optional arguments
    if ~isempty(varargin)
        for arg = 1:length(varargin)
            if strcmp('align', varargin{arg}); ao_align = false; end
            if strcmp('slide', varargin{arg}); ao_slide = true; end
            if strcmp('plot', varargin{arg}); display = true; end
        end
    end
    
    % default arguments
    if ~exist('ao_align', 'var'); ao_align = true; end
    if ~exist('ao_slide', 'var'); ao_slide = false; end
    if ~exist('display', 'var'); display = false; end
    
    %%% --------------------------------------------------------
    %%% Initialization 
    %%% --------------------------------------------------------
    
    % turn tolerances from ms to samples
    ao_tol = round(ao_tol*fs/1000);
    ac_tol = round(ac_tol*fs/1000);
    ac_temp_tol = round(ac_temp_tol*fs/1000);
    
    % save dimensions of scg 
    num_beats = size(scg, 2);
    num_points = size(scg, 1);
    
    % beat buffer for template creation
    ao_buf = zeros(size(scg, 1), nst)*nan;
    ac_buf = zeros(size(scg, 1), nst)*nan;
    
    % ao/ac features 
    ao = zeros(size(scg, 1), 1);
    ac = zeros(size(scg, 1), 1);
    
    % ao/ac features in templates
    ao_temp = ao;
    ac_temp = ac;
    
    % templates
    templates_ao = zeros(size(scg));
    templates_ac = zeros(size(scg));
    
    % initialize beat buffer 
    ao_buf(:, 1:nst) = scg(:, 1:nst);
  
    % create initial template and find AO/AC
    init_temp = mean(ao_buf, 2);
    ac_start = 250*fs/1000;
    ao_init = ...
        find(init_temp(1:ac_start) == max(init_temp(1:ac_start)), 1, 'first');
    ac_init = ac_start + ...
        find(init_temp(ac_start:end) == max(init_temp(ac_start:end)), 1, 'first');   
    
    if display
        plot(init_temp); hold on;
        scatter([ao_init], init_temp([ao_init]), 'filled');
        scatter([ac_init], init_temp([ac_init]), 'filled');
        title('Initial Template with AO/AC Indicated')
        legend('Template', 'AO', 'AC')
    end
    
    % normalize from AO peak
    if ao_align
        ac_init = ac_init - ao_init + 1;  
    end
    
    % initialize starting AO and AC points 
    prev_ao_temp = ao_init; 
    prev_ac_temp = ac_init;
    
    % get the initial AO points to set up the AO aligned SCG beat buffer
    for beat = 1:nst
        
        % first find all ao in the beats used in 
        ao(beat) = find_max(scg(:, beat), prev_ao_temp, ao_tol, ao_tol); 
        
        % update the ac buffer
        ac_buf(1:num_points-ao(beat)+1, beat) = ao_buf(ao(beat):end, beat);

        
    end

    % Decide whether to use R peak aligned signal buffer 
    if ~ao_align
        ac_buf = ao_buf;
    end
    
    %%% -------------------------------------------------------------
    %%% AO/AC detection
    %%% -------------------------------------------------------------
    
    for beat = 1:num_beats
        
        % Make the AO/AC template (reference from R/AO peak) 
        ao_slide_temp = mean(ao_buf, 2, 'omitnan'); 
        ac_slide_temp = mean(ac_buf, 2, 'omitnan');
        
        
        if ao_slide
            
            % find AO in the current template using AO from the previous template
            curr_ao_temp = find_max(ao_slide_temp, prev_ao_temp, ao_tol, ao_tol);
            
        else
            
            % otherwise find AO in the current template using AO from the
            % previous beat
            if beat == 1
                curr_ao_temp = prev_ao_temp;
            else
                curr_ao_temp = ao(beat-1);
            end
            
        end
        
        % find AO in SCG beat using AO from the past beat / current template
        ao(beat) = find_max(scg(:, beat), curr_ao_temp, ao_tol, ao_tol);
        
        % save R peak aligned template and corresponding AO point
        templates_ao(:, beat) = ao_slide_temp;
        ao_temp(beat) = curr_ao_temp;
        
        if ao_align
            
            % Create the AO aligned SCG beat using found AO
            ac_scg = zeros(size(scg, 1), 1)*nan;
            ac_scg(1:num_points-ao(beat)+1) = scg(ao(beat):end, beat);
            
        else
            
            % otherwise just used R aligned SCG beat 
            ac_scg = scg(:, beat);
            
        end
       
        % find AC in the current template using AC from the previous template
        curr_ac_temp = find_max(ac_slide_temp, prev_ac_temp, ac_temp_tol, ac_temp_tol);
        
        % save AO aligned beat and corresponding AC 
        templates_ac(:, beat) = ac_slide_temp;
        ac_temp(beat) = curr_ac_temp;
        
        % find AC point on actual AO aligned beat 
        ac(beat) = find_max(ac_scg, curr_ac_temp, ac_tol, ac_tol);
      
        
        if ao_align 
            
            % correct AO to be from R peak
            ac(beat) = ac(beat) + ao(beat)-1;
            
        end
        
        % update the new reference AO/AC point 
        prev_ao_temp = curr_ao_temp;
        prev_ac_temp = curr_ac_temp;
        
        % update the signal buffers 
        ao_buf(:, 1:nst) = [ao_buf(:, 2:nst), scg(:, beat)];
        ac_buf(:, 1:nst) = [ac_buf(:, 2:nst), ac_scg];
        
    end
end


function new_point = find_max(target_signal, ref_point, lw_tol, hi_tol)
% Function Description:
% 
% Find the location of the largest peak in the target signal that is within
% the bounds [ref_point - lw_tol, ref_point + hi_tol]. If no valid
% locations are found, then the new location is set as the old location

    % indicate the target signal bounds 
    bnds = [max(ref_point - lw_tol, 0), min(length(target_signal), ref_point+hi_tol)];
    
    % find peaks in contstrained target signal
    [peaks, locs] = findpeaks(target_signal(bnds(1):bnds(2)));
   
    % if it's empty then set to ref value 
    if isempty(locs)
        
        % if no valid locations, set to the reference point
        new_point = ref_point;
        
    else
        
        % find the largest peak location
        new_point = locs(find(peaks == max(peaks), 1, 'first'))+bnds(1)-1;
        
    end
    
end