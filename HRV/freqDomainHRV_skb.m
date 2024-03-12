function [freqD] = freqDomainHRV_skb(RR,tRR,fs)

%% This fs is the resampling freq. which is usually 4 Hz
freqD=struct();

RR=RR(:);
tRR=tRR(:);


t2 = tRR(1):1/fs:tRR(end); %time values for interp.
y=interp1(tRR,RR,t2','spline')'; %cubic spline interpolation

% y = detrend_smoothness_prior(y.'); %% smoothness prior detrending
y=detrend(y,0); %% mean removal

[PSD1,F1] = pwelch(y,[],[],2^nextpow2(length(y)),fs);

% [PSD1,F1] = pwelch(y,256,128,256,fs); 

plot_on=0; % for plotting purpose
FT = CalcLfHfParams_trapz(PSD1,F1, plot_on);

freqD.VLF=FT(1);
freqD.LF=FT(2);
freqD.HF=FT(3);
freqD.LFHF=FT(4);
freqD.LFn=FT(5);
freqD.HFn=FT(6);
freqD.POW=FT(7);


end



function [z_stat] = detrend_smoothness_prior(z)
% codes from "An Advanced Detrending Method With Application to
% HRV Analysis"

T=length(z);
lambda=10; %% previously=10
I=speye(T);
D2=spdiags(ones(T-2,1)*[1 -2 1],[0:2], T-2, T);
z_stat=(I-inv(I+lambda^2*D2'*D2))*z;

end



function FT = CalcLfHfParams_trapz(PSD,F, plot_on)
% [ulf, vlf, lf, hf, lfhf, ttlpwr] = CalcLfHfParams(PSD, F, limits,plot_on)
%
%   OVERVIEW: Compute the frequency domain features for a given PSD and
%             frequency bans limits
%         
%   INPUT:      
%        PSD     - power spectral density 
%        F       - frequency vector
%        limits  - frequency domain analysis limits
%        plot_on - 
%
%   OUTPUT:     
%	- ulf     : (ms^2) Power in the ultra low frequency range (default < 0.003 Hz)
%	- vlf     : (ms^2) Power in very low frequency range (default 0.003 <= vlf < 0.04 Hz)
%	- lf      : (ms^2) Power in low frequency range (default 0.04Hz  <= lf < 0.15 Hz)
%	- hf      : (ms^2) Power in high frequency range (default 0.15 <= hf < 0.4 Hz)
%	- lfhf    : Ratio LF [ms^2]/HF [ms^2]
%	- ttlpwr  : (ms^2) Total spectral power (approximately <0.4 Hz)
%     
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Main script written by Adriana N. Vest
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%
% if nargin <3
    ULF = [0 .0033]; %% don't care about ULF. Taskforce paper
    VLF = [0 .04];
    
    LF = [.04 .15]; %% HRV range
    HF = [0.15 0.4];
    
%     LF = [.04 .2]; %% self-custom PSD range
%     HF = [0.2 0.4];
    
    limits = [ULF; VLF; LF; HF];
% end
if nargin < 3
    plot_on =0;
end

Indx_ULF = find( (limits(1,1) <= F) & (F <= limits(1,2)) );
Indx_VLF = find( (limits(2,1) <= F) & (F <= limits(2,2)) );
Indx_LF = find( (limits(3,1) <= F) & (F <= limits(3,2)) );
Indx_HF = find( (limits(4,1) <= F) & (F <= limits(4,2)) );
% space = F(2)-F(1);

% ulf = sum(PSD(Indx_ULF)*space) * 1e6; % convert to ms^2
% vlf = sum(PSD(Indx_VLF)*space) * 1e6; % convert to ms^2
% lf = sum(PSD(Indx_LF)*space) * 1e6;   % convert to ms^2
% hf = sum(PSD(Indx_HF)*space) * 1e6;   % convert to ms^2


if length(Indx_ULF)==1
    ulf=0;
else
ulf = trapz(F(Indx_ULF),PSD(Indx_ULF))* 1e6; % Bashar edited, trapz.
end

vlf = trapz(F(Indx_VLF),PSD(Indx_VLF))* 1e6; % Bashar edited, trapz.
lf = trapz(F(Indx_LF),PSD(Indx_LF))* 1e6; % Bashar edited, trapz.
hf = trapz(F(Indx_HF),PSD(Indx_HF))* 1e6; % Bashar edited, trapz.


% ttlpwr = sum([ulf vlf lf hf]);


%% used this in OUD
ttlpwr = sum([vlf lf hf]); %% using the VLF from 0 (Ref: taskforce)


% ttlpwr = sum([lf hf])-vlf; %% subtracting the VLF (Ref: taskforce)

ttlpwr_norm=ttlpwr-vlf;
lf_n = lf/ttlpwr_norm; % normalized
hf_n = hf/ttlpwr_norm;

% lf_n = lf/ttlpwr; % normalized
% hf_n = hf/ttlpwr;


lfhf = round(lf_n/hf_n*100)/100; % lf/hf ratio


FT=[vlf, lf, hf, lfhf, lf_n, hf_n, ttlpwr];





if plot_on
    figure
    % plot PSD
%     plot(F,10*log10(PSD),'b','linewidth',2)
    plot(F,PSD,'linewidth',5.0,'Color', "#B266FF")
    
    hold on
    % plot limits on graph for lf and hf
%     plot([F(Indx_LF(1)) F(Indx_LF(1))],[-80 40],'k:')
%     hold on
%     plot([F(Indx_LF(end)) F(Indx_LF(end))],[-80 40],'k:')
%     hold on
%     plot([F(Indx_HF(end)) F(Indx_HF(end))],[-80 40],'k:')

    % labelsc
    text(0.07,30,'LF','Fontname','Times New Roman','Fontsize',10)
    text(0.25,30,'HF','Fontname','Times New Roman','Fontsize',10)
    %text(0.15, 35, 'Power Spectral Density','Fontname','Times New Roman','Fontsize',10)
    text(0.3, -60, strcat('LF/HF=',num2str(lfhf)),'Fontname','Times New Roman','Fontsize',10)
    ylabel('Normalized PSD (db/Hz)','Fontname','Times New Roman','fontsize',10)
    xlabel('Frequency (Hz)','Fontname','Times New Roman','fontsize',10)
%     axis([0 .45 -80 40]);
    xlim([0 1])
    box off
end % end plot



end % end function



