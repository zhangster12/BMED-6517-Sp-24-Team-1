function HRV = RRtoHRVDemo (RR_org, fs)

%%% RR_org is the RR interval series. 1-D data containing the  RR intervals
% in second
%%% fs is the sampling frequency of the ECG signal

RR_org=RR_org(:)'; %% making RR_org a row vector


[RR, tRR] = RRIntervalPreprocess_skb(RR_org,fs);

out = timeDomainHRV_skb(RR);

freqD = freqDomainHRV_skb(RR,tRR,4); %% 4 Hz is the resampling frequency


HRV(1)=out.HRmean;  %% Time domain
HRV(2)=out.NNmean;
HRV(3)=out.NNmedian;
HRV(4)=out.NNmode;
HRV(5)=out.NNvar;
HRV(6)=out.NNskew;
HRV(7)=out.NNkurt;
HRV(8)=out.NNiqr;
HRV(9)=out.SDNN;
HRV(10)=out.RMSSD;
HRV(11)=out.pnn50;
HRV(12)=out.nn50;    %%%%%% Time domain
HRV(13)=out.SD1;    %% Poincare
HRV(14)=out.SD2;
HRV(15)=out.SDratio;
HRV(16)=out.sampen;  %% Non-linear
HRV(17)=out.apen;
HRV(18)=out.alpha1;  %% DFA
HRV(19)=out.alpha2;
HRV(20)=freqD.VLF; %%%% FREQ domain
HRV(21)=freqD.LF;
HRV(22)=freqD.HF;
HRV(23)=freqD.LFHF;
HRV(24)=freqD.LFn;
HRV(25)=freqD.HFn;
HRV(26)=freqD.POW;

HRV=HRV(:);

end
