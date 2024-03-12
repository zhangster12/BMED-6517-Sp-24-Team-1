function [out] = timeDomainHRV_skb(RR)
% Adopted from  https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
% and from SKB's previous codes

% RR is the processed RR intervals in seconds
% Processed means:spikes and absurd beats are removed. But not detrended


out=struct();

            out.HRmean=mean(60./RR); % hear rate in BPM. Then taking the mean HR 
            out.NNmean = mean(RR.* 1000); % compute and convert to ms
            out.NNmedian = median(RR.* 1000); % compute and convert to ms
            out.NNmode= mode(RR.* 1000); % compute and convert to ms
            out.NNvar = var(RR); % compute and convert to ms^2
            out.NNskew = skewness(RR.*1000); 
            out.NNkurt = kurtosis(RR.*1000); 
            out.NNiqr = iqr(RR.* 1000); % compute and convert to ms
            out.SDNN = std(RR.* 1000); % compute and convert to ms % SDNN should only be done on longer data segments

            % RMSSD
            out.RMSSD = runrmssd(RR.* 1000); % compute and convert to ms

            % pNN50
            [p,n] = pNNx(RR,50); % pNN50 calculation
            out.pnn50 = p; % 
            out.nn50 = n; % 

            % Poincare
            [SD1, SD2, SD_rat]=poincare(RR);

            out.SD1=SD1;
            out.SD2=SD2;
            out.SDratio=SD_rat;

%             m=3; r=0.2*std(RR);
%             m=1; r=0.2*std(RR); %%% apen sig. for VNS
            m=2; r=0.2*std(RR); %%% apen sig. for VNS
            
            
            out.sampen= fastSampen(RR, m, r); % m=3, r=0.2*std(RR)
            out.apen =  ApproxEntropy(RR, m, r);

            % DFA: alpha1 and alpha2
            out.alpha1 = dfaScalingExponent(RR, 4, 15, 0); % minBox=4, midBox=15
            out.alpha2 = dfaScalingExponent(RR, 15); % midbox to maxbox (maxbox length is the entire RR)

%             out.alpha1 = dfaScalingExponent(RR, 4, 16, 0); % minBox=4, midBox=15
%             out.alpha2 = dfaScalingExponent(RR, 16, 64,0); % midbox to maxbox (maxbox length is the entire RR)

end



function rm = runrmssd(rr)
dif = diff(rr);
rm = sqrt(mean(dif.*dif));

end



function [p,n] = pNNx(RR,x)
%pNNx: percentage of successive/adjacent NN intervals differing by x (ms) 
%or more

    df=diff(RR); %successive ibi diffs (ms)    
    n=sum(abs(df)>=x/1000);
    p=(n/length(df))*100;
end


%% ADD poincare, non-linear and then freqDomain (separately)

function [SD1, SD2, SD_rat]=poincare(RR)

SDSD = std(diff(RR));
SDRR = std(RR);
SD1 = (1 / sqrt(2)) * SDSD; % measures the width of poincare cloud
SD2 = sqrt((2 * SDRR^2) - (0.5 * SDSD^2)); % measures the length of the poincare cloud

SD_rat = SD1/SD2;

% Convert to ms
SD1 = SD1 * 1000;
SD2 = SD2 * 1000;

end





function alpha = dfaScalingExponent(x, minBoxSize, maxBoxSize, pflag)
%
% varargout = dfaScalingExponent(xminBoxSize, midBoxSize, maxBoxSize, pflag) 
% calculates the detrended fluctuation analysis estimate of the scaling 
% exponent alpha. 
%
% INPUTS
%         x          : A Nx1 vector containing the series to be analyzed
%         minBoxSize : Smallest box width (default: 4)
%         maxBoxSize : Largest box width (default: N/4)
%         pflag      : (Optional) pflag=1 plot,  pflag=0 
% OUTPUTS     
%         alpha      : estimate of scaling exponent, +
%                      minBoxSize <= n <= maxBoxSize
%
% The raw time series x(i) is first integrated to give y(i); i=1,...,N. 
% For each length scale, n, y(i) is divided into segments of equal length, n.
% In each segment, the data is detrended by subtracting the local linear least 
% squares fit, yn(k).  The root-mean-square fluctuation of this integrated 
% and detrended time series is given by 
% F(n) = sqrt( (1/N) sum_{k=1}^N [y(k) - yn(k)]^2 )
% We calculate the average fluctuation F(n) for each segment n. 
% If the scaling approximately given by F(n) = c n^alpha, 
% we can estimate alpha by calculating the slope of log F(n) versus log n.
% Such a linear relationship on a log-log plot indicates the presence of 
% power law (fractal) scaling. 
% A log-log plot of F(n) against n is provided when pflag=1.  Default: plag=0.
% Peng C-K, Buldyrev SV, Havlin S, Simons M, Stanley HE, Goldberger AL. 
% Mosaic organization of DNA nucleotides. Phys Rev E 1994;49:1685-1689.
%
%
% 09-20-2017 Modified by Giulia Da Poian (GDP) to be included in the Physionet 
%            HRV Toolkit for Matlab. (Original function name: dfa)
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
% Copyright (c) 2005 Patrick E. McSharry (patrick@mcsharry.net)
%
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information

if nargin < 2 || isempty(minBoxSize)
    minBoxSize = 4;
end
if nargin < 3 || isempty(maxBoxSize)
    maxBoxSize = length(x)/4;
end
if nargin < 4
   pflag = 0;
end

if size(x,1)<size(x,2)
    x=x';
end

N = length(x);     
y = cumsum(x);

n1 = round(log2(minBoxSize)); % modified GDP, was 3
n2 = round(log2(maxBoxSize)); % modified GDP, was n2 = round(log2(N/2))
ns = (2.^(n1:n2))';           % modified GDP, was ns =[2.^[n1:n2] N]' 

nn = length(ns);
F = zeros(nn,1);
for n=1:nn
   t = trend(y, ns(n));
   z = y - t;
   F(n) = sqrt(mean(z.^2));
 
end

lns = log10(ns);
lF = log10(F);
A = ones(nn,2);
A(:,2) = lns;
a = pinv(A)*lF;
alpha = a(2);  
lFpred = A*a;

 

if pflag == 1
    figure;
    loglog(10.^lns, 10.^lF,'b.-','MarkerSize',16);
    hold on;
    loglog(10.^[lns(1) lns(nn)], 10.^[lFpred(1) lFpred(nn)],'k');
    xlabel('n');
    ylabel('F(n)');
    title(['F(n) ~ n^{\alpha} with \alpha = ' num2str(a(2)) ]);
end

end % dfaScalingExponent function



function t = trend(y, n)
    N = length(y);
    t = zeros(N,1);
    r = floor(N/n);
    for i=1:r 
       v = y((i-1)*n+1:i*n);
       t((i-1)*n+1:i*n) = linfit(v); 
    end 
    v = y(r*n+1:N);
    t(r*n+1:N) = linfit(v);
end % trend function
   
function up = linfit(v)
    k = length(v);
    A = ones(k,2);
    u = [1:k]';
    A(:,2) = u;
    a = pinv(A)*v;
    up = A*a;
end % linfit function



