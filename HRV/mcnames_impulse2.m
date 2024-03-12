function [s_ret] = mcnames_impulse2(x,tau,N)

% code trying to implement the J. McNames impulse rejection paper, EMBC
% 2004
% Impulse Rejection Filter for Artifact Removal in Spectral Analysis of Biomedical Signals

% tau is the threshold
% window length is fixed 5 samples
x=x(:).';

x=[x(1:N),x,x(end-N+1:end)]; %% adding 1-2 and last+1+1 just for coding simplicity to avoid size issue


xm=median(x(N+1:end-N));

if median(abs(x-xm))==0
D=abs(x-xm)./(1.483*mean(abs(x-xm)));
else
D=abs(x-xm)./(1.483*median(abs(x-xm)));
end
% s_hat=x; 


for n=N+1:length(x)-N
if D(n)<tau
    s_hat(n)=x(n);

else
% n
   s_hat(n)= median(x(n-N:n+N)); 
%    s_hat(n)= mean(x(n-N:n+N)); 


%    s_hat(n)= NaN; %median(x(n-N:n+N)); 
% s_hat(n)=mean([x(n),x(n-1),x(n-2)]);


%    s_hat(n)= median(s_hat(n-N:n+N)); %% for replacing
   
    
end
    
    
end

s_ret=s_hat(N+1:end);
% s_ret=s_hat(N+1:end-N);

end

