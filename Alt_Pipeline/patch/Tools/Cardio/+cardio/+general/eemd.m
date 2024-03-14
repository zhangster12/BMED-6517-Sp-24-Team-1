function [imf,res] = eemd(sig, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'MaxNumIMF', default 10
% 'Beta', default 0.2
% 'NumEnsemble, default 100
% Implementation: David Lin
% Function 
%   Implements Ensemble Empirical Mode Decomposition (EEMD) for signal
%   decomposition. The algorithm applies various white noise series onto
%   the inputted data and applies emd on each noisy series to create 
%   different scaled signals called intrinsic mode functions (IMFs). 
%   The final IMFs from each series are then averaged together to form the
%   final IMFs.
% Required Arguments 
%   sig:        [Nx1]       signal to decompose
%
% Option Arguments 
%   MaxNumIMF   [dbl]       number of imfs to extract 
%   Beta        [dbl]       scaling constant for white noise to add
%   NumEnsemble [dbl]       number of ensembles of emd to run 
%   Custom      [dbl]       if you want to specify your own std of
%                                   white noise to add 
%   Verbose                 Output progression of algorithm
%   ShiftTol    [dbl]       shift tolerance to 
%   MaxSift     [dbl]       Max number of shifts 
%   CEEMD                   Apply complete ensembe empirical mode
%                               decomposition
%   Plot                    Plot results 
%
% Outputs
%   imf         [NxM]       Outputs M imfs of length N
%   res         [Nx1]       Ouputs residual of length N
%
% case 'MaxNumIMF'; numIMF = varargin{arg + 1};
% case 'Beta'; Beta = varargin{arg+1};
% case 'NumEnsemble'; ensemble = varargin{arg+1};
% case 'Custom'; custom = 1; customstd = varargin{arg+1};
% case 'Verbose'; display_text = 1;
% case 'SifftTol'; sifttol = varargin{arg+1};
% case 'MaxSift'; siftnum = varargin{arg+1};
% CEEMD
% Usage:
%       [imf, residual] = eemd(signal)
%       [imf, residual] = eemd(signal, Name, Value)
%                       = eemd(signal, MaxNumIMF, 7, 'NumEnsemble', 7)
% Reference
% https://pdfs.semanticscholar.org/a97e/e1d4a15c04160c323bd650e9cb9dff9dfced.pdf?_ga=2.201857681.1270078362.1582756436-1721872012.1579012722uniformly distributed white noise Z.Wu 92
% Paper: http://geogin.narod.ru/hht/link03/eemd2005.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % if sig is a rwo vec change to column vec
    if size(sig, 1) < size(sig, 2)
        sig = sig';
    end
    
    % arguments 
    numIMF = 10;
    Beta = 0.2;
    ensemble = 100;
    custom = 0;
    display_text = 0;
    sifttol = 0.2;
    ceemd = 0;
    siftnum = 100;
    Plot = 0;
    if ~isempty(varargin)
        for arg = 1:length(varargin)
            switch varargin{arg}
                case 'MaxNumIMF'; numIMF = varargin{arg + 1};
                case 'Beta'; Beta = varargin{arg+1};
                case 'NumEnsemble'; ensemble = varargin{arg+1};
                case 'Custom'; custom = 1; customstd = varargin{arg+1};
                case 'Verbose'; display_text = 1;
                case 'SifftTol'; sifttol = varargin{arg+1};
                case 'MaxSift'; siftnum = varargin{arg+1};
                case 'CEEMD'; ceemd = 1; 
                case 'Plot'; Plot = 1;
            end
        end
    end
    
    % grab a scaled std of the signal for setting the power of white noise
    std_sig = Beta * std(sig);
    if custom
        std_sig = customstd;
    end
    imf = zeros(length(sig), numIMF);
    res = zeros(size(sig));
    if ceemd 
        imf_neg = zeros(length(sig), numIMF);
        res_neg = zeros(length(sig), numIMF);
    end
    for i = 1:ensemble
        if mod(i, 10) == 0 && display_text
            fprintf('Adding Ensemble %d\n', i);
        end
        
        % set the white noise to add to the signal 
        gaus_noise = std_sig*randn(size(sig));
        
        % apply vanilla emd 
        [imfe, rese] = emd(sig+gaus_noise,...
            'MaxNumIMF',numIMF, 'SiftRelativeTol', sifttol, 'SiftMaxIterations', siftnum, 'Display', 0);
        
        % if emd decides that it has converged and the extracted imfs are fewer than anticipated
        if size(imf,2) > size(imfe, 2)
            imfe = [imfe,zeros(size(imfe, 1), size(imf, 2)-size(imfe,2))];
            rese = [rese,zeros(size(rese, 1), size(res, 2)-size(rese,2))];
        end
        
        % keep a running sum of all imfs
        imf = imf + imfe;
        res = res + rese;
        
        % apply emd on with negative noise for 0 
        if ceemd
            
            [imfe_neg, rese_neg] = emd(sig-gaus_noise,...
                'MaxNumIMF',numIMF, 'SiftRelativeTol', sifttol, 'SiftMaxIterations', siftnum, 'Display', 0);
            
            % if emd decides that it has converged and the extracted imfs are fewer than anticipated
            if size(imf_neg,2) > size(imfe_neg, 2)
                imfe_neg = [imfe_neg,zeros(size(imfe_neg, 1), size(imf_neg, 2)-size(imfe_neg,2))];
                rese_neg = [rese_neg,zeros(size(rese_neg, 1), size(res_neg, 2)-size(rese_neg,2))];
            end
            
            % keep a running sum of all subtracted noisy imfs
            imf_neg = imf_neg + imfe_neg;
            res_neg = res_neg + rese_neg;
            
        end
        
    end


    % calculate the average imfs 
    if ceemd
        imf = (imf + imf_neg)/(2*ensemble);
        res = (res + res_neg)/(2*ensemble);
    else
        imf = imf./ensemble;
        res = res./ensemble;
    end
    
    % plot the signal, imfs, and residual 
    if Plot
        figure;
        tiledlayout(ceil((size(imf, 2)+2)/2), 2, 'Padding', 'none', 'TileSpacing', 'compact'); 
        ax(1) = nexttile;
        plot(sig);
        title('Original')
        for i=1:size(imf, 2)   
            ax(i+1) = nexttile;
            plot(imf(:, i)); 
            title(['IMF ', num2str(i)])
        end
        ax(i + 1) = nexttile;
        plot(res);
        title('Residual')
        linkaxes(ax, 'x');
    end
                    
        
end