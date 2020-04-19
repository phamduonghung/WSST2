% fig56 : draws Figure 5 and 6
% It compares the ideal TF representation with the one given by different
% time-frequency representations, using the earth-mover distance.
% 
% This script needs some functions, including the fast EMD and the
% synchrosqueezed wavelet packet transform, that can be downloaded from 
% https://github.com/HaizhaoYang/SST_compare/blob/master/comparison/

%function Theorie_EMD()

clc; clear all; close all;
set(0,'DefaultAxesFontSize',16);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Theorie_WSST2_2017/elsarticle/figures';

%% Parameters
nv = 32;
mywav = 'cmor6-1';
%gamma = 0.001;
sigma = 0.05;
%% set up data
N  = 1024;
t  = (0:N-1)/N; t = t(:);

index = round(0.2*N):round(0.8*N);
%index = 2*N/8+1:6*N/8;
% Choice of time and frequency bins
ft =1:N/2;bt=1:N;
 
 %% Test signal
 iff=zeros(N,3);
 ss = zeros(N,3);
[a1,a2,a3,iff(:,1),iff(:,2),iff(:,3),ss(:,1),ss(:,2),ss(:,3),sss] = signal_test(t,1);
for i=1:3
    s = ss(:,i);
    % True IF
    InstFreq = iff(:,i);
    % Maps into scale space
    InstScale =log2(InstFreq)*nv;
 
    % simulation parameters
    numTest = 2; % Nb realizations (20)
    NMvec = -5:5:30;

    % store results
    resSTFT = zeros(length(NMvec),1);
    resFSST = zeros(length(NMvec),1);
    resFSST2 = zeros(length(NMvec),1);
    resFSST2_old = zeros(length(NMvec),1);
    
    resWT = zeros(length(NMvec),1);
    resWSST = zeros(length(NMvec),1);
    resWSST2 = zeros(length(NMvec),1);
    %resWSST2_a = zeros(length(NMvec),1);
    resWSST2_old= zeros(length(NMvec),1);
    %% parameters for SSWPT

    for cntNM = 1:length(NMvec) 
        % loop on noise
        for cntt = 1:numTest
            % different realizations
            % set noise
             b =randn(N,1)+1i*randn(N,1);
             %(std(b))^2
             [sb] = sigmerge(s,b,NMvec(cntNM));
             %sb=s;
             round(snr(s,sb-s))
             gamma = std(real(sb-s))
             
            %% Computes synchrosqueezing transforms 

%             [WT, WSST, WSST2, ~, as, ~, ~, ~,norm2psi] = Wsst2_new(sb-s,gamma,mywav,nv);
%             [STFT,FSST,FSST2,omegaF,tauF,omega2F,g] = sst2_new(sb-s,1/sigma^2/N,N,gamma);
%             
%             figure()
%             plot(1:size(omegaF,1),gamma(1)*norm(g)/N*ones(1,size(omegaF,1)))
%             hold on
%             plot(1:size(omegaF,1),std(transpose(real(STFT))),'--')
%                         
%             figure()
%             plot(1:length(as),gamma(1)*norm2psi/sqrt(2*N))
%             hold on
%             plot(1:length(as),std(transpose(real(WT))),'--')
%            pause     
           % [~, ~, FRM, ~, ~, ~, ~, ~, ~,] = sst2_Thomas(sb,gamma,sigma); 
            [STFT,FSST,FSST2,~,~,~] = sst2_new(sb,1/sigma^2/N,N,gamma);
            [STFT1,FSST1,FSST12,~,~,~] = sst2_omega_old(sb,1/sigma^2/N,N,gamma);
            
            STFT = STFT(1:N/2,:); 
            FSST = FSST(1:N/2,:); 
            FSST2 = FSST2(1:N/2,:); 
            
            STFT1 = STFT1(1:N/2,:); 
            FSST1 = FSST1(1:N/2,:); 
            FSST12 = FSST12(1:N/2,:); 
            
            [Cs12_F, Es12] = exridge_mult_Noise(FSST12, 1,0,10);
            [Cs2_F, Es2]   = exridge_mult_Noise(FSST2, 1,0,10);
            [Cs_F, Es]     = exridge_mult_Noise(FSST, 1,0,10);
 
            [WT, WSST, WSST2, fs, as, omega, omega2, tau, phipp, norm2psi]=...
             Wsst2_new(sb,gamma,mywav,nv);
            %[WT, WSST, WSST2, ~, as, ~, ~, ~,~] = Wsst2_new(sb,gamma,mywav,nv);
            [WT1, WSST1, WSST12, fs, as, omega, omega2_old, tau, phipp, norm2psi]...
            = Wsst2_omega_old(sb,gamma,mywav,nv);
            
            [Cs12, Es12] = exridge_mult_Noise(WSST12, 1,0,10);
            [Cs2, Es2]   = exridge_mult_Noise(WSST2, 1,0,10);
            [Cs, Es] = exridge_mult_Noise(WSST, 1,0,10);
 
            
%             plot(1:length(as),omega2_old(:,642),1:length(as),omega2(:,642),'--',1:length(as),omega(:,642),'-.')
%             pause
%             figure
%             imagesc(log(abs(WSST2)+10^-6))
%             set(gca,'ydir','normal');
%             figure
%             imagesc(log(abs(WSST12)+10^-6))
%             set(gca,'ydir','normal');
%             pause
            %[~, ~, WSST2_a, ~, ~, ~, ~, ~] = Wsst2_a(sb,gamma,mywav,nv);
                
             %% inspects results
            if 0
                figure;
                subplot(221);imagesc(log(0.1+abs(FSST(:,index))));title('FSST');set(gca,'YDir','normal');set(gca, 'xtick', []) ;
                subplot(222);imagesc(log(0.1+abs(FSST2(:,index))));title('FSST2'); set(gca,'YDir','normal');set(gca, 'xtick', []) ;
                subplot(223);imagesc(log(0.1+abs(WSST(:,index))));title('WSST'); set(gca,'YDir','normal');set(gca, 'xtick', []) ;
                subplot(224);imagesc(log(0.1+abs(WSST2(:,index))));title('WSST2'); set(gca,'YDir','normal');set(gca, 'xtick', []) ;
                pause
            end
        
            %% computes EMD
            %STFT-based
            % resWSST(cntNM) = resWSST(cntNM) + EMDMatbis(abs(WSST(:,index)),InstScale(index),length(as))/numTest;
             ind = zeros(length(ft),length(index));
             for kk=1:length(index),
              ind(Cs_F(1,index(kk)),kk) = 1;
             end
             resSTFT(cntNM) = resSTFT(cntNM) + EMDMatbis(abs(STFT(:,index)),InstFreq(index),N/2)/numTest;
             resFSST(cntNM) = resFSST(cntNM) + EMDMatbis(abs(FSST(:,index)).*ind,InstFreq(index),N/2)/numTest;
             ind = zeros(length(ft),length(index));
             for kk=1:length(index),
              ind(Cs2_F(1,index(kk)),kk) = 1;
             end
             resFSST2(cntNM) = resFSST2(cntNM) + EMDMatbis(abs(FSST2(:,index)).*ind,InstFreq(index),N/2)/numTest;
             ind = zeros(length(ft),length(index));
             for kk=1:length(index),
              ind(Cs12_F(1,index(kk)),kk) = 1;
             end
             resFSST2_old(cntNM) = resFSST2_old(cntNM) + EMDMatbis(abs(FSST12(:,index)).*ind,InstFreq(index),N/2)/numTest;
    
            % wavelet-based
            resWT(cntNM) = resWT(cntNM) + EMDMatbis(abs(WT(:,index)),InstScale(index),length(as))/numTest;
            resWSST(cntNM) = resWSST(cntNM) + EMDMatbis(abs(WSST(:,index)),InstScale(index),length(as))/numTest;
            ind = zeros(length(as),length(index));
            for kk=1:length(index),
             ind(Cs2(1,index(kk)),kk) = 1;
            end
            resWSST2(cntNM) = resWSST2(cntNM) + EMDMatbis(abs(WSST2(:,index)).*ind,InstScale(index),length(as))/numTest;
            ind = zeros(length(as),length(index));
            for kk=1:length(index),
             ind(Cs12(1,index(kk)),kk) = 1;
            end
            resWSST2_old(cntNM) = resWSST2_old(cntNM) + EMDMatbis(abs(WSST12(:,index)).*ind,InstScale(index),length(as))/numTest;
            %resWSST2_a(cntNM) = resWSST2_a(cntNM) + EMDMatbis(abs(WSST2_a(:,index)),InstScale(index),length(as))/numTest;

            %% inspects results
            if 0
                figure;
                subplot(221);imagesc(log(0.1+abs(FSST(:,index))));title('FSST');set(gca,'YDir','normal');set(gca, 'xtick', []) ;
                subplot(222);imagesc(log(0.1+abs(FSST2(:,index))));title('FSST2'); set(gca,'YDir','normal');set(gca, 'xtick', []) ;
                subplot(223);imagesc(log(0.1+abs(WSST(:,index))));title('WSST'); set(gca,'YDir','normal');set(gca, 'xtick', []) ;
                subplot(224);imagesc(log(0.1+abs(WSST2(:,index))));title('WSST2'); set(gca,'YDir','normal');set(gca, 'xtick', []) ;
                pause
            end
        end
    end

    set(findall(0,'type','axes'),'xtick',[],'ytick',[]);

    %% figure();
%     figure()
%     plot(NMvec,resSTFT,'b-.d',NMvec,resFSST,'k-.*',NMvec,resFSST2,'r-.s',NMvec,resWT,'b-d',NMvec,resWSST,'k-*',NMvec,resWSST2,'r-s');
%     legend('STFT','FSST','FSST2','CWT','WSST','WSST2','Location','northeast'); xlim([min(NMvec), max(NMvec)]);
%     xlabel('input SNR');ylabel('EMD'); 
% 
%     figure();
%     %set(FigHandle7, 'Position', [400,400, 800, 450]);
%     plot(NMvec,resFSST,'b*-',NMvec,resFSST2,'ro-',NMvec,resWSST,'gs-',NMvec,resWSST2,'kd-',NMvec,resWSST2_a,'m+-');
%     legend('FSST','FSST2','WSST','WSST2','WSST2_a','Location','northeast'); xlim([min(NMvec), max(NMvec)]);
%     xlabel('input SNR');ylabel('EMD'); 

    FigHandle(1) = figure();
    plot(NMvec,resFSST,'k*-',NMvec,resFSST2,'bo-',NMvec,resFSST2_old,'m.-',NMvec,resWSST,'gs-',NMvec,resWSST2,'rd-',NMvec,resWSST2_old,'cs-');
    legend('FSST','FSST2','FSST2_o','WSST','WSST2','WSST2_o','Location','northeast'); xlim([min(NMvec), max(NMvec)]);
    set(gca,'XTick',min(NMvec):5:max(NMvec))
    xlabel('input SNR');ylabel('EMD'); 
    
    FigHandle(2) = figure();
    plot(NMvec,resFSST,'k*-',NMvec,resFSST2,'bo-',NMvec,resWSST,'gs-',NMvec,resWSST2,'rd-');
    legend('FSST','FSST2','WSST','WSST2','Location','northeast'); xlim([min(NMvec), max(NMvec)]);
    set(gca,'XTick',min(NMvec):5:max(NMvec))
    xlabel('input SNR');ylabel('EMD'); 
    %%
    %%%%%%%%%%%%%%%%%%%%print EDM
    pause
    for j=1:2
        export_fig(FigHandle(j), ... % figure handle
            sprintf('%s/EMD_Gaussian_Amplitude_%d_%d', chemin0,i,j),... % name of output file without extension
            '-painters', ...      % renderer
            '-transparent', ...   % renderer
            '-pdf', ...         % file format
            '-r5000' );             % resolution in dpi
    end
end 
%end
