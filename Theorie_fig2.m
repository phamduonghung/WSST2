% fig 2 : draws Figures of paper "Second-order Synchrosqueezing Transform: The
%Wavelet Case and Comparison, by PHAM and Meignen.

close all; clc; clear all;
set(0,'DefaultAxesFontSize',16);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Theorie_WSST2_2017/elsarticle/figures';

%% Parameters
N = 1024;
gamma = 0.001;
sigma = 0.05;
nv = 32;
mywav = 'cmor6-1';
nmodes = 3;

%% set up data
t  = (0:N-1)/N; t = t(:);
%fs = (0:N-1)/N;
index = round(0.2*N):round(0.8*N);
% Choice of time and frequency bins
ft =1:N/2;bt=1:N;

%% Test signal 
[a1,a2,a3,if1,if2,if3,s1,s2,s3,s,iff1,iff2,iff3] = signal_test(t,1);

%% noise
SNR = 10;
noise = randn(1,N);
sb = sigmerge(s(:),noise(:),SNR);
disp([num2str(round(20*log10(norm(s)/norm(sb-s(:))))) ' dB for input SNR']);

%% noisefree signal TFRs
%STFT-based
%[STFT,FSST,FSST2,FSST3,FSST4,omegaF,omega2F,omega3F,tau2F,tau3F,phi22pF,phi33pF,phi44pF] = sstn(s,gamma,sigma,ft,bt);
[STFT,FSST,FSST2,omegaF,tauF,omega2F] = sst2_new(s,1/sigma^2/N,N,gamma);
STFT = STFT(1:N/2,:);
FSST = FSST(1:N/2,:);
FSST2 = FSST2(1:N/2,:);

%[~, ~, FRM, ~, ~, ~, ~, ~, ~,] = sst2_Thomas(s,gamma,sigma); 

%WT-based
[WT, WSST, WSST2, fs, as, omega, tau, phipp] = Wsst2_new(s,gamma,mywav,nv);
%[~, ~, WRM, ~, ~, ~, ~, ~, ~, ~] = Wsst2_old(s,gamma,mywav,nv);
%[~, ~,WRM, ~, ~, ~, ~] = Wsst(s,gamma,mywav,nv);

%% Display TFRs - Figure 1
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.05], [0.05 0.1], [0.05 0.05]);
% if ~make_it_tight,  clear subplot;  end
% 
% FigHandle1(1) = figure; %colormap(1-gray);
% set(FigHandle1(1),'units','normalized','outerposition',[0 0 1 1]);

% subplot(2,3,1);
FigHandle(1) = figure(); 
imagesc(t,fs,abs(STFT)); %title('(a) STFT');
set(gca,'YDir','normal'); %set(gca, 'xtick', []) ;
%set(gca,'xtick',[0:0.2:1],'xlim',[0.2 .8]);

% subplot(2,3,2);
FigHandle(2) = figure();
imagesc(t,fs,abs(FSST)); %title('(b) FSST');
set(gca,'YDir','normal');%set(gca, 'xtick', []) ; 

% subplot(2,3,3);
FigHandle(3) = figure();
imagesc(t,fs,abs(FSST2)); %title('(c) FSST2');
set(gca,'YDir','normal'); %set(gca, 'xtick', []) ;

% subplot(2,3,3);
%FigHandle(7) = figure();
%imagesc(t,fs,abs(FRM)); %title('(c) FSST2');
%set(gca,'YDir','normal'); %set(gca, 'xtick', []) ;


%subplot(2,3,4);
FigHandle(4) = tfplot(WT,'log',fs,t);set(gca,'ylim',[nv*4 nv*9+1]); %title('(d) CWT');
set(gca, 'xtick', []) ;
%subplot(2,3,5);
FigHandle(5) = tfplot(WSST,'log',fs,t);set(gca,'ylim',[nv*4 nv*9+1]); %title('(e) WSST');
set(gca, 'xtick', []) ;
%subplot(2,3,6);
FigHandle(6) =tfplot(WSST2,'log',fs,t);set(gca,'ylim',[nv*4 nv*9+1]); %title('(f) WSST2');
set(gca, 'xtick', []) ;

%subplot(2,3,6);
%FigHandle(8) = tfplot(WRM,'log',fs,t);set(gca,'ylim',[nv*4 nv*9+1]); %title('(g) WRM');
%set(gca, 'xtick', []) ;

set(findall(0,'type','axes'),'xlim',[0.2 0.8],'xtick', []);
%set(findall(0,'type','axes'),'xtick',[0:0.2:0.8],'xlim',[0.2 0.8]);

for i = 1:6
    export_fig(FigHandle(i), ... % figure handle
        sprintf('%s/TFRs_%d', chemin0,i),... % name of output file without extension
        '-painters', ...      % renderer
        '-transparent', ...   % renderer
        '-pdf', ...         % file format
        '-r5000' );             % resolution in dpi
end