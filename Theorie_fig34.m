function Theorie_fig34()
% fig34: plots figure 3 and 5 of the paper
clc; clear all; close all;
set(0,'DefaultAxesFontSize',16);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Theorie_WSST2_2017/elsarticle/figures';


N = 1024;
t = (0:N-1)/N; t  = t(:);
fs = 0:N/2;
nv = 32;
mywav = 'cmor6-1';
%% Tests signals

[a1,a2,a3,if1,if2,if3,s1,s2,s3,ss,iff1,iff2,iff3]  = signal_test(t,2);
st1 = s1;

% Parameters
gamma = 0.001
sigma = 0.05;
%index = 2*N/8+1:6*N/8;
index = round(0.2*N):round(0.8*N);
% Choice of time and frequency bins
ft =1:N/2;bt=1:N;

%% Sparsity (no noise)
nc = 6;
es = 150;

s = st1;
%[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3] = sstn(s,gamma,sigma,ft,bt);
[STFT,FSST,FSST2,omegaF,tauF,omega2F,g] = sst2_new(s,1/sigma^2/N,N,gamma);
% FSST = FSST(1:N/2,:); 
% FSST2 = FSST2(1:N/2,:); 
[WT, WSST, WSST2, ~, ~, ~, ~, ~] = Wsst2_new(s,gamma,mywav,nv);

plotsparse_new(FSST(:,index),FSST2(:,index),WSST(:,index),WSST2(:,index),nc,es);
set(gca,'XTick',0:0.5:nc)
legend('FSST','FSST2', 'WSST','WSST2','Location','southeast');
xlabel('number of coefficients / M');
ylabel('normalized energy');
%title('chirp');
%%%%%%%%%%%%%%%%%%%%%%print SS
export_fig(gca, ... % figure handle
    sprintf('%s/SS_1', chemin0),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...         % file format
    '-r5000' );             % resolution in dpi

%%
s = s2;
%[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3] = sstn(s,gamma,sigma,ft,bt);
[STFT,FSST,FSST2,omegaF,tauF,omega2F,g] = sst2_new(s,1/sigma^2/N,N,gamma);
% FSST = FSST(1:N/2,:); 
% FSST2 = FSST2(1:N/2,:);
[WT, WSST, WSST2, ~, ~, ~, ~, ~] = Wsst2_new(s,gamma,mywav,nv);

plotsparse_new(FSST(:,index),FSST2(:,index),WSST(:,index),WSST2(:,index),nc,es);
set(gca,'XTick',0:0.5:nc)
legend('FSST','FSST2', 'WSST','WSST2','Location','southeast');
xlabel('number of coefficients / M');
ylabel('normalized energy');
%title('mode 2');

%%%%%%%%%%%%%%%%%%%%%%print SS
export_fig(gca, ... % figure handle
    sprintf('%s/SS_2', chemin0),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...         % file format
    '-r5000' );             % resolution in dpi

%%
s = s3;
%[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3] = sstn(s,gamma,sigma,ft,bt);
[STFT,FSST,FSST2,omegaF,tauF,omega2F,g] = sst2_new(s,1/sigma^2/N,N,gamma);
% FSST = FSST(1:N/2,:); 
% FSST2 = FSST2(1:N/2,:);
[WT, WSST, WSST2, ~, ~, ~, ~, ~] = Wsst2_new(s,gamma,mywav,nv);

plotsparse_new(FSST(:,index),FSST2(:,index),WSST(:,index),WSST2(:,index),nc,es);
set(gca,'XTick',0:0.5:nc)
legend('FSST','FSST2', 'WSST','WSST2','Location','southeast');
xlabel('number of coefficients / M');
ylabel('normalized energy');
%title('mode 3');

%%%%%%%%%%%%%%%%%%%%%%print SS
export_fig(gca, ... % figure handle
    sprintf('%s/SS_3', chemin0),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...         % file format
    '-r5000' );             % resolution in dpi
%% Spartsity (noise)
nc = 6;
es = 3000;

% noise
SNR = 0;
%noise =hilbert(randn(N,1));
noise = randn(N,1)+1i*randn(N,1);
s1=s1(:);
[sn] = sigmerge(s1,noise,SNR);
round(snr(s1,sn-s1))
gamma = std(real(sn-s1))

%[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3] = sstn(sn,gamma,sigma,ft,bt);
[STFT,FSST,FSST2,omegaF,tauF,omega2F,g] = sst2_new(sn,1/sigma^2/N,N,gamma);
% FSST = FSST(1:N/2,:); 
% FSST2 = FSST2(1:N/2,:);
[WT, WSST, WSST2, ~, ~, ~, ~, ~] = Wsst2_new(sn,gamma,mywav,nv);

plotsparse_new(FSST(:,index),FSST2(:,index),WSST(:,index),WSST2(:,index),nc,es);
set(gca,'XTick',0:0.5:nc)
legend('FSST','FSST2', 'WSST','WSST2','Location','southeast');
xlabel('number of coefficients / M');
ylabel('normalized energy');
%title('noisy signal');

%%%%%%%%%%%%%%%%%%%%%%print SS
export_fig(gca, ... % figure handle
    sprintf('%s/SS_noise_1', chemin0),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...         % file format
    '-r5000' );             % resolution in dpi


% noise
s2=s2(:);
[sn] = sigmerge(s2,noise,SNR);
round(snr(s2,sn-s2))
gamma = std(real(sn-s2))

%[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3] = sstn(sn,gamma,sigma,ft,bt);
[STFT,FSST,FSST2,omegaF,tauF,omega2F,g] = sst2_new(sn,1/sigma^2/N,N,gamma);
% FSST = FSST(1:N/2,:); 
% FSST2 = FSST2(1:N/2,:);
[WT, WSST, WSST2, ~, ~, ~, ~, ~] = Wsst2_new(sn,gamma,mywav,nv);

plotsparse_new(FSST(:,index),FSST2(:,index),WSST(:,index),WSST2(:,index),nc,es);
set(gca,'XTick',0:0.5:nc)
legend('FSST','FSST2', 'WSST','WSST2','Location','southeast');
xlabel('number of coefficients / M');
ylabel('normalized energy');
%title('noisy signal');

%%%%%%%%%%%%%%%%%%%%%%print SS
export_fig(gca, ... % figure handle
    sprintf('%s/SS_noise_2', chemin0),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...         % file format
    '-r5000' );             % resolution in dpi

% noise
s3=s3(:);
[sn] = sigmerge(s3,noise,SNR);
round(snr(s3,sn-s3))
gamma = std(real(sn-s3))

%[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3] = sstn(sn,gamma,sigma,ft,bt);
[STFT,FSST,FSST2,omegaF,tauF,omega2F,g] = sst2_new(sn,1/sigma^2/N,N,gamma);
% FSST = FSST(1:N/2,:); 
% FSST2 = FSST2(1:N/2,:);
[WT, WSST, WSST2, ~, ~, ~, ~, ~] = Wsst2_new(sn,gamma,mywav,nv);

plotsparse_new(FSST(:,index),FSST2(:,index),WSST(:,index),WSST2(:,index),nc,es);
set(gca,'XTick',0:0.5:nc)
legend('FSST','FSST2', 'WSST','WSST2','Location','southeast');
xlabel('number of coefficients / M');
ylabel('normalized energy');
%title('noisy signal');

%%%%%%%%%%%%%%%%%%%%%%print SS
export_fig(gca, ... % figure handle
    sprintf('%s/SS_noise_3', chemin0),... % name of output file without extension
    '-painters', ...      % renderer
    '-transparent', ...   % renderer
    '-pdf', ...         % file format
    '-r5000' );             % resolution in dpi

end
