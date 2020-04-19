function Theorie_fig1()
% fig1 : draws Figure 1 of paper "Second-order Synchrosqueezing Transform: The
%Wavelet Case and Comparisons", by PHAM and Meignen.
%

clc; clear all; close all;
set(0,'DefaultAxesFontSize',14);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Theorie_WSST2_2017/elsarticle/figures';

% Parameters
N = 1024;

t = (0:N-1)/N; t  = t(:);

%% Test signal 
[a1,a2,a3,if1,if2,if3,s1,s2,s3,s,iff1,iff2,iff3] = signal_test(t,1);
 
%% display waveform signal - Figure 1 
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.05], [0.1 0.1], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end

FigHandle(1) = figure; 
set(FigHandle(1),'units','normalized','outerposition',[0 0 1 1]);
ha = subplot(4,3,1) ;
plot(t, real(s1), 'm'); 
hold on; plot(t, a1, 'k', 'linewidth', 2) ;
set(ha,'ylim',[-3 3]);
legend('Re(f_1)','A_{1}','Orientation','horizontal') ; 

FigHandle(2) = figure; 
set(FigHandle(2),'units','normalized','outerposition',[0 0 1 1]);
ha = subplot(4,3,1) ;
plot(t, real(s2), 'g'); 
hold on; plot(t, a2, 'k', 'linewidth', 2) ;
set(ha,'ylim',[-3 3]);
legend('Re(f_2)','A_{2}','Orientation','horizontal') ; 

FigHandle(3) = figure; 
set(FigHandle(3),'units','normalized','outerposition',[0 0 1 1]);
ha = subplot(4,3,1) ;
plot(t, real(s3), 'c'); 
hold on; plot(t, a3, 'k', 'linewidth', 2) ;
set(ha,'ylim',[-3 3]);
legend('Re(f_3)','A_{3}','Orientation','horizontal') ;

FigHandle(4) = figure; 
set(FigHandle(4),'units','normalized','outerposition',[0 0 1 1]);
ha = subplot(4,3,1) ;
plot(t, real(s), 'b'); 
%hold on; plot(t, a3, 'k', 'linewidth', 2) ;
set(ha,'ylim',[-3 3]);
legend('Re(f)','Orientation','horizontal') ;

% FigHandle(4) = figure; 
% set(FigHandle(4),'units','normalized','outerposition',[0 0 1 1]);
% ha = subplot(4,3,1) ;
% plot(t, real(s), 'b') ;
% legend('Re(f)') ; 
% set(ha,'ylim',[-3 3]);

%%%%%%%%%%%%%%%%%%%%%% print Figure 1

for i=1:4
    export_fig(FigHandle(i), ... % figure handle
        sprintf('%s/MCS_%d', chemin0,i),... % name of output file without extension
        '-painters', ...      % renderer
        '-transparent', ...   % renderer
        '-pdf', ...         % file format
        '-r5000' );             % resolution in dpi
end
end
