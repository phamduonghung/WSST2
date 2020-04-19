%function accuracy_3modes()

%clc; 
clear all; close all;
set(0,'DefaultAxesFontSize',14);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Theorie_WSST2_2017/elsarticle/figures';

% Parameters
N = 1024;
gamma = 0.001;
sigma = 0.05;
%index = 2*N/8+1:6*N/8;
index = round(0.2*N):round(0.8*N);

d = 0;
clwin = 10;
lambda = 0; beta = 0;
nv = 32;
mywav = 'cmor6-1';
nmodes = 3;

t = (0:N-1)/N; t=t(:);
frecv = 1:N/2;

% Choice of time and frequency bins
ft =1:N/2;bt=1:N;

%% Test signal 
[a1,a2,a3,if1,if2,if3,s1,s2,s3,s,iff1,iff2,iff3] = signal_test(t,2);
Omega = [if3,if2,if1]';
Phipp = [iff3,iff2,iff1]';
%% TF presentations of signal
%[STFT,FSST,FSST2,omegaF,tauF,omega2F] = sst2_omega_old(s,1/sigma^2/N,N,gamma);
[STFT,FSST,FSST2,omegaF,tauF,omega2F] = sst2_new(s,1/sigma^2/N,N,gamma);

STFT = STFT(1:N/2,:);
FSST = FSST(1:N/2,:);
FSST2 = FSST2(1:N/2,:);

%[WT, WSST, WSST2, fs, as, omega, tau, phipp, ~]  = Wsst2_omega_old(s,gamma,mywav,nv);
[WT, WSST, WSST2, fs, as, omega, tau, phipp, ~]  = Wsst2_new(s,gamma,mywav,nv);
omega2 = omega - phipp.*tau;

 %% inspects results
if 0
    figure;
    subplot(221);imagesc(log(0.1+abs(FSST(:,index))));title('FSST');set(gca,'YDir','normal');set(gca, 'xtick', []) ;
    subplot(222);imagesc(log(0.1+abs(FSST2(:,index))));title('FSST2'); set(gca,'YDir','normal');set(gca, 'xtick', []) ;
    subplot(223);imagesc(log(0.1+abs(WSST(:,index))));title('WSST'); set(gca,'YDir','normal');set(gca, 'xtick', []) ;
    subplot(224);imagesc(log(0.1+abs(WSST2(:,index))));title('WSST2'); set(gca,'YDir','normal');set(gca, 'xtick', []) ;
    pause
end
%% Ridge extraction FSST

[Cs1, Es1] = exridge_mult_Noise(FSST, nmodes, lambda, clwin);
%[Cs1, Es1] = brevridge_mult(FSST, frecv, nmodes, lambda, clwin);
imf = recmodes(FSST,Cs1,d); imf=transpose(imf);

% disp('SNRout by FSST')
% disp(snr(s1(index),imf(index,3)- s1(index)));
% disp(snr(s2(index),imf(index,2)- s2(index)));
% disp(snr(s3(index),imf(index,1)- s3(index)));
% disp(snr(s(index),sum(imf(index,:),2)- s(index)));

figure(); 
subplot(3,1,1); 
plot(t(index),real(s1(index)),'--',t(index),real(imf(index,3)));
subplot(3,1,2)
plot(t(index),real(s2(index)),'--',t(index),real(imf(index,2)));
subplot(3,1,3)
plot(t(index),real(s3(index)),'--',t(index),real(imf(index,1)));


%% Ridge extraction FSST2

[Cs2, Es2] = exridge_mult_Noise(FSST2, nmodes, lambda, clwin);
%[Cs2, Es2] = brevridge_mult(FSST2, frecv, nmodes, lambda, clwin);
if 0 % inspects results
    for p=1:nmodes
        for q = 1:N
            IFest1(p,q) = round(real((omega2F(Cs2(p,q),q))));
        end
    end
    figure(); 
    subplot(4,1,1)
    imagesc(t,ft,abs(FSST2)); set(gca,'YDir','normal'); title('Ridge extraction FSST2') 
    hold on; plot(t,Cs2,'Linewidth',2); 
    subplot(4,1,2)
    [X,Y] = meshgrid(ft,t);
    plot(t,X); 
    plot(t,if1,'-*',t,if2,'-*',t,if3,'-*','Linewidth',1); hold on; plot(t,IFest1,'-+'); title('Case STFT: True IF and estimated IF') ; 
    subplot(4,1,3)
    title('Their diffrence') ; 
    plot(t,if1'-IFest1(3,:),'--',t,if2'-IFest1(2,:),'-.',t,if3'-IFest1(1,:),'-','Linewidth',1); 
    legend('mode 1','mode 2','mode 3'); title('Case STFT: Difference betwen true IF and estimated IF')
    
    subplot(4,1,4)
    plot(t,iff1,'-.',t,iff2,'--',t,iff3); title('Case STFT: True Frequency Modulation'); 
    %pause
end

imf2 = recmodes(FSST2,Cs2,d); imf2=transpose(imf2);
% disp('SNRout by FSST2')
% disp(snr(s1(index),imf2(index,3)- s1(index)));
% disp(snr(s2(index),imf2(index,2)- s2(index)));
% disp(snr(s3(index),imf2(index,1)- s3(index)));
% disp(snr(s(index),sum(imf2(index,:),2)- s(index)));

figure(); 
subplot(3,1,1); 
plot(t(index),real(s1(index)),'--',t(index),real(imf2(index,3)));
subplot(3,1,2)
plot(t(index),real(s2(index)),'--',t(index),real(imf2(index,2)));
subplot(3,1,3); title('Reconstruction FSST2') 
plot(t(index),real(s3(index)),'--',t(index),real(imf2(index,1)));

 
%% Ridge extraction WSST

[C1, Es3] = exridge_mult_Noise(WSST, nmodes, lambda, clwin);
[C2, Es4] = exridge_mult_Noise(WSST2, nmodes, lambda, clwin);
%[C1, Es3] = brevridge_mult(WSST, frecv, nmodes, lambda, clwin);
%[C2, Es4] = brevridge_mult(WSST2, frecv, nmodes, lambda, clwin);


if 0 % inspects results
    for p=1:nmodes
        for q = 1:N
            IFest2(p,q) = round(real((omega2(C2(p,q),q))));
        end
    end
    
    figure();
    subplot(2,1,1)
    [X,Y] = meshgrid(fs,t);
    plot(t,X); hold on;
    plot(t,if1,'-*',t,if2,'-*',t,if3,'-*','Linewidth',1); hold on; plot(t,IFest2,'--+','Linewidth',1.5); title('True IF and estimated IF in fs grid') ; 
    subplot(2,1,2)
    title('Their diffrence') ; 
    plot(t,if1'-IFest2(3,:),'--',t,if2'-IFest2(2,:),'-.',t,if3'-IFest2(1,:),'-','Linewidth',1); 
    legend('mode 1','mode 2','mode 3'); title('Difference betwen true IF and estimated IF')
  
    %tfplot(WT,'log',fs,t); 
    tfplot(WT,'log',fs,t); hold on;   plot(t,X); hold on;  
%     
%    tfplot(WSST2,'log',fs,t);
%    colridge(gca,round(1+nv*log2(Omega))); %set(gca,'ylim',[nv*4 nv*9+1]); 
%       
    
%     tfplot(real(omega),'log',fs,t); 
%     colridge(gca,C1); set(gca,'ylim',[nv*4 nv*9+1]); 

%     figure(); 
%     imagesc(t,fs,real(omega2)); set(gca,'YDir','normal'); title('omega2 WT'); 
%     %hold on; plot(t,if1,'-.',t,if2,'--',t,if3); 
% 
%     tfplot(real(omega2),'log',fs,t); 
%     colridge(gca,C2); set(gca,'ylim',[nv*4 nv*9+1]); 
    
%     
%     figure(); 
%     imagesc(t,fs, 1./(real(phipp)./real(omega2.^2))); set(gca,'YDir','normal'); title('\phi''^2/\phi" ')  
%     
%     figure(); 
%     imagesc(t,fs,real(phipp)); set(gca,'YDir','normal'); title('Frequency Modulation \phi" ')  


end

%% Normalized coefficient reconstruction WSST2
if  strncmp(mywav,'cmor',4)
    [v1 v2] = regexp(mywav,'[0-9.]*-[0-9.]');
    Fb = str2num(mywav(v1:v2-2));
    Fc= str2num(mywav(v2:end));
    filt = @(a) cmor(Fb,Fc,a*xi);
    func = @(x) sqrt(Fb) * exp(-Fb^2*pi*(x-Fc).^2);
    func1 = @(x) conj(sqrt(Fb) * exp(-Fb^2*pi*(x-Fc).^2))./x;
    Cpsi1 = quadgk(func1,0,inf); 
    func2 = cell(nmodes,N);
    Cpsi2 = Cpsi1*ones(nmodes,N); 
    RaPO = zeros(nmodes,N); 
    
%    for q = 1:N
%        phipp(C2(1,q),q) = iff1(q);
%        phipp(C2(2,q),q) = iff3(q);
%        phipp(C2(3,q),q) = iff2(q);
%        omega2(C2(1,q),q) = if1(q);
%        omega2(C2(2,q),q) = if3(q);
%        omega2(C2(3,q),q) = if2(q);
%    end

    
 for p=1:nmodes  
     for q = 1:N
         RaPO(p,q) = real(phipp(C2(p,q),q))./(real(omega2(C2(p,q),q))).^2;
         if ~isnan(RaPO(p,q)) 
             func2{p,q} = @(x) conj(1/sqrt(Fb).*sqrt(1./(1./Fb^2+1i.*RaPO(p,q).*x.^2)).*exp(-pi.*(x-Fc).^2./(1./Fb^2+1i.*RaPO(p,q).*x.^2)))./x;
             Cpsi2(p,q) = quadgk(func2{p,q},0,inf);   
         end
     end
 end     

end


%%
srec1 = (1./Cpsi1).*recmodes(WSST,C1,d); srec1=transpose(srec1);

% disp('SNRout by WSST')
% disp(snr(s1(index),srec1(index,3)- s1(index)));
% disp(snr(s2(index),srec1(index,2)- s2(index)));
% disp(snr(s3(index),srec1(index,1)- s3(index)));
% disp(snr(s(index),sum(srec1(index,:),2)- s(index)));

figure(); 
subplot(3,1,1)
plot(t(index),real(s1(index)),'--',t(index),real(srec1(index,3)));
subplot(3,1,2)
plot(t(index),real(s2(index)),'--',t(index),real(srec1(index,2)));
subplot(3,1,3)
plot(t(index),real(s3(index)),'--',t(index),real(srec1(index,1)));

%% Ridge extraction WSST2
srec2 = (1./Cpsi2).*recmodes(WSST2,C2,d); srec2=transpose(srec2);

% disp('SNRout by WSST2')
% disp(snr(s1(index),srec2(index,3)- s1(index)));
% disp(snr(s2(index),srec2(index,2)- s2(index)));
% disp(snr(s3(index),srec2(index,1)- s3(index)));
% disp(snr(s(index),sum(srec2(index,:),2)- s(index)));
% 
figure() 
subplot(3,1,1)
plot(t(index),real(s1(index)),'--',t(index),real(srec2(index,3)));
subplot(3,1,2)
plot(t(index),real(s2(index)),'--',t(index),real(srec2(index,2)));
subplot(3,1,3)
plot(t(index),real(s3(index)),'--',t(index),real(srec2(index,1)));

%set(findall(0,'type','axes'),'xlim',[0.2 0.8]); 

%% Mode reconstruction accuracy
ac = cell(5,5);  ac{1,2} = {}; 
ac{1,2} = 'FSST'; ac{1,3} = 'FSST2'; ac{1,4} = 'WSST'; ac{1,5} = 'WSST2'; 
ac{2,1} = 'Mode f_1'; ac{3,1} = 'Mode f_2'; ac{4,1} = 'Mode f_3'; ac{5,1} = 'Mode f'; 
%FSST
ac{2,2} = snr(s1(index),imf(index,3)- s1(index));
ac{3,2} = snr(s2(index),imf(index,2)- s2(index));
ac{4,2} = snr(s3(index),imf(index,1)- s3(index));
ac{5,2} = snr(s(index),sum(imf(index,:),2)- s(index));
%FSST2
ac{2,3} = snr(s1(index),imf2(index,3)- s1(index));
ac{3,3} = snr(s2(index),imf2(index,2)- s2(index));
ac{4,3} = snr(s3(index),imf2(index,1)- s3(index));
ac{5,3} = snr(s(index),sum(imf2(index,:),2)- s(index));

%WSST
ac{2,4} = snr(s1(index),srec1(index,3)- s1(index));
ac{3,4} = snr(s2(index),srec1(index,2)- s2(index));
ac{4,4} = snr(s3(index),srec1(index,1)- s3(index));
ac{5,4} = snr(s(index),sum(srec1(index,:),2)- s(index));
%WSST2
ac{2,5} = snr(s1(index),srec2(index,3)- s1(index));
ac{3,5} = snr(s2(index),srec2(index,2)- s2(index));
ac{4,5} = snr(s3(index),srec2(index,1)- s3(index));
ac{5,5} = snr(s(index),sum(srec2(index,:),2)- s(index));

disp(ac)
%digits(3); %this changes the output precision
%ac=; %the ???d??? flag makes sure the sym output is in decimal form
%latex(sym(ac,'d'))
%end
