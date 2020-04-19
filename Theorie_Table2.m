%function accuracy_3modes()

clear all; close all;
set(0,'DefaultAxesFontSize',14);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Theorie_WSST2_2017/elsarticle/figures';


% Parameters
N = 1024;
%gamma = 0.001;
sigma = 0.05;
%index = 2*N/8+1:6*N/8;
index = round(0.2*N):round(0.8*N);

d = 0;
clwin = 10;
lambda = 0.001; beta = 0;
nv = 32;
mywav = 'cmor6-1';
nmodes = 3;

t = (0:N-1)/N; t=t(:);
frecv = 1:N/2;
% Choice of time and frequency bins
ft =1:N/2;bt=1:N;

%% Test signal 
[a1,a2,a3,if1,if2,if3,s1,s2,s3,s] = signal_test(t,1);

%% noise
SNR = 0;
noise = randn(N,1)+1i*randn(N,1);
sb1 = sigmerge(s1,noise,SNR);
noise = randn(N,1)+1i*randn(N,1);
sb2 = sigmerge(s2,noise,SNR);
noise = randn(N,1)+1i*randn(N,1);
sb3 = sigmerge(s3,noise,SNR);

sb = sb1+sb2+sb3;
disp([num2str(round(20*log10(norm(s)/norm(sb-s(:))))) 'dB for input SNR']);

gamma = std(sb-s)

%% TF presentations of signal
%[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3] = sstn(s,gamma,sigma,ft,bt);
[STFT,FSST,FSST2,omegaF,tauF,omega2F] = sst2_new(sb,1/sigma^2/N,N,gamma);
% STFT = STFT(1:N/2,:);
% FSST = FSST(1:N/2,:);
% FSST2 = FSST2(1:N/2,:);

[WT, WSST, WSST2, fs, as, omega, tau, phipp, ~]  = Wsst2_new(sb,gamma,mywav,nv);
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
% disp(snr(s(index),imf(index,3)+imf(index,2)+imf(index,1)- s(index)));

figure() 
subplot(3,1,1)
plot(t(index),real(s1(index)),'--',t(index),real(imf(index,3)))
subplot(3,1,2)
plot(t(index),real(s2(index)),'--',t(index),real(imf(index,2)))
subplot(3,1,3)
plot(t(index),real(s3(index)),'--',t(index),real(imf(index,1)))

%% Ridge extraction FSST2

[Cs2, Es2] = exridge_mult_Noise(FSST2, nmodes, lambda, clwin);
%[Cs2, Es2] = brevridge_mult(FSST2, frecv, nmodes, lambda, clwin);

imf2 = recmodes(FSST2,Cs2,d); imf2=transpose(imf2);

% disp('SNRout by FSST2')
% disp(snr(s1(index),imf2(index,3)- s1(index)));
% disp(snr(s2(index),imf2(index,2)- s2(index)));
% disp(snr(s3(index),imf2(index,1)- s3(index)));
% disp(snr(s(index),imf2(index,3)+imf2(index,2)+imf2(index,1)- s(index)));

figure() 
subplot(3,1,1)
plot(t(index),real(s1(index)),'--',t(index),real(imf2(index,3)))
subplot(3,1,2)
plot(t(index),real(s2(index)),'--',t(index),real(imf2(index,2)))
subplot(3,1,3)
plot(t(index),real(s3(index)),'--',t(index),real(imf2(index,1)))

%% Normalized coefficient reconstruction WSST2
%[C1, Es3] = brevridge_mult(WSST, frecv, nmodes, lambda, clwin);
%[C2, Es4] = brevridge_mult(WSST2, frecv, nmodes, lambda, clwin);
[C1, Es] = exridge_mult_Noise(WSST, nmodes, lambda, clwin);
[C2, Es] = exridge_mult_Noise(WSST2, nmodes, lambda, clwin);


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

%% Ridge extraction WSST
srec1 = 1./Cpsi1.*recmodes(WSST,C1,d); srec1=transpose(srec1);

% disp('SNRout by WSST')
% disp(snr(s1(index),srec1(index,3)- s1(index)));
% disp(snr(s2(index),srec1(index,2)- s2(index)));
% disp(snr(s3(index),srec1(index,1)- s3(index)));
% disp(snr(s(index),srec1(index,3)+srec1(index,2)+srec1(index,1)- s(index)));

figure() 
subplot(3,1,1)
plot(t(index),real(s1(index)),'--',t(index),real(srec1(index,3)))
subplot(3,1,2)
plot(t(index),real(s2(index)),'--',t(index),real(srec1(index,2)))
subplot(3,1,3)
plot(t(index),real(s3(index)),'--',t(index),real(srec1(index,1)))
%% Ridge extraction WSST2
srec2 = 1./Cpsi2.*recmodes(WSST2,C2,d); srec2=transpose(srec2);

% disp('SNRout by WSST2')
% disp(snr(s1(index),srec2(index,3)- s1(index)));
% disp(snr(s2(index),srec2(index,2)- s2(index)));
% disp(snr(s3(index),srec2(index,1)- s3(index)));
% disp(snr(s(index),srec2(index,3)+srec2(index,2)+srec2(index,1)- s(index)));

figure() 
subplot(3,1,1)
plot(t(index),real(s1(index)),'--',t(index),real(srec2(index,3)))
subplot(3,1,2)
plot(t(index),real(s2(index)),'--',t(index),real(srec2(index,2)))
subplot(3,1,3)
plot(t(index),real(s3(index)),'--',t(index),real(srec2(index,1)))

%pause
%% Mode reconstruction accuracy
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

% ac = zeros(4,4);
% %FSST
% ac(1,1) = snr(s1(index),imf(index,3)- s1(index));
% ac(2,1) = snr(s2(index),imf(index,2)- s2(index));
% ac(3,1) = snr(s3(index),imf(index,1)- s3(index));
% ac(4,1) = snr(s(index),imf(index,3)+imf(index,2)+imf(index,1)- s(index));
% %FSST2
% ac(1,2) = snr(s1(index),imf2(index,3)- s1(index));
% ac(2,2) = snr(s2(index),imf2(index,2)- s2(index));
% ac(3,2) = snr(s3(index),imf2(index,1)- s3(index));
% ac(4,2) = snr(s(index),imf2(index,3)+imf2(index,2)+imf2(index,1)- s(index));
% 
% %WSST
% ac(1,3) = snr(s1(index),srec1(index,3)- s1(index));
% ac(2,3) = snr(s2(index),srec1(index,2)- s2(index));
% ac(3,3) = snr(s3(index),srec1(index,1)- s3(index));
% ac(4,3) = snr(s(index),srec1(index,3)+srec1(index,2)+srec1(index,1)- s(index));
% %WSST2
% ac(1,4) = snr(s1(index),srec2(index,3)- s1(index));
% ac(2,4) = snr(s2(index),srec2(index,2)- s2(index));
% ac(3,4) = snr(s3(index),srec2(index,1)- s3(index));
% ac(4,4) = snr(s(index),srec2(index,3)+srec2(index,2)+srec2(index,1)- s(index));
% 
% digits(3); %this changes the output precision
% %ac=; %the ???d??? flag makes sure the sym output is in decimal form
% latex(sym(ac,'d'))
% %end
