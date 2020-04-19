function [STFT,SST,VSST,omega,tau,omega2,g] = sst2_new(s,aa,Nfft,gamma)

 %% sst2_new : computes the STFT of a signal and different versions of synchrosqueezing
 %
 % INPUTS:   
 %   s: real or complex signal
 %   aa: the parameter a in the amgauss function of the Gaussian window
 %   Nfft: number of frequency bins
 %   gamma: threshold on the STFT for reassignment
 % OUTPUTS:   
 %   STFT : the short-time Fourier transform
 %   SST  : standard synchrosqueezing
 %   VSST : vertical second-order synchrosqueezing 
 % REFERENCES
 % [1] Pham, D-H., Meignen, S. Second-order Synchrosqueezing Transform: The
 % Wavelet Case and Comparisons
 
 s = s(:);
           
 ft   = 1:Nfft;
 bt   = 1:length(s);
 nb   = length(bt);
 neta = length(ft);
        
 prec = 10^(-6) ;
 L = sqrt(Nfft/aa);
 l = floor(sqrt(-Nfft*log(prec)/(aa*pi)))+1;
 g = amgauss(2*l+1,l+1,L);       
  %figure(); plot(g)
 % Window definition
 n   = (0:2*l)'-l;
 t0  = n/Nfft;
 t0  = t0(:);
 a   = aa*Nfft*pi; 
 gp  = -2*a*t0.*g; 
 gpp = (-2*a+4*a^2*t0.^2).*g; % g''
  
 % Initialization
 STFT  = zeros(neta,nb);
 SST   = zeros(neta,nb);
 VSST  = zeros(neta,nb);
 
 omega  = zeros(neta,nb);
 tau    = zeros(neta,nb);
 omega2 = zeros(neta,nb);
 phipp  = zeros(neta,nb);
             
 %% Computes STFT and reassignment operators
        
 for b=1:nb
 	% STFT, window g  
 	time_inst = -min([l,bt(b)-1]):min([l,nb-bt(b)]); 
    tmp(1:length(time_inst)) = s(bt(b)+time_inst).*g(l+1+time_inst);
    A = fft(tmp(:),Nfft);
    STFT(:,b) = A.*exp(-2/Nfft*pi*1i*(0:Nfft-1)'*time_inst(1))/Nfft*g(l+1);%renormalized so that it fits with recmodes
 	vg  = A;
     
 	% STFT, window xg
    tmp(1:length(time_inst)) = s(bt(b)+time_inst).*(time_inst)'/Nfft.*g(l+1+time_inst);
    vxg = fft(tmp(:),Nfft);
  
    % operator Lx (dtau)
	tau(:,b)  = vxg./vg;
 	
    % STFT, window gp
    tmp(1:length(time_inst))= s(bt(b)+time_inst).*gp(l+1+time_inst);
    vgp = fft(tmp(:),Nfft);
 
            
    omega(:,b) = (ft-1)'- real(vgp/2/1i/pi./vg);
 	
    % STFT, window gpp
 	tmp(1:length(time_inst)) = s(bt(b)+time_inst).*gpp(l+1+time_inst);
    vgpp = fft(tmp(:),Nfft);
 
    %STFT, windox xgp
 	tmp(1:length(time_inst)) = s(bt(b)+time_inst).*(time_inst)'/Nfft.*gp(l+1+time_inst);
 	vxgp = fft(tmp(:),Nfft);
      
 	%computation of the two different omega 
        
    phipp(:,b) = 1/2/1i/pi*(vgpp.*vg-vgp.^2)./(vxg.*vgp-vxgp.*vg);
       
    %new omega2
    %omega2(:,b) = omega(:,b) - real(phipp(:,b)).*real(tau(:,b))...
    %                          + imag(phipp(:,b)).*imag(tau(:,b)); 
    %omega2(:,b) = omega(:,b) - real(phipp(:,b)).*real(tau(:,b)); 

 end

 %% reassignment step
 phipp(abs(real(tau)*nb)<1)=0;
 %omega2 = real(omega - phipp.*tau);
 omega2 = real(omega) - real(phipp).*real(tau);

 for b=1:nb
   
     for eta=1:neta
        if (abs(STFT(eta,b))> 2*sqrt(2)*gamma*norm(g)/length(s))
            k = 1+round(omega(eta,b));
            if (k >= 1) && (k <= neta)
                % original reassignment
                SST(k,b) = SST(k,b) + STFT(eta,b);
            end
           
            k = 1+round(omega2(eta,b));
            if k>=1 && k<=neta
                % second-order Vertical reassignment: VSST
                VSST(k,b) = VSST(k,b) + STFT(eta,b);
            end 
        end
    end
 end
end

