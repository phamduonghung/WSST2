function [a1,a2,a3,if1,if2,if3,s1,s2,s3,s,iff1,iff2,iff3] = signal_test(t,n0)

if n0==1 % Gaussian amplitudes
    a1 =  exp(-5*(t-0.5).^2);
    a2 =  exp(-7*(t-0.4).^2);
    a3 = exp(-10*(t-0.6).^2);
    
    cof1 = 47;
    cof2 =1.1;
    cof3 = 1.77;

    phi1 =28*t+cof1*t.^2;
    phi2 = -83*log(cof2-t);
    phi3 = 62.*exp(cof3.*t);

    if1 = 28*ones(length(t),1)+cof1*2*t;
    if2 = 83./(cof2-t);
    if3 = 62*cof3*exp(cof3.*t);

    iff1 = cof1*2*ones(length(t),1);
    iff2 = 83./(cof2-t).^2;
    iff3 = 62*cof3*cof3*exp(cof3.*t);
elseif n0==2 % constant amplitudes
    a1 = 1;
    a2 = 1;
    a3 = 1;
    
    cof1 = 47;
    cof2 =1.1;
    cof3 = 1.77;

    phi1 =28*t+cof1*t.^2;
    phi2 = -83*log(cof2-t);
    phi3 = 62.*exp(cof3.*t);

    if1 = 28*ones(length(t),1)+cof1*2*t;
    if2 = 83./(cof2-t);
    if3 = 62*cof3*exp(cof3.*t);

    iff1 = cof1*2*ones(length(t),1);
    iff2 = 83./(cof2-t).^2;
    iff3 = 62*cof3*cof3*exp(cof3.*t);
else  % constant chirps
    a1 = 1;
    a2 = 1;
    a3 = 1;
    
    cof1 = 0;
    
    phi1 =60.5*t+cof1*t.^2;
    phi2 = 200*t+cof1*t.^2;
    phi3 = 350*t+cof1*t.^2;

    if1 = 60.5*ones(length(t),1)+cof1*2*t;
    if2 = 200*ones(length(t),1)+cof1*2*t;
    if3 = 350*ones(length(t),1)+cof1*2*t;

    iff1 = cof1*2*ones(length(t),1);
    iff2 = cof1*2*ones(length(t),1);
    iff3 = cof1*2*ones(length(t),1);
end


% 
% phi1 = 50*t+30*t.^3-20*(1-t).^4;
% phi2 = 340*t-2.*exp(-2*(t-0.2)).*sin(14*pi.*(t-0.2));
% phi3 = 200*t;
% if1 = 50+90*t.^2+80*(1-t).^3; 
% if2 = 340+4*exp(-2*(t-0.2)).*sin(14*pi.*(t-0.2))-28*pi.*exp(-2*(t-0.2)).*cos(14*pi.*(t-0.2)); 
% if3 = 200*ones(length(t),1);

% phi1 = 60*t-100*log(1.05-t);
% phi2 = 70*t+70*t.^2;
% phi3 = 30*t;
% 
% if1 = 60+100./(1.05-t); 
% if2 = 70+140*t; 
% if3 = 30*ones(length(t),1); 

% phi1 = 25*exp(3*t);
% phi2 = 50*t+70*t.^2;
% phi3 = 40*t;
% 
% if1 = 75*exp(3*t); 
% if2 = 50+140*t; 
% if3 = 40*ones(length(t),1); 

s1 = a1.*exp(2*pi*1i*(phi1));
s2 = a2.*exp(2*pi*1i*(phi2));
s3 = a3.*exp(2*pi*1i*(phi3));

s = s1+s2+s3;

end
