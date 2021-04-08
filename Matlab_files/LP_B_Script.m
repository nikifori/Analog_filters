close all; 
clear all;

%% AEM
% AEM=9084

a1=9;
a2=0;
a3=8;
a4=4;

%% Prodiagrafes filtrou

mi=1; %a3=8
f_p=(3+mi)*600; %Hz
f_s=2*f_p;      %Hz

amin=17.5+(max(1,a4)-5)*0.5;     %dB
amax=0.6+(max(1,a3)-3)/10; %dB

omegaP=2*pi*f_p;
omegaS=2*pi*f_s;

n_arxiko=1/2*log((10^(amin/10)-1)/(10^(amax/10)-1))/log(omegaS/omegaP);
n=4; % afou n_arxiko=3.706405445004520
omega0=omegaP/(10^(amax/10)-1)^(1/(2*n));

Q1=0.54;
Q2=1.31;

%% Monades filtrou
% ylopoioyme 2 monades sallen key
% 1h Monada me stratigiki 2 logw a2=0
M1_Q=Q1;

M1_k=1; % gia omega0=1
M1_R1=1;
M1_R2=1;
M1_C1=2*M1_Q;
M1_C2=1/(2*M1_Q);

% klimakopoihsh
M1_kf=omega0;
M1_km=M1_C1/(M1_kf*10^(-7));  % thelw 0.1 uF=10^-7

M1_R1_k=M1_km*M1_R1;
M1_R2_k=M1_km*M1_R2;
M1_C1_k=(M1_C1/(M1_kf*M1_km));
M1_C2_k=(M1_C2/(M1_kf*M1_km));

% sinartisi metaforas ths monadas

T1=tf([0 0 M1_k*omega0^2],[1 omega0/M1_Q omega0^2]);

% ------------------------------------------------------------------------------
% 1h Monada me stratigiki 2 logw a2=0
M2_Q=Q2;

M2_k=1; % gia omega0=1
M2_R1=1;
M2_R2=1;
M2_C1=2*M2_Q;
M2_C2=1/(2*M2_Q);

% klimakopoihsh
M2_kf=omega0;
M2_km=M2_C1/(M2_kf*10^(-7));  % thelw 0.1 uF=10^-7

M2_R1_k=M2_km*M2_R1;
M2_R2_k=M2_km*M2_R2;
M2_C1_k=(M2_C1/(M2_kf*M2_km));
M2_C2_k=(M2_C2/(M2_kf*M2_km));

% sinartisi metaforas ths monadas

T2=tf([0 0 M2_k*omega0^2],[1 omega0/M2_Q omega0^2]);
% -----------------------------------------------------------------------------------
% Rithmisi kerdous
% kerdos stis xamhles sixnothtes

Total_gain=M1_k*M2_k; %=1;

% rithmish kerdous sta 10 dB enisxish epeidh a2=0
deka_db=10^(10/20); %=3.162277660168380------>3.2 gia na eimaste sigouroi oti tha exoyme 10 db kerdos
deka_db=3.2;
% rithmish kerdous me telestiko enisxith xwris anastrofi
k_telestikou=deka_db/Total_gain;
T_R1=1000;
T_R2=(k_telestikou-1)*T_R1;

%% Transfer fuction
k3=3.2;
f0=omega0/(2*pi);
sys1=tf([0 omega0^2],[1 omega0/Q1 omega0^2]);
sys2=tf([0 omega0^2],[1 omega0/Q2 omega0^2]);
sysall=k3*series(sys1,sys2);
sysall_inv_0=inv(sysall/k3);
sysall_inv=inv(sysall);
sys1_kanonik=tf([0 1],[1 1/Q1 1]);
sys2_kanonik=tf([0 1],[1 1/Q2 1]);
sysall_kanonik=series(sys1_kanonik,sys2_kanonik);


plot_transfer_function(sys1, [f_p,f0,f_s]);
plot_transfer_function(sys2, [f_p,f0,f_s]);
plot_transfer_function(sysall, [f_p,f0,f_s]);
ltiview({'bode'}, sys1,sys2,sysall);
plot_transfer_function(sysall_inv_0, [f_p,f0,f_s]);
plot_transfer_function(sysall_inv, [f_p,f0,f_s]);
ltiview({'bode'}, sysall_kanonik);

%% Fourier analysis


%Sawtooth input signal
T=5*(1/2000);
fs=50000;
t=0:1/fs:0.005;


in=sawtooth(2*pi*2000*t);
in=(in+1)/1.95;
figure;
plot(t,in);
title('Sawtooth input signal Matlab');
xlabel('Time (s)');
ylabel('Voltage (V)');

% fourier
% single-sided spectrum of input signal 
infour=fft(in);
Lin=length(in)+1;
P2=abs(infour/Lin);
P1=P2(1:Lin/2+1);
P1(2:end-1)=2*P1(2:end-1);
f = fs*(0:(Lin/2))/Lin;

figure;
plot(f,P1);
title('Sawtooth input fourier analysis Matlab (Single-sided)');
xlabel('Frequency (Hz)');

% output signal 
out=lsim(sysall,in,t);
figure;
plot(t,out);
title('Output signal of sawtooth input Matlab');
xlabel('Time (s)');
ylabel('Voltage (V)');

% in-out signal plots
figure;
plot(t,in);
hold on;
plot(t,out);
hold off;
title('Input and output of sawtooth input Matlab');
xlabel('Time (s)');
ylabel('Voltage (V)');

% fourier
% single-sided spectrum of output signal
outfour=fft(out);
Lout=length(out)+1;
P4=abs(outfour/Lout);
P3=P4(1:Lout/2+1);
P3(2:end-1)=2*P3(2:end-1);
f = fs*(0:(Lout/2))/Lout;

figure;
plot(f,P3);
title('Output of sawtooth input fourier analysis Matlab (Single-sided)');
xlabel('Frequency (Hz)');
