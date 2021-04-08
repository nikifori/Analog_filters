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
fp=(3+mi)*1000; %Hz
fs=fp/1.8;      %Hz

amin=25+a3*4/9;     %dB
amax=0.5+a4*0.25/9; %dB

omegap=2*pi*fp;
omegas=2*pi*fs;

%% Ypologismoi

% Prototipo katwdiavato
Omega_p=1;
Omega_s=omegap/omegas;

n_arxiko=acosh(sqrt((10^(amin/10)-1)/(10^(amax/10)-1)))/acosh(Omega_s);
n=5; % afou n_arxiko=4.128473397115852

epsilon=sqrt(10^(amax/10)-1);
alpha=(1/n)*asinh(1/epsilon);
Omega_hp=cosh(acosh(1/epsilon)/n);

%----------------------------------------------------------------------------
% Gwnies twn polwn Butterworth gia n=5
psi1=0;
psi23=36;    %+-
psi45=72;    %+-

%Poloi gia LP Chebysev n=5
% 0
sigma1=-sinh(alpha)*cosd(psi1);
pole_omega1=cosh(alpha)*sind(psi1);   %+- =0
% 36
sigma23=-sinh(alpha)*cosd(psi23);
pole_omega23=cosh(alpha)*sind(psi23);   %+-
% 72
sigma45=-sinh(alpha)*cosd(psi45);
pole_omega45=cosh(alpha)*sind(psi45);   %+-

% (wo,Q) gia Chebysev
w0_1=sqrt(sigma1^2+pole_omega1^2);
Q_1=1/(2*cos(atan(pole_omega1/sigma1)));

w0_23=sqrt(sigma23^2+pole_omega23^2);
Q_23=1/(2*cos(atan(pole_omega23/sigma23)));

w0_45=sqrt(sigma45^2+pole_omega45^2);
Q_45=1/(2*cos(atan(pole_omega45/sigma45)));
%--------------------------------------------------------------------------------
% antistrofi twn polwn 
% sixnotita hmiseias isxios HighPass
omega_hp_highpass=omegap/Omega_hp;

% poloi highpass
P1_S1=omegap/abs(sigma1);   %rad/s
P23_w23=omegap/w0_23;       %rad/s
P45_w45=omegap/w0_45;        %rad/s

%% Ylopoihsh ths sinartishs metaforas

% 1h monada CR tou pragmatikou polou
% P1_S1

M1_R1=1;
M1_C1=1;
M1_k=1; % kerdos monada 1

% klimakopoihsh
M1_kf=P1_S1;
M1_km=M1_C1/(M1_kf*10^(-6));  % thelw 1.0 uF=10^-6

M1_R1_k=M1_km*M1_R1;
M1_C1_k=(M1_C1/(M1_kf*M1_km));

% sinartisi metaforas ths monadas

T1=tf([0 M1_k 0],[0 1 P1_S1]);

%-----------------------------------------------------------------------------------
% 2h monada Sallen-Key tou P23_w23 tou 
% sxhmatos 6.6 me allagh theseis twn C-R
% Stratigikh 2 logw a3=8

M2_k=1; %kerdos 1 logw stratigikhs 2
M2_Q=Q_23;

M2_C1=1;
M2_C2=1;
M2_R2=2*M2_Q;
M2_R1=1/(2*M2_Q);

% klimakopoihsh
M2_kf=P23_w23;
M2_km=M2_C1/(M2_kf*10^(-6));  % thelw 1.0 uF=10^-6

M2_R1_k=M2_km*M2_R1;
M2_R2_k=M2_km*M2_R2;
M2_C1_k=(M2_C1/(M2_kf*M2_km));
M2_C2_k=(M2_C2/(M2_kf*M2_km));

% sinartisi metaforas ths monadas

T2=tf([M2_k 0 0],[1 P23_w23/M2_Q P23_w23^2]);

% -------------------------------------------------------------------------
% 3h monada Sallen-Key tou P45_w45 tou 
% sxhmatos 6.6 me allagh theseis twn C-R
% Stratigikh 2 logw a3=8

M3_k=1; %kerdos 1 logw stratigikhs 2
M3_Q=Q_45;

M3_C1=1;
M3_C2=1;
M3_R2=2*M3_Q;
M3_R1=1/(2*M3_Q);

% klimakopoihsh
M3_kf=P45_w45;
M3_km=M3_C1/(M3_kf*10^(-6));  % thelw 1.0 uF=10^-6

M3_R1_k=M3_km*M3_R1;
M3_R2_k=M3_km*M3_R2;
M3_C1_k=(M3_C1/(M3_kf*M3_km));
M3_C2_k=(M3_C2/(M3_kf*M3_km));

% sinartisi metaforas ths monadas

T3=tf([M3_k 0 0],[1 P45_w45/M3_Q P45_w45^2]);

% -------------------------------------------------------------------------

% Rithmish kerdous
% kerdos stis ipsiles sixnothtes

Total_gain=M1_k*M2_k*M3_k; %=1;

% rithmish kerdous sta 10 dB enisxish epeidh a2=0
deka_db_=10^(10.611111111111111/20); % vazoyme kati parapanw apo 10 db dioti stis ipsiles
% sixnotites epeidh einai Chebyshev tha exei diakimansi amax=0.6111
% kai afou theloume toulaxiston 10 db se oles tis ipsiles sixnothtes tha
% epilejoume enisxisi sta 10.6111 db, wste se kamia ipsili sixnothta na mhn
% peftw katw apo 10 db.
deka_db=3.392774536572461;
% rithmish kerdous me telestiko enisxith xwris anastrofi
k_telestikou=deka_db/Total_gain;
T_R1=1000;
T_R2=(k_telestikou-1)*T_R1;

% ---------------------------------------------------------------------------

% Sinartiseis metaforas
% gia thn klimakopoihsh
M1_omega0=P1_S1/omega_hp_highpass;
M2_omega0=P23_w23/omega_hp_highpass;
M3_omega0=P45_w45/omega_hp_highpass;


T12=series(T1,T2);
T_all=series(T12,T3);
T_hp=k_telestikou*T_all;
invT_all=inv(T_all);
invT_hp=inv(T_hp);
T1_kanonik=tf([0 M1_k 0],[0 1 M1_omega0]);
T2_kanonik=tf([M2_k 0 0],[1 M2_omega0/M2_Q M2_omega0^2]);
T3_kanonik=tf([M3_k 0 0],[1 M3_omega0/M3_Q M3_omega0^2]);
T12_kanonik=series(T1_kanonik,T2_kanonik);
T_kanonik=series(T12_kanonik,T3_kanonik);

plot_transfer_function(T1, [fs fp 10^5]);
plot_transfer_function(T2, [fs fp 10^5]);
plot_transfer_function(T3, [fs fp 10^5]);
plot_transfer_function(T_all, [fs fp 10^5]);
plot_transfer_function(T_hp, [fs fp 10^5]);
ltiview({'bode'}, T1,T2,T3,T_hp)
plot_transfer_function(invT_all, [fs fp 10^5]);
plot_transfer_function(invT_hp, [fs fp 10^5]);
ltiview({'bode'}, T_kanonik);




%% Fourier analysis

% a4=4


in_freq1=0.5*omegas/(2*pi);
in_freq2=0.8*omegas/(2*pi);
in_freq3=1.2*omegap/(2*pi);
in_freq4=2.4*omegap/(2*pi);
in_freq5=3.5*omegap/(2*pi);


T=10*(1/100);  % mporw na valw kai prscn me (1/100) pou fainontai kalitera alla den exw megali akriveia.
             
f_s=1000000; % sixnothta deigmatolipsias
dt=1/f_s;
t=0:dt:T-dt;

in=cos(0.5*omegas*t)+0.6*cos(0.8*omegas*t) ...
+cos(1.2*omegap*t)+0.8*cos(2.4*omegap*t) ...
+0.4*cos(3.5*omegap*t);

figure;
plot(t,in);
title('Input signal Matlab with multiple frequencies');
xlabel('Time (s)');
ylabel('Voltage (V)');
axis([0 0.01 -14 14]);


% fourier
% single-sided spectrum of input signal 
infour=fft(in);
Lin=length(in);
P2=abs(infour/Lin);
P1=P2(1:Lin/2+1);
P1(2:end-1)=2*P1(2:end-1);
f = f_s*(0:(Lin/2))/Lin;

figure;
plot(f,P1);
title('Multiple frequencies input fourier analysis Matlab (Single-sided)');
xlabel('Frequency (Hz)');
axis([0 15000 0 2]);

% output signal 
out=lsim(T_hp,in,t);
figure;
plot(t,out);
title('Output signal of multiple frequencies input Matlab');
xlabel('Time (s)');
ylabel('Voltage (V)');
axis([0 0.01 -14 14]);

% in-out signal plots
figure;
plot(t,in);
hold on;
plot(t,out);
hold off;
title('Input and output of multiple frequencies input Matlab');
xlabel('Time (s)');
ylabel('Voltage (V)');
axis([0 0.01 -14 14]);

% fourier
% single-sided spectrum of output signal
outfour=fft(out);
Lout=length(out);
P4=abs(outfour/Lout);
P3=P4(1:Lout/2+1);
P3(2:end-1)=2*P3(2:end-1);
f = f_s*(0:(Lout/2))/Lout;

figure;
plot(f,P3);
title('Output of multiple frequencies input fourier analysis Matlab (Single-sided)');
xlabel('Frequency (Hz)');
axis([0 15000 0 4]);













