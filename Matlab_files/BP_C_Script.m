close all; 
clear all;

%% AEM
% AEM=9084

a1=9;
a2=0;
a3=8;
a4=4;

%% Prodiagrafes filtrou

f0=650; %Hz
f1=400+25*a3;
f2=(f0^2)/f1;
D=2.3*((f0^2-f1^2)/f1);
f3=(-D+sqrt(D^2+4*f0^2))/2;
f4=f0^2/f3;

amin=27.5+a4; %dB
amax=0.5+(a3-5)/10; %dB

omega0=2*pi*f0;  %rad/s
omega1=2*pi*f1;
omega2=2*pi*f2;
omega3=2*pi*f3;
omega4=2*pi*f4;

bw=omega2-omega1;
qc=omega0/bw;
%% Ypologismoi

Omega_p=1;
Omega_s=(omega4-omega3)/(omega2-omega1);

n_arxiko=acosh(sqrt((10^(amin/10)-1)/(10^(amax/10)-1)))/acosh(Omega_s);
%tipos 9.83 gia thn taxh Chebysev 

n=4;  %afou n_arxiko=3.470008857513355

epsilon=sqrt(10^(amax/10)-1);
alpha=(1/n)*asinh(1/epsilon);

Omega_hp=cosh(acosh(1/epsilon)/n); %sixnotita hmisias isxios
% --------------------------------------------------------------------------

%gwnies Butterworth me n=4

psi12=22.5;   %+-              % polosX=sigmaX+j*pole_omegaX
psi34=67.5;   %+-         


%Poloi gia LP Chebysev n=4
%+-22.5
sigma12=-sinh(alpha)*cosd(psi12);
pole_omega12=cosh(alpha)*sind(psi12);   %+-

%+-67.5
sigma34=-sinh(alpha)*cosd(psi34);
pole_omega34=cosh(alpha)*sind(psi34);   %+-

% (wo,Q) gia Chebysev
w0_12=sqrt(sigma12^2+pole_omega12^2);
Q_12=1/(2*cos(atan(pole_omega12/sigma12)));

w0_34=sqrt(sigma34^2+pole_omega34^2);
Q_34=1/(2*cos(atan(pole_omega34/sigma34)));

% metasxhmatismos 1ou Migadikou Polou(1)
P1_S12=abs(sigma12);
P1_W12=abs(pole_omega12);

P1_C=P1_S12^2+P1_W12^2;
P1_D=2*P1_S12/qc;
P1_E=4+P1_C/(qc^2);
P1_G=sqrt(P1_E^2-4*P1_D^2);
P1_Q=(1/P1_D)*sqrt((1/2)*(P1_E+P1_G));
P1_k=P1_S12*P1_Q/qc;
P1_W=P1_k+sqrt(P1_k^2 - 1);
P1_w02=P1_W*omega0;
P1_w01=omega0/P1_W;

P1_psi=acos(1/(2*P1_Q));  %rad

% o prwtos migadikos polos metasximatizetai
% se 2 zeugoi migadikwn polwn 

% 1o zeugos
P1_sigma1=P1_w01/(2*P1_Q);
P1_omega1=P1_w01*sin(P1_psi); %+-

% 2o zeugos
P1_sigma2=P1_w02/(2*P1_Q);
P1_omega2=P1_w02*sin(P1_psi); %+-
% kai pairnw kai alla 2 midenika sto s=0
% --------------------------------------------------------

% metasxhmatismos 2ou Migadikou Polou(2)
P2_S34=abs(sigma34);
P2_W34=abs(pole_omega34);

P2_C=P2_S34^2+P2_W34^2;
P2_D=2*P2_S34/qc;
P2_E=4+P2_C/(qc^2);
P2_G=sqrt(P2_E^2-4*P2_D^2);
P2_Q=(1/P2_D)*sqrt((1/2)*(P2_E+P2_G));
P2_k=P2_S34*P2_Q/qc;
P2_W=P2_k+sqrt(P2_k^2 - 1);
P2_w02=P2_W*omega0;
P2_w01=omega0/P2_W;

P2_psi=acos(1/(2*P2_Q));

% o deuteros migadikos polos metasximatizetai
% se 2 zeugoi migadikwn polwn 

% 1o zeugos
P2_sigma1=P2_w01/(2*P2_Q);
P2_omega1=P2_w01*sin(P2_psi); %+-

% 2o zeugos
P2_sigma2=P2_w02/(2*P2_Q);
P2_omega2=P2_w02*sin(P2_psi); %+-
% kai pairnw kai alla 2 midenika sto s=0

%% Ylopoihsh ths synarthshs metaforas
% ylopoioyme 4 vathmides Delyiannis-Fried 
% stratigikis (1)

% 1h Monada P2_w01
M1_w0=P2_w01;
M1_wz=0;
% klimakopoiw
M1_W0=1;
M1_omega0=omega0/M1_w0;
% sinexizw
M1_Q=P2_Q;
M1_R1=1;
M1_R2=4*M1_Q^2;
M1_C=1/(2*M1_Q);
M1_C1=M1_C;
M1_C2=M1_C;

% sinartisi metaforas tou kanonikopoihmenou kiklwmatos
T1_k=tf([0 -2*M1_Q 0],[1 1/M1_Q 1]);

M1_z=j*M1_omega0;
M1_k_w0=abs(evalfr(T1_k,M1_z)); % vriskw kerdos sto M1_omega0 
% kai thelw na kanw aposvesi ara
% antikathistw thn R1 me thn Z2 kai thn Z3 (diaireths tashs) 
M1_Z2=M1_k_w0;                % wste na exw kerdos 0 dB                          
M1_Z3=M1_k_w0/(M1_k_w0-1);   % se ayth thn monada  sto omega0 apo tis prodiagrafes
 
% klimakopoiw
M1_kf=M1_w0;
M1_km=M1_C1/(M1_kf*10^(-8));  % thelw 0.01 uF=10^-8

M1_R1_k=M1_km*M1_R1;
M1_R2_k=M1_km*M1_R2;
M1_C1_k=(M1_C1/(M1_kf*M1_km));
M1_C2_k=(M1_C2/(M1_kf*M1_km));
M1_Z2_k=M1_km*M1_Z2;
M1_Z3_k=M1_km*M1_Z3;

% sinartisi metaforas ths monadas

T1=tf([0 -(2*M1_Q^2/M1_k_w0)*(M1_w0/M1_Q) 0],[1 (M1_w0/M1_Q) M1_w0^2]);
%plot_transfer_function(T1,[f0 f1 f2 f3 f4]);

% ---------------------------------------------------------------------------

% 2h Monada P1_w01
M2_w0=P1_w01;
M2_wz=0;
% klimakopoiw
M2_W0=1;
M2_omega0=omega0/M2_w0;
% sinexizw
M2_Q=P1_Q;
M2_R1=1;
M2_R2=4*M2_Q^2;
M2_C=1/(2*M2_Q);
M2_C1=M2_C;
M2_C2=M2_C;

% sinartisi metaforas tou kanonikopoihmenou kiklwmatos
T2_k=tf([0 -2*M2_Q 0],[1 1/M2_Q 1]);

M2_z=j*M2_omega0;
M2_k_w0=abs(evalfr(T2_k,M2_z)); % vriskw kerdos sto M2_omega0 
% kai thelw na kanw aposvesi ara
% antikathistw thn R1 me thn Z2 kai thn Z3 (diaireths tashs) 
M2_Z2=M2_k_w0;                % wste na exw kerdos 0 dB                          
M2_Z3=M2_k_w0/(M2_k_w0-1);   % se ayth thn monada  sto omega0 apo tis prodiagrafes

% klimakopoiw
M2_kf=M2_w0;
M2_km=M2_C1/(M2_kf*10^(-8));  % thelw 0.01 uF=10^-8

M2_R1_k=M2_km*M2_R1;
M2_R2_k=M2_km*M2_R2;
M2_C1_k=(M2_C1/(M2_kf*M2_km));
M2_C2_k=(M2_C2/(M2_kf*M2_km));
M2_Z2_k=M2_km*M2_Z2;
M2_Z3_k=M2_km*M2_Z3;

% sinartisi metaforas ths monadas

T2=tf([0 -(2*M2_Q^2/M2_k_w0)*(M2_w0/M2_Q) 0],[1 (M2_w0/M2_Q) M2_w0^2]);
%plot_transfer_function(T2,[f0 f1 f2 f3 f4]);

% ---------------------------------------------------------------------------------

% 3h Monada P1_w02
M3_w0=P1_w02;
M3_wz=0;
% klimakopoiw
M3_W0=1;
M3_omega0=omega0/M3_w0;
% sinexizw
M3_Q=P1_Q;
M3_R1=1;
M3_R2=4*M3_Q^2;
M3_C=1/(2*M3_Q);
M3_C1=M3_C;
M3_C2=M3_C;

% sinartisi metaforas tou kanonikopoihmenou kiklwmatos
T3_k=tf([0 -2*M3_Q 0],[1 1/M3_Q 1]);

M3_z=j*M3_omega0;
M3_k_w0=abs(evalfr(T3_k,M3_z)); % vriskw kerdos sto M3_omega0 
% kai thelw na kanw aposvesi ara
% antikathistw thn R1 me thn Z2 kai thn Z3 (diaireths tashs) 
M3_Z2=M3_k_w0;                % wste na exw kerdos 0 dB                          
M3_Z3=M3_k_w0/(M3_k_w0-1);   % se ayth thn monada  sto omega0 apo tis prodiagrafes   

% klimakopoiw
M3_kf=M3_w0;
M3_km=M3_C1/(M3_kf*10^(-8));  % thelw 0.01 uF=10^-8

M3_R1_k=M3_km*M3_R1;
M3_R2_k=M3_km*M3_R2;
M3_C1_k=(M3_C1/(M3_kf*M3_km));
M3_C2_k=(M3_C2/(M3_kf*M3_km));
M3_Z2_k=M3_km*M3_Z2;
M3_Z3_k=M3_km*M3_Z3;

% sinartisi metaforas ths monadas

T3=tf([0 -(2*M3_Q^2/M3_k_w0)*(M3_w0/M3_Q) 0],[1 (M3_w0/M3_Q) M3_w0^2]);
%plot_transfer_function(T3,[f0 f1 f2 f3 f4]);

% ---------------------------------------------------------------------------

% 4h Monada P2_w01
M4_w0=P2_w02;
M4_wz=0;
% klimakopoiw
M4_W0=1;
M4_omega0=omega0/M4_w0;
% sinexizw
M4_Q=P2_Q;
M4_R1=1;
M4_R2=4*M4_Q^2;
M4_C=1/(2*M4_Q);
M4_C1=M4_C;
M4_C2=M4_C;

% sinartisi metaforas tou kanonikopoihmenou kiklwmatos
T4_k=tf([0 -2*M4_Q 0],[1 1/M4_Q 1]);

M4_z=j*M4_omega0;
M4_k_w0=abs(evalfr(T4_k,M4_z)); % vriskw kerdos sto M4_omega0 
% kai thelw na kanw aposvesi ara
% antikathistw thn R1 me thn Z2 kai thn Z3 (diaireths tashs) 
M4_Z2=M4_k_w0;                % wste na exw kerdos 0 dB                          
M4_Z3=M4_k_w0/(M4_k_w0-1);   % se ayth thn monada  sto omega0 apo tis prodiagrafes   

% klimakopoiw
M4_kf=M4_w0;
M4_km=M4_C1/(M4_kf*10^(-8));  % thelw 0.01 uF=10^-8

M4_R1_k=M4_km*M4_R1;
M4_R2_k=M4_km*M4_R2;
M4_C1_k=(M4_C1/(M4_kf*M4_km));
M4_C2_k=(M4_C2/(M4_kf*M4_km));
M4_Z2_k=M4_km*M4_Z2;
M4_Z3_k=M4_km*M4_Z3;

% sinartisi metaforas ths monadas

T4=tf([0 -(2*M4_Q^2/M4_k_w0)*(M4_w0/M4_Q) 0],[1 (M4_w0/M4_Q) M4_w0^2]);
%plot_transfer_function(T4,[f0 f1 f2 f3 f4]);

% ---------------------------------------------------------------------------
% sinoliko kerdos ths sinolikhs sinartisis metaforas giro apo thn omega0


Total_gain_w0=1; %  0 db afou kathe  monada sthn sixnotita 
%                   omega0 = 4.084070449666731e+03 exei kerdos 0 db, giati etsi prosarmosame 
%                   tin aposvesi se kathe monada

% rithmish kerdous sta 10 dB enisxish sthn zwnh dieleushs
deka_db=10^(10/20); %=3.162277660168380------>3.17 gia na eimaste sigouroi oti tha exoyme 10 db kerdos
deka_db=3.17;
% rithmish kerdous me telestiko enisxith xwris anastrofi
k_telestikou=deka_db/Total_gain_w0;
T_R1=1000;
T_R2=(k_telestikou-1)*T_R1;

% -------------------------------------------------------------------------
% Sinartiseis Metaforas

% gia thn kanonikopoihsh
M1_omega0_kanonik=M1_w0/omega0;
M2_omega0_kanonik=M2_w0/omega0;
M3_omega0_kanonik=M3_w0/omega0;
M4_omega0_kanonik=M4_w0/omega0;
 


T12=series(T1,T2);
T34=series(T3,T4);
T_all=series(T12,T34); % xwris rithmish kerdous
T_bp=k_telestikou*T_all;
invT_all=inv(T_all);
invT_bp=inv(T_bp);
T1_kanonik=tf([0 -(2*M1_Q^2/M1_k_w0)*(M1_omega0_kanonik/M1_Q) 0],[1 (M1_omega0_kanonik/M1_Q) M1_omega0_kanonik^2]);
T2_kanonik=tf([0 -(2*M2_Q^2/M2_k_w0)*(M2_omega0_kanonik/M2_Q) 0],[1 (M2_omega0_kanonik/M2_Q) M2_omega0_kanonik^2]);
T3_kanonik=tf([0 -(2*M3_Q^2/M3_k_w0)*(M3_omega0_kanonik/M3_Q) 0],[1 (M3_omega0_kanonik/M3_Q) M3_omega0_kanonik^2]);
T4_kanonik=tf([0 -(2*M4_Q^2/M4_k_w0)*(M4_omega0_kanonik/M4_Q) 0],[1 (M4_omega0_kanonik/M4_Q) M4_omega0_kanonik^2]);
T12_kanonik=series(T1_kanonik,T2_kanonik);
T34_kanonik=series(T3_kanonik,T4_kanonik);
T_kanonik=series(T12_kanonik,T34_kanonik);


%{
T1 % gia na parw tis sinartiseis metaforas apo to command.
T2
T3
T4
T_bp
%}

plot_transfer_function(T1, [f0 f1 f2 f3 f4]);
plot_transfer_function(T2, [f0 f1 f2 f3 f4]);
plot_transfer_function(T3, [f0 f1 f2 f3 f4]);
plot_transfer_function(T4, [f0 f1 f2 f3 f4]);
plot_transfer_function(T_all, [f0 f1 f2 f3 f4]);
plot_transfer_function(T_bp, [f0 f1 f2 f3 f4]);
ltiview({'bode'}, T1,T2,T3,T4,T_bp);
plot_transfer_function(invT_all, [f0 f1 f2 f3 f4]);
plot_transfer_function(invT_bp, [f0 f1 f2 f3 f4]);
ltiview({'bode'}, T_kanonik);




%% Fourier analysis
% a4=4


in_freq1=(omega0-(omega0-omega1)/2)/(2*pi);
in_freq2=(omega0+(omega0+omega1)/2)/(2*pi);
in_freq3=0.5*omega3/(2*pi);
in_freq4=2.4*omega4/(2*pi);
in_freq5=3.5*omega4/(2*pi);


T=20*(1/100);  % mporw na valw kai prscn me (1/100) pou fainontai kalitera alla den exw megali akriveia.
             
f_s=1000000; % sixnothta deigmatolipsias
dt=1/f_s;
t=0:dt:T-dt;

in=cos((omega0-(omega0-omega1)/2)*t) ...
+0.6*cos((omega0+(omega0+omega1)/2)*t) ...
+cos(0.5*omega3*t)+0.8*cos(2.4*omega4*t)+0.4*cos(3.5*omega4*t);

figure;
plot(t,in);
title('Input signal Matlab with multiple frequencies');
xlabel('Time (s)');
ylabel('Voltage (V)');
axis([0 0.05 -4 4]);


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
axis([0 10000 0 2]);

% output signal 
out=lsim(T_bp,in,t);
figure;
plot(t,out);
title('Output signal of multiple frequencies input Matlab');
xlabel('Time (s)');
ylabel('Voltage (V)');
axis([0 0.05 -4 4]);

% in-out signal plots
figure;
plot(t,in);
hold on;
plot(t,out);
hold off;
title('Input and output of multiple frequencies input Matlab');
xlabel('Time (s)');
ylabel('Voltage (V)');
axis([0 0.05 -4 4]);

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
axis([0 10000 0 4]);






