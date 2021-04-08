close all;
clear all;

%% AEM
% AEM=9084

a1=9;
a2=0;
a3=8;
a4=4;

%% Prodiagrafes filtrou

f0=1800;                    %Hz
f1=1200+25*(9-a4);
f2=f0^2/f1;
D=(1/1.8)*((f0^2-f1^2)/f1);
f3=(-D+sqrt(D^2+4*f0^2))/2;
f4=f0^2/f3;

amin=30-a3; %dB
amax=0.5+a4/18; %dB

omega0=2*pi*f0;  %rad/s
omega1=2*pi*f1;
omega2=2*pi*f2;
omega3=2*pi*f3;
omega4=2*pi*f4;

%% Ypologismoi 

bw=omega2-omega1;

Omega_p=1;
Omega_s=(omega2-omega1)/(omega4-omega3);

n_arxiko=acosh(sqrt((10^(amin/10)-1)/(10^(amax/10)-1)))/acosh(Omega_s);
%tipos 9.83 afou h taxi chebysev idia me IC. Me swsth kanonikopoihsh.  

n=4;  %afou n_arxiko=3.4180

epsilon=1/sqrt(10^(amin/10)-1);
alpha=(1/n)*asinh(1/epsilon);

Omega_hp=1/cosh((1/n)*acosh(1/epsilon)); %sixnotita hmisias isxios

%---------------------------------------------------------------------------
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

%--------------------------------------------------------------------------

% (wo,Q) gia Chebysev
w0_12=sqrt(sigma12^2+pole_omega12^2);
Q_12=1/(2*cos(atan(pole_omega12/sigma12)));

w0_34=sqrt(sigma34^2+pole_omega34^2);
Q_34=1/(2*cos(atan(pole_omega34/sigma34)));

% antristrofi twn pollwn gia na vrw polous tou IC
% kai klimakopoihsh
% ta Q paramenoun idia
w0_12_ic=Omega_s/w0_12;
w0_34_ic=Omega_s/w0_34;

% midenika
Wz1=sec(pi/8);  %k=1
Wz2=sec(3*pi/8); %k=3

% klimakopoihsh mhdenikwn
Wz1_ic=Omega_s*Wz1;
Wz2_ic=Omega_s*Wz2;

% antistrefoume tous polous tou ic
W_12_bar=1/w0_12_ic;     %ta Q menoun idia
W_34_bar=1/w0_34_ic;

% antistrefoume ta midenika tou ic
Wz1_ic_bar=1/Wz1_ic;
Wz2_ic_bar=1/Wz2_ic;

% poloi anodiavaths sinarthshs
hp_S12=-W_12_bar/(2*Q_12);
hp_W12=sqrt(W_12_bar^2-hp_S12^2); 
hp_S34=-W_34_bar/(2*Q_34);
hp_W34=sqrt(W_34_bar^2-hp_S34^2);

qc=omega0/bw;

% metasxhmatismos 1ou Migadikou Polou(1)
P1_S12=abs(hp_S12);
P1_W12=abs(hp_W12);

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
P2_S34=abs(hp_S34);
P2_W34=abs(hp_W34);

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

% metasxhmatismos 1ou fantastikou midenikou
Z1_K=2+(Wz1_ic_bar^2 / qc^2);
Z1_x=(Z1_K + sqrt(Z1_K^2 - 4))/2;
% 1o fantasiko mideniko
Z1_wz1=omega0*sqrt(Z1_x);
% 2o fantastiko mideniko    + 2 poloi sto 0
Z1_wz2=omega0/(sqrt(Z1_x));

% metasxhmatismos 2ou fantastikou midenikou
Z2_K=2+(Wz2_ic_bar^2 / qc^2);
Z2_x=(Z2_K + sqrt(Z2_K^2 - 4))/2;
% 1o fantasiko mideniko
Z2_wz1=omega0*sqrt(Z2_x);
% 2o fantastiko mideniko    + 2 poloi sto 0
Z2_wz2=omega0/(sqrt(Z2_x));

%% Ylopoihsh ths synarthshs metaforas
% ylopoioyme 4 vathmides----> 2 LPN kai 2 HPN

% 1h Monada LPN boctor
M1_w0=P2_w01;
M1_wz=Z1_wz2;
% Thetw-klimakopoiw
M1_W0=1;
M1_Wz=M1_wz/M1_w0;
% sinexizw
M1_Q=P2_Q;
% parametros k1
M1_W0_Wz=(M1_W0/M1_Wz)^2; %=0.897306096769159
% epilegw k1=0.95
M1_k1=0.95;

M1_R1=2/(M1_k1*M1_Wz^2 - 1);
M1_R2=1/(1- M1_k1);
M1_R3=0.5*(M1_k1/M1_Q^2+M1_k1*M1_Wz^2-1);
M1_R4=1/M1_k1;
M1_R5=1;
M1_R6=M1_R5;
M1_C1=M1_k1/(2*M1_Q);
M1_C2=2*M1_Q;
% kerdos k
M1_k=2/((M1_k1/M1_Q^2)+(M1_k1*M1_Wz^2)+1);

% klimakopoiw
M1_kf=M1_w0;
M1_km=M1_C1/(M1_kf*10^(-7));

M1_R1_k=M1_km*M1_R1;
M1_R2_k=M1_km*M1_R2;
M1_R3_k=M1_km*M1_R3;
M1_R4_k=M1_km*M1_R4;
M1_R5_k=M1_km*M1_R5;
M1_R6_k=M1_km*M1_R6;
M1_C1_k=(M1_C1/(M1_kf*M1_km));
M1_C2_k=(M1_C2/(M1_kf*M1_km));

% sinartisi metaforas ths monadas
T1=tf([M1_k 0 M1_k*M1_wz^2],[1 M1_w0/M1_Q M1_w0^2]);

% -----------------------------------------------------------------

% 2h Monada LPN boctor
M2_w0=P1_w01;
M2_wz=Z2_wz2;
% Thetw-klimakopoiw
M2_W0=1;
M2_Wz=M2_wz/M2_w0;
% sinexizw
M2_Q=P1_Q;
% parametros k1
M2_W0_Wz=(M2_W0/M2_Wz)^2; %=0.954030281555381
% epilegw k1=0.98
M2_k1=0.98;

M2_R1=2/(M2_k1*M2_Wz^2 - 1);
M2_R2=1/(1- M2_k1);
M2_R3=0.5*((M2_k1/M2_Q^2)+(M2_k1*M2_Wz^2)-1);
M2_R4=1/M2_k1;
M2_R5=1;
M2_R6=M2_R5;
M2_C1=M2_k1/(2*M2_Q);
M2_C2=2*M2_Q;
% kerdos k
M2_k=2/((M2_k1/M2_Q^2)+(M2_k1*M2_Wz^2)+1);

% klimakopoiw
M2_kf=M2_w0;
M2_km=M2_C1/(M2_kf*10^(-7));

M2_R1_k=M2_km*M2_R1;
M2_R2_k=M2_km*M2_R2;
M2_R3_k=M2_km*M2_R3;
M2_R4_k=M2_km*M2_R4;
M2_R5_k=M2_km*M2_R5;
M2_R6_k=M2_km*M2_R6;
M2_C1_k=(M2_C1/(M2_kf*M2_km));
M2_C2_k=(M2_C2/(M2_kf*M2_km));

% sinartisi metaforas ths monadas
T2=tf([M2_k 0 M2_k*M2_wz^2],[1 M2_w0/M2_Q M2_w0^2]);
% ------------------------------------------------------------------------

% 3h monada HPN boctor

M3_w0=P1_w02;
M3_wz=Z2_wz1;
M3_W0=1;
M3_Wz=M3_wz/M3_w0;
M3_Q=P1_Q;
M3_C2= 10^(-7);

M3_Wz_W0=1/(1-(M3_Wz^2 / M3_W0^2)); %=21.753450615642212
% ara isxiei M3_Q < 21.753450615642212

BoctorHighPass(M3_wz, M3_w0, M3_Q, 2000, 10^-7);

M3_R1=circuit.R_1;
M3_R2=circuit.R_2;
M3_R3=circuit.R_3;
M3_R4=circuit.R_4;
M3_R5=circuit.R_5;
M3_R6=circuit.R_6;
M3_C1=circuit.C_1;
M3_H=circuit.H;

T3=tf([M3_H 0 M3_H*M3_wz^2],[1 M3_w0/M3_Q M3_w0^2]);

% ------------------------------------------------------------------

% 4h monada HPN boctor

M4_w0=P2_w02;
M4_wz=Z1_wz1;
M4_W0=1;
M4_Wz=M4_wz/M4_w0;
M4_Q=P2_Q;
M4_C2= 10^(-7);

M4_Wz_W0=1/(1-(M4_Wz^2 / M4_W0^2)); %=9.737676420304558
% ara isxiei M4_Q < 9.737676420304558

BoctorHighPass(M4_wz, M4_w0, M4_Q, 1000, 10^-7);

M4_R1=circuit.R_1;
M4_R2=circuit.R_2;
M4_R3=circuit.R_3;
M4_R4=circuit.R_4;
M4_R5=circuit.R_5;
M4_R6=circuit.R_6;
M4_C1=circuit.C_1;
M4_H=circuit.H;

T4=tf([M4_H 0 M4_H*M4_wz^2],[1 M4_w0/M4_Q M4_w0^2]);

% --------------------------------------------------------------------

% sinoliko kerdos ths sinolikhs sinarthshs metaforas
Total_gain=M1_k*M2_k*M3_H*M4_H; %=3.665967518293267
% rithmish kerdous- thelw aposvesi sta 5db
pente_db=1.78; %1.778279410038923---> 1.78 gia na eimaste sigouroi
k_telestikou1=pente_db/Total_gain; % telestikos se anastrefon sindesmologia
T1_R1=1000;
T1_R2=T1_R1*k_telestikou1;
k_telestikou2=1;  % gia na exw to arxiko shma kai oxi to anastrofo
T2_R1=1000;       % vazw allon enan telestiko se anastrefon sindesmologia
T2_R2=1000;       % me kerdos 1 (0db)
%---------------------------------------------------------------------------
% Sinartiseis metaforas

% gia thn kanonikopoihmenh
M1_omega0=M1_w0/omega0;
M2_omega0=M2_w0/omega0;
M3_omega0=M3_w0/omega0;
M4_omega0=M4_w0/omega0;
M1_omegaZero=M1_wz/omega0;
M2_omegaZero=M2_wz/omega0;
M3_omegaZero=M3_wz/omega0;
M4_omegaZero=M4_wz/omega0;



T12=series(T1,T2);
T34=series(T3,T4);
T_all=series(T12,T34); % xwris rithmish kerdousS
T_be=k_telestikou1*k_telestikou2*T_all;
invT_all=inv(T_all);
invT_be=inv(T_be);
T1_kanonik=tf([M1_k 0 M1_k*M1_omegaZero^2],[1 M1_omega0/M1_Q M1_omega0^2]);
T2_kanonik=tf([M2_k 0 M2_k*M2_omegaZero^2],[1 M2_omega0/M2_Q M2_omega0^2]);
T3_kanonik=tf([M3_H 0 M3_H*M3_omegaZero^2],[1 M3_omega0/M3_Q M3_omega0^2]);
T4_kanonik=tf([M4_H 0 M4_H*M4_omegaZero^2],[1 M4_omega0/M4_Q M4_omega0^2]);
T12_kanonik=series(T1_kanonik,T2_kanonik);
T34_kanonik=series(T3_kanonik,T4_kanonik);
T_kanonik=series(T12_kanonik,T34_kanonik);


plot_transfer_function(T1, [f0 f1 f2 f3 f4]);
plot_transfer_function(T2, [f0 f1 f2 f3 f4]);
plot_transfer_function(T3, [f0 f1 f2 f3 f4]);
plot_transfer_function(T4, [f0 f1 f2 f3 f4]);
plot_transfer_function(T_all, [f0 f1 f2 f3 f4]);
plot_transfer_function(T_be, [100 f0 f1 f2 f3 f4]);
ltiview({'bode'}, T1,T2,T3,T4,T_be)
plot_transfer_function(invT_all*Total_gain, [f0 f1 f2 f3 f4]);
plot_transfer_function(invT_be, [100 f0 f1 f2 f3 f4]);
ltiview({'bode'}, T_kanonik/Total_gain);



%% Fourier analysis

% a4=4


in_freq1=(omega0-(omega0-omega3)/2)/(2*pi);
in_freq2=(omega0+(omega0+omega3)/3)/(2*pi);
in_freq3=0.4*omega1/(2*pi);
in_freq4=2.5*omega2/(2*pi);
in_freq5=3*omega2/(2*pi);


T=20*(1/100);  % mporw na valw kai prscn me (1/100) pou fainontai kalitera alla den exw megali akriveia.
             
f_s=1000000; % sixnothta deigmatolipsias
dt=1/f_s;
t=0:dt:T-dt;

in=0.5*cos((omega0-(omega0-omega3)/2)*t) ...
+0.8*cos((omega0+(omega0+omega3)/3)*t) ...
+0.8*cos(0.4*omega1*t)+0.6*cos(2.5*omega2*t)+1.2*cos(3*omega2*t);

figure;
plot(t,in);
title('Input signal Matlab with multiple frequencies');
xlabel('Time (s)');
ylabel('Voltage (V)');
axis([0 0.03 -10 10]);


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
out=lsim(T_be,in,t);
figure;
plot(t,out);
title('Output signal of multiple frequencies input Matlab');
xlabel('Time (s)');
ylabel('Voltage (V)');
axis([0 0.03 -10 10]);

% in-out signal plots
figure;
plot(t,in);
hold on;
plot(t,out);
hold off;
title('Input and output of multiple frequencies input Matlab');
xlabel('Time (s)');
ylabel('Voltage (V)');
axis([0 0.03 -10 10]);

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
