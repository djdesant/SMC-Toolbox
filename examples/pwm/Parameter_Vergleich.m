clc; clear all;

%% Streckenparameter
R=50;
L=200e-3;
C=500e-6;
w_gR=2*pi*35;
T_s=1.5e-3;
w_s = 2*pi*1/T_s;
U_e=12;
d=0.5;

%% SPAM
A_strecke=[0 -1/L;
    1/C -1/(R*C)];
B_strecke=[1/L;
    0];
B_z_strecke=[0;
    1/C];
C_strecke=[0 1];
D_strecke=0;

%% Sanders/Verghese
A_sv = [0 0 0 -1/L 0 0;...
        0 0 w_s 0 -1/L 0;...
        0 -w_s 0 0 0 -1/L;...
        1/C 0 0 -1/R/C 0 0;...
        0 1/C 0 0 -1/R/C w_s;...
        0 0 1/C 0 -w_s -1/R/C];

B_sv = [2/L*d;...
        0;...
        -2/pi/L*sin(pi*d);...
        0;...
        0;...
        0];
    
C_sv = [1/2 1 1 0 0 0;...
        0 0 0 1/2 1 1];

%% PI-Zustandsregler
A_R=[0 1;
    -1/(L*C) -1/(R*C)];
B_R=[0;
    1];
C_R=[1/(L*C) 0];
T_R=[0 L*C;
    L -L/R];
r_1=3*w_gR^2-1/(L*C)-w_gR^2;
r_2=3*w_gR-1/(R*C);
R_ZR=[r_1 r_2];
Q=L*C;
r_I=w_gR^3;
r_P=w_gR^2;