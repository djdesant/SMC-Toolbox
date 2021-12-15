clc
clear all
close all hidden

%% Parameter Tiefsetzsteller
R=4; % [Ohm]
L=0.1; % [H]
R_s=0.1;
tau_s=1e-2; % [s] Rising time of RC
C_s=tau_s/R_s;

U_E=12; % [V]

T_s=1/150; %1.5e-3;

%% Parameter SMC
c2=1;
wg=1/(C_s*R_s);
c1=wg;

u_stern = 10;
i_stern = u_stern/R;

u_vektor = [i_stern; u_stern];

%Matrizen SPAM
A=[-1/(R_s*C_s) 1/C_s; -1/L -R/L];

B=[0; 1/L];

B_z = [0; R/L];

C_t = [c1 c2];

T_R = [C_s*L 0; -L/R L];

u_vektor_R = T_R*u_vektor;


return

sim('SMC_sim', 200e-3);
figure(1)
x1R = xR.signals(1).values - u_vektor_R(1);
x2R = xR.signals(2).values - u_vektor_R(2);
% 
plot(x1R, x2R, 'b', 'LineWidth',2);
hold on;
x1R_sw = -2e-3:1e-5:0;
x2R_sw = -c1/c2*x1R_sw;

plot(x1R_sw, x2R_sw, 'r--', 'LineWidth',2);
grid on;
legend('transf. Zustandsgrößen', 'Schaltkennlinie (transf)');
xlabel('x_{1R}');
ylabel('x_{2R}');



% return

figure(2)
% x1 = ScopeData.signals(2).values - u_vektor(1);
% x2 = ScopeData.signals(3).values - u_vektor(2);
x1 = ScopeData.signals(2).values;
x2 = ScopeData.signals(3).values;
plot(x1, x2, 'b', 'LineWidth',2);
hold on;

invTR = inv(T_R);
vec = invTR*[x1R_sw; x2R_sw];
% x1_sw = vec(1,:);
% x2_sw = vec(2,:);
% x1_sw = 0:1e-3:0.25;
% x2_sw = -c1/c2*x1_sw;
% 
% plot(x1_sw, x2_sw, 'r--', 'LineWidth',2);
% axis([0 1.2e-3 -0.05 0.05]);
grid on;
legend('Zustandsgrößen', 'Schaltkennlinie');

xlabel('i_{L}');
ylabel('u_{C}');