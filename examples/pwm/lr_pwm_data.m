clc; clear all;

%% Plant Parameter
R_l=4; % [Ohm]
L_l=0.150; % [H]
%T_s=1e-2; % [s]
T_s=1e-2; % [s]
w_s = 1/T_s; % [Hz]
f_s = w_s/(2*pi);
%f_s=20; % [Hz] filter cut-off frequency
R_s=0.1;
C_s=T_s/R_s;
w_gR=2*pi*10; % [rad/s]
f_pwm=125; % [Hz]
T_pwm=1/f_pwm; % [s]
%T_pwmsim=T_pwm/10; % [s]
T_sim=T_pwm/200; % [s] min. less then 1/100 of pwm period
%w_s = 2*pi*1/T_s;
U_e=12.5; % [V]
d=0.1;

%% State Space Matrices
% Define default Storage-Data-Types
%
ntemp = 16;
ShortType    = sfix(1*ntemp);
BaseType     = sfix(2*ntemp);
LongType     = sfix(4*ntemp);
%ExtendedType = sfix(8*ntemp);

%
% Define default Sample Time
%
SampleTime = 0.001;   % sec

Ac_s=[-1/(C_s*R_s)    1/C_s;    % U_c(t)
           -1/L_l       -R_l/L_l];   % i_L(t)
Bc_s=[0;1/L_l];
Cc_s=[1/R_s 0];
Ob = obsv(Ac_s,Cc_s);
% Number of unobservable states
unob = length(Ac_s)-rank(Ob);
%Cc_s=[0 1];
[n,m]=size(Bc_s);
[p,n]=size(Cc_s);
Dc_s=zeros(p,m);

% convert to discrete time
sys_s = ss(Ac_s,Bc_s,Cc_s,Dc_s);
sys_z = c2d( sys_s, SampleTime, 'tustin');

% convert to canonical state space form
%[Ac_s, Bc_s, Cc_s, Dc_s ] = ssdata(sys_s);
[Ac_z, Bc_z, Cc_z, Dc_z ] = ssdata(sys_z);

% convert to a balanced realization
sysb_z = balreal(sys_z);

[Ab_z, Bb_z, Cb_z, Db_z ] = ssdata(sysb_z);

% specify initial input
u0 = zeros(1,m);

% derive initial conditions for balanced realization
xb0_z  = ( eye(length(Ab_z)) - Ab_z ) \ Bb_z * u0;

% derive initial conditions for canonical realization
xc0_z  = ( eye(length(Ac_z)) - Ac_z ) \ Bc_z * u0;
xc0_s  = ( -Ac_s ) * Bc_s .\ u0;

%%%%%%%%%%%%%%%%%%%%%%%%%% Observer Design %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ao_s=Ac_s;
Bo_s=Bc_s;
Ao_s=[-1/(C_s*R_s)    1/C_s;    % U_c(t)
           -1/L_l       -1.0*R_l/L_l];   % i_L(t)
Bo_s=[0;1/(L_l)];
To_inv=eye(n);
%Co_s=Cc_s(2,:);
%Co_s=[1/R_s 0];
%Co_s=[1/R_s 0; 0 1];
Co_s=[0 1]; % for getting rank(CB)=m (Markov parameter has full rank)
[p,n]=size(Co_s);
%[At, Bv, Tr, CA]=SMCCanForm(Ao_s, Bo_s);
% Ct=Co_s*inv(Tr);
% calculate eigenvalues by minimizing quadratic cost function
Q=diag([1 1]);
Q=eye(n);
[~,E]=lqcfCA(Ao_s,Bo_s,Q);

% % Design the dynamic compensator
% ped = -10; %  desired pole(s) for error dynamics (Lo)
% ped=[];
% psm = E'; % reduced order sliding motion poles
% prsd = -10; % desired pole(s) for range space dynamics
% 
% % design compensator with reduced observer
% [Hhat,Dhat,S,L,P,Lam,T]=comrobs(Ao_s,Bo_s,Co_s,ped,psm,prsd);

% move the observer sliding mode poles further to the left
% multiplied by factor
psm=[];
po_factor=5;
Eid=1;
if ~isempty(E)
    for i=1:(n-p) % ToDo check nn-pp
        % check if the poles are conjugated complex
        if ~isreal(E(Eid))
            psm(i)=po_factor*real(E(Eid));
        else
            psm(i)=po_factor*E(Eid);
        end

        if length(E)>(n-p)
            Eid=Eid+2;
        else
            Eid=Eid+1;
        end
    end
end

% check the length of psm
[~,~,~,~,r]=outfor(Ao_s,Bo_s,Co_s);
if (n-p-r)>0
    assert(length(psm)==(n-p-r));
else
    assert(length(psm)==(n-p)); % ToDo check
end

% set the estimation error poles
per=-1*ones(1,p);
assert(length(per)==p);

[G,F]=wzobs(Ao_s,Bo_s,Co_s,psm,per);

%%%%%%%%%%%%%%%%% Design of PI State Space Controller %%%%%%%%%%%%%%%%%%%%%
[Ai,Bi]=intac(Ao_s,Bo_s,Co_s);

% weighting matrix for the PI controller
Qi = eye(length(Ai));
Qi=diag([50 5 1]);

[S,E]=lqcfCA(Ai,Bi,Qi);

% Eopt=E;
% Eopt(Eopt>-50)=-50;
% S=rpp(Ai,Bi,Eopt);

% design the dynamics of the PI controller
Phi=-150*eye(size(S,1));

[L,Lr,Lrdot,Sr,Lam,P]=contliaCA(Ao_s,Bo_s,Co_s,S,Phi);
Lam_inv=inv(Lam);

rho=1;
delta=0.001;

FR=-diag([1/0.01]);

z0=zeros(1,n);
z0(2)=1;
rhoo=5;
deltao=0.001;

%% Deviated Parameters
% big load
% U_e = 10; % [V]
% R_l = 5.6;
% L_l = 300e-3;

% small load
% U_e = 15; % [V]
% R_l = 2.8;
% L_l = 150e-3;