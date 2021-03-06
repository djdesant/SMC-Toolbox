%% State feedback regulation simulation
% Demonstration of using the function "contl"

%% clear the workspace
clear all;

%% Load the system data
% mat_name = 'aircraft';
%mat_name = 'compped';
mat_name = 'dcmotor';
% mat_name = 'furnobs'; % Chap. 7.4: EXAMPLE: A Temperature Control Scheme
% mat_name = 'hecka';
% mat_name = 'helcopt';
% mat_name = 'huiex';
% mat_name = 'invpen';
% mat_name = 'l1011r';
% mat_name = 'vertint';
load([mat_name,'.mat']);

[n,m]=size(B);
[p,n]=size(C);

%% Design of the compensator

% hyperlane design with desired poles for system dynamics
wn=2.0; %  [rad/s] cut-off angular frequency
zeta=0.95; % [-] damping ratio less 1 --> complex poles
sp=roots([1 2*zeta*wn wn^2]); % second order for first 2 poles

% additional poles (real) by system order greater than 3
if n-m>2
    % to avoid multiple poles iterate through and increase each value
    for i=1:(n-m-2)
        sp=[sp; -i];  % [-1+i1; -1-i1; -1; -2; -3 ...]
    end
end

% Calculate the switching surface using robust eigenstructure assignment
%S=rpp(A,B,sp);

% Design diagonal weighting matrix for the state vector
Q=diag([1 1 10]);

% design of sliding mode poles via LQR
[S,E]=lqcf(A,B,Q);

[G,F]=wzobs(A,B,C,E',-20*ones(1,3));

% Calculate L and Lam required for the linear component of the control law
% Phi=-5*eye(size(S,1));
% [F,P,Lam]=contl(A,B,S,Phi);

%[Hhat,Dhat,S,F,P]=comrobs(A,B,C);

gamma0=-0.5;

rho=1.0; % user defined parameter

SimStopTime=5;

% check roots of closed loop
% The condition number of the eigenvectors is an indicator of the
% robustness with respect to unmatched uncertainty (smaller => better)
M=inv(S*B)*S;
[V,D] = eig(A-B*M);
cond(V)

%% Set initial conditions
x0=zeros(1,n);
x0(1)=1; % set first state initial value
x0(2)=1;
x0(3)=1;

%% Simulate the model
mdl_name='output_mdl';
if ~bdIsLoaded(mdl_name)
    open([mdl_name,'.slx']);
end
set_param(mdl_name,'StopTime',sprintf('%d',SimStopTime));
sim(mdl_name);

%% Plot the results
subplot(2,2,1)
plot(t.Data,u.Data);
grid on;
legend(get_legend('u'))

subplot(2,2,2)
plot(t.Data,s.Data);
grid on;
legend(get_legend('s'));

subplot(2,2,3)
plot(t.Data,y.Data);
grid on;
legend(get_legend('y'));

sgtitle([mdl_name,' - ',mat_name],'Interpreter','None');
