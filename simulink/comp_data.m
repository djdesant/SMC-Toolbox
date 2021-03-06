%% Compensator simulation
% Demonstration of using the function "comrobs"

%% clear the workspace
clear all;

%% Load the system data
mat_name = 'compped';  % Chapter 5.6.2: Design Example 1
mat_name = 'invpen';  % Chap. 5.6.3 Design Example 2: Inverted Pendulum
mat_name = 'l1011r';  % Chap. 5.7.1 Aircraft Example
mat_name = 'dcmotor'; % Chap. 3.6.4 Example: Control of a DC Motor
%mat_name = 'bwb'; % Blended Wing Body

% mat_name = 'vertint'; % Chap 4.5 Design Study: Pitch-Pointing Flight Controller
% mat_name = 'furnobs'; % Chap. 7.4: EXAMPLE: A Temperature Control Scheme
% mat_name = 'hecka';
% mat_name = 'helcopt';
% mat_name = 'huiex';

load([mat_name,'.mat']);

if strcmp(mat_name,'bwb')
    A=BWB.lat.A;
    B=BWB.lat.B;
    C=eye(length(A));
    D=zeros(size(C,1),size(B,2));
end

[n,m]=size(B); % n-state vector length, m-number of inputs
[p,n]=size(C); % n-state vector length, p-number of outputs

switch mat_name
    case 'compped'
        CA=eye(m);

        % define desired dynamics
        ped=-2.5;       % poles for error dynamics
        psm=[-1 -1.5];  % poles for sliding motion
        prsd=-5*ones(1,m);        % poles for range space dynamics

        % design compensator with reduced observer
        [Hhat,Dhat,S,L,P,Lam,T]=comrobs(A,B,C,ped,psm,prsd);

        x0=zeros(1,n); % [xr, x11, x12]
        x0(2)=1;
        xc0=zeros(1,n-p);

        rho=1;
        delta=0.001;

        SimStopTime=10;

    case 'invpen'
        CA=eye(m);

        % design of sliding mode poles via LQR
        Q=diag([10 1 1 0.1]);
        [~,E]=lqcf(A,B,Q);
        
        % Design the dynamic compensator
        ped = -10; %  desired pole(s) for error dynamics (Lo)
        psm = E; % reduced order sliding motion poles
        prsd = -6*ones(1,m); % desired pole(s) for range space dynamics

        % design compensator with reduced observer
        [Hhat,Dhat,S,L,P,Lam,T]=comrobs(A,B,C,ped,psm,prsd);

        x0=zeros(1,n);
        x0(2)=0.1;
        xc0=zeros(1,n-p);
        xc0(1)=-0.9855;

        rho=1;
        delta=0.001;

        SimStopTime=5;

    case 'l1011r'
        CA=eye(m);

        % Design diagonal weighting matrix for the state vector
        Q=diag([5 1 1 5 5]);

        % design of sliding mode poles via LQR
        [~,E]=lqcf(A,B,Q);

        % Design the dynamic compensator
        ped = -5; %  desired pole(s) for error dynamics (Lo)
        psm = E'; % reduced order sliding motion poles
        prsd = -5*ones(1,m); % desired pole(s) for range space dynamics

        % design compensator with reduced observer
        [Hhat,Dhat,S,L,P,Lam,T]=comrobs(A,B,C,ped,psm,prsd);

        x0=zeros(1,n);
        x0(1)=1;
        xc0=zeros(1,n-p);
        xc0(1)=0;

        rho=1;
        delta=0.001;

        SimStopTime=10;

    case 'dcmotor'
        C=C([1,3],:); % partial output matrix to allow a reduced observer design
        [p,n]=size(C);

        % canonical form
        [At, Bv, Tr, CA]=SMCCanForm(A, B);
        Ct=C*inv(Tr);
        [nn,ll]=size(Bv);
        [pp,nn]=size(Ct);

        % Design diagonal weighting matrix for the state vector
        Q=diag([5 5 1]);

        % design of sliding mode poles via LQR
        [~,E]=lqcfCA(At,Bv,Q);

        % Design the dynamic compensator
        ped = -10; %  desired pole(s) for error dynamics (Lo)
        psm = E'; % reduced order sliding motion poles
        prsd = -20*ones(1,m); % desired pole(s) for range space dynamics

        % design compensator with reduced observer
        [Hhat,Dhat,S,L,P,Lam,T]=comrobs(At,Bv,Ct,ped,psm,prsd);

        x0=zeros(1,n);
        x0(1)=0.1;
        xc0=zeros(1,n-p);
        xc0(1)=0;

        rho=1;
        delta=0.001;

        SimStopTime=5;

    case 'bwb'
        C=C([1,2,3],:);
        [p,n]=size(C); % n-state vector length, p-number of outputs
        D=zeros(p,m);
        
        % canonical form of the overatuated system
        [At, Bv, Tr, CA]=SMCCanForm(A, B);
        Ct=C*inv(Tr);
        [nn,ll]=size(Bv);
        [pp,nn]=size(Ct);

        % Design diagonal weighting matrix for the state vector
        Qdiag=linspace(1,nn,nn);
        %Qdiag(1:2)=5;
        Q=diag(Qdiag);
        
        % design of sliding mode poles via LQR
        [S,E]=lqcfCA(At,Bv,Q);

        % Design the dynamic compensator
        ped = -10; %  desired pole(s) for error dynamics (Lo)
        %ped=[];
        psm = E'; % reduced order sliding motion poles
        prsd = -10*ones(1,ll); % desired pole(s) for range space dynamics

        % design compensator with reduced observer
        [Hhat,Dhat,S,L,P,Lam,T]=comrobs(At,Bv,Ct,ped,psm,prsd);

        x0=zeros(1,n);
        x0(1)=0.1;
        xc0=zeros(1,n-p);
        xc0(1)=0;

        rho=1;
        delta=0.001;

        SimStopTime=20;
end

Lo=ped;

%% simulate the model
mdl_name='comp_mdl';
if ~bdIsLoaded(mdl_name)
    open([mdl_name,'.slx']);
end
set_param(mdl_name,'StopTime',sprintf('%d',SimStopTime));
sim(mdl_name);

%% plot simulation results
subplot(2,2,1)
plot(t.Data,s.Data);
grid on;
legend(get_legend('s'));

xc=xhat.Data(:,2);
x11=x.Data(:,2);
x12=x.Data(:,3);
ec=xc-x11-Lo*x12;

xr=x.Data(:,1);
zr=xhat.Data(:,1);
er=zr-xr;

subplot(2,2,2)
plot(t.Data,ec,t.Data,er);
grid on;
legend([get_legend('ec');get_legend('er')]);

subplot(2,2,3)
plot(t.Data,x.Data);
grid on;
legend(get_legend('x'));

subplot(2,2,4)
plot(t.Data,xhat.Data);
grid on;
legend(get_legend('xhat'));

sgtitle([mdl_name,' - ',mat_name],'Interpreter','None');
