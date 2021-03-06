%% Example originally created by Liam Vile
% Demonstration of using the function "contl" and "contlia"
% with special canonical representation

%% clear the workspace
clear all;

%% Load the system data

mat_name = 'compped';  % Chapter 5.6.2: Design Example 1
mat_name = 'invpen';  % Chap. 5.6.3 Design Example 2: Inverted Pendulum
mat_name = 'l1011r';  % Chap. 5.7.1 Aircraft Example
mat_name = 'dcmotor'; % Chap. 3.6.4 Example: Control of a DC Motor
mat_name = 'bwb'; % Belnded Wing Body

% mat_name = 'vertint'; % Chap 4.5 Design Study: Pitch-Pointing Flight Controller
% mat_name = 'furnobs'; % Chap. 7.4: EXAMPLE: A Temperature Control Scheme
% mat_name = 'hecka';
% mat_name = 'helcopt';
% mat_name = 'huiex';

load([mat_name,'.mat']);

% [n,m]=size(B); % n-state vector length, m-number of inputs
% [p,n]=size(C); % n-state vector length, p-number of outputs

switch mat_name
    case 'compped'
        % define desired dynamics
        ped=-2.5;       % poles for error dynamics
        psm=[-1 -1.5];  % poles for sliding motion
        prsd=-5*ones(1,m);        % poles for range space dynamics

        x0=zeros(1,n); % [xr, x11, x12]
        x0(2)=1;
        xc0=zeros(1,n-p);

        rho=1;
        delta=0.001;

        SimStopTime=10;

    case 'invpen'
        % design of sliding mode poles via LQR
        Q=diag([10 1 1 0.1]);
        [~,E]=lqcf(A,B,Q);
        
        % Design the dynamic compensator
        ped = -10; %  desired pole(s) for error dynamics (Lo)
        psm = E; % reduced order sliding motion poles
        prsd = -6*ones(1,m); % desired pole(s) for range space dynamics

        x0=zeros(1,n);
        x0(2)=0.1;
        xc0=zeros(1,n-p);
        xc0(1)=-0.9855;

        rho=1;
        delta=0.001;

        SimStopTime=5;

    case 'l1011r'
        % Design diagonal weighting matrix for the state vector
        Q=diag([5 1 1 5 5]);

        Qi=eye(ll);
        FR=-eye(ll);

        x0=zeros(1,n);
        x0(1)=1;
        xc0=zeros(1,n-p);
        xc0(1)=0;

        rho=1;
        delta=0.001;

        rhoo=5;
        deltao=0.001;

        SimStopTime=10;

    case 'dcmotor'
        Cc=C([3],:); % partial output matrix to allow a reduced observer design
        [n,m]=size(B); % n-state vector length, m-number of inputs
        [p,n]=size(C);
        D=zeros(p,m);
        ADPT=0; PC=1; PK=1; C0=0; K0=0;

        [At, Bv, Tr, CA]=SMCCanForm(A, B);
        [nn,ll]=size(Bv);

        % Design diagonal weighting matrix for the state vector
        Q=diag([5 5 1]);

        [Ai,Bvi]=intac(At,Bv,Cc*inv(Tr));

        ii=size(Bvi,1)-nn;
        Qi=eye(ii);
        FR=-eye(ii);

%         [L, P, Lam_inv, S, E, CA,...
%             Ao, Bo, Co, To_inv, G, F, Gl, Gn, P2]=SMCDesign_init(A, B, Cc, Q);
%         ll=size(S,1);
%         %%%%%%C=zeros(ll, nn);
%         L=[zeros(ll, ll) L];
%         S=[zeros(ll, ll) S];
%         Sr=zeros(ll, ll);
%         Lr=[zeros(ll, ll)];
%         Lrdot=[zeros(ll, ll)];

        FR=-eye(ll);

        QQ=[Qi zeros(ii, nn)
            zeros(nn, ii) Q];
        [L, Lr, Lrdot, P, Lam_inv, S, Sr, E, CA,...
            Ao, Bo, Co, To_inv, G, F, Gl, Gn, P2]=SMCDesign_initIA(A, B, Cc, QQ);

        x0=zeros(1,n);
        x0(1)=0.1;
        z0=zeros(1,n);
        z0(1)=0;

        rho=1;
        delta=0.001;

        rhoo=5;
        deltao=0.001;

        SimStopTime=5;

    case 'bwb'
        A=BWB.lat.A;
        B=BWB.lat.B;
        [n,m]=size(B); % n-state vector length, m-number of inputs
        C=eye(n);
        [p,n]=size(C); % n-state vector length, p-number of outputs
        D=zeros(p,m);

        ADPT=0; PC=1; PK=1; C0=0; K0=0;

        Cc=[1 0 0 0
            0 1 0 0];
        [p,n]=size(Cc);

        [At, Bv, Tr, CA]=SMCCanForm(A, B);
        [nn,ll]=size(Bv);

        %ii=size(Bvi,1)-nn;
%         Q=eye(n);
%         FR=-eye(n);

%         [L, P, Lam_inv, S, E, CA,...
%             Ao, Bo, Co, To_inv, G, F, Gl, Gn, P2]=SMCDesign_init(A, B, Cc, Q);
%         ll=size(S,1);
%         %%%%%%C=zeros(ll, nn);
%         L=[zeros(ll, ll) L];
%         S=[zeros(ll, ll) S];
%         Sr=zeros(ll, ll);
%         Lr=[zeros(ll, ll)];
%         Lrdot=[zeros(ll, ll)];

        Q=eye(nn);
        Qi=eye(ll);
        QQ=[Qi zeros(ll, nn)
            zeros(nn, ll) Q];
        [L, Lr, Lrdot, P, Lam_inv, S, Sr, E, CA,...
            Ao, Bo, Co, To_inv, G, F, Gl, Gn, P2]=SMCDesign_initIA(A, B, Cc, QQ);

        %Qi=eye(n-p);
        FR=-eye(ll);

        x0=zeros(1,n);
        x0(1)=0.1;

        z0=zeros(1,n);
        z0(1)=0;

        rho=5;
        delta=0.01;

        rhoo=15;
        deltao=0.01;

        SimStopTime=40;
end


%% simulate the model
mdl_name='Example';
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



