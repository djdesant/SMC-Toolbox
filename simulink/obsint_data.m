%% Observer simulation with integral action
% observer based controller design with integral action
% Demonstration of using the function "contlia"

%% clear the workspace
clear all;
clc

%% Load the system data

mat_name = 'furnobs';  % Chap. 7.4: EXAMPLE: A Temperature Control Scheme
%mat_name = 'vertint';  % Chap. 4.5 Design Study: Pitch-Pointing Flight Controller
                       % Chap. 4.5.2 Integral Action Based Design

% mat_name = 'aircraft';
% mat_name = 'compped';
mat_name = 'dcmotor';
% mat_name = 'hecka';
% mat_name = 'helcopt';
% mat_name = 'huiex';
% mat_name = 'invpen';
%mat_name = 'l1011r';

load([mat_name,'.mat']);

[n,m]=size(B);
[p,n]=size(C);
%sys_info(ss(A,B,C,zeros(p,m)),mat_name);

switch mat_name
    case 'vertint'

%         C=C(2,:); % partial output matrix to allow a reduced observer design
%         C=[C; [0,1,0,0,0]];
%         [p,n]=size(C);
% 
%         sys_info(ss(A,B,C,zeros(p,m)),mat_name);

        [Aa,Ba]=intac(A,B,C);

        %
        nocomp=1;
        lambda=[-5.6 4.2 -1 -0.4 -0.7];
        %specent=eigbuil(7,2,nocomp)
        % 4 Eigenvectors (only the first one is complex)
        %        1re 1im   2   3   4
        specent=[NaN NaN NaN NaN NaN;
                 NaN NaN NaN NaN NaN;
                   0   0   1   0   1;
                   1 NaN   0   1   0;
                 NaN   1   0 NaN   0;
                 NaN NaN NaN NaN NaN;
                 NaN NaN NaN NaN NaN];

        [S,V]=dea(Aa,Ba,specent,lambda,nocomp);

        Phi=-20*eye(size(Ba,2));

        S=inv(S*Ba)*S;

        % Calculate the switching surface using robust eigenstructure assignment
        Sob=rpp(A,B,[-5.6+1i*4.2,-5.6-1i*4.2,-1]);
        [F,~,~]=contl(A,B,Sob,Phi);

        G=-inv(C*inv(A+B*F)*B);
        %C=[C; [0,1,0,0,0];[0,0,1,0,0];[0,0,0,1,0]];
        %[G,F]=wzobs(A,B,C,[-5.6+1i*4.2,-5.6-1i*4.2,-1,-0.4,-0.7],-10);

        % Design of observer based controller with integral action
        [L,Lr,Lrdot,Sr,Lam,P]=contlia(A,B,C,S,Phi);

        x0=zeros(1,n); % [xr, x11, x12]
        x0(2)=1;
        z0=zeros(1,n);

        Gamma=diag([-0.9 -0.7]);

        rho=1;
        delta=0.001;

        SimStopTime=10;

    case 'furnobs'

        % desired sliding motion poles
        psm=[-0.025,-0.03+1i*0.025,-0.03-1i*0.025];

        [~,~,~,~,r]=outfor(A,B,C);
        if n-p-r>0
            [G,F]=wzobs(A,B,C,psm,-0.2);
        else
            [G,F]=wzobs(A,B,C,[],-0.2);
        end

        Gamma=diag([-0.025]);

        Phi=-0.1;

        nn=n;pp=p;mm=m;
        %----------------------------------------------------------------------------%
        % Augment the statespace with integral action states
        %----------------------------------------------------------------------------%
        AT=[zeros(pp,pp) -C; zeros(nn,pp) A];
        A11a=AT(1:nn,1:nn);
        A12a=AT(1:nn,nn+1:nn+pp);

        Co=ctrb(A11a,A12a);
        unco = length(A) - rank(Co);

        if unco==0
            Tr=eye(nn);
        else
            [~,~,~,~,~,Tr]=regfor(A,B);
            %BT=[zeros(pp,mm); B];
        end

        %----------------------------------------------------------------------------%
        % Change coordinates so that regular form is achieved
        %----------------------------------------------------------------------------%
        Trtilde=[eye(pp) zeros(pp,nn);zeros(nn,pp) Tr];
        ATnew=Trtilde*AT*Trtilde';
        A11a=ATnew(1:nn,1:nn);
        A12a=ATnew(1:nn,nn+1:nn+pp);

        M=place(A11a,A12a,psm);
        S=[M eye(mm)]*Trtilde';

        % Design of observer based controller with integral action
        [L,Lr,Lrdot,Sr,Lam,P]=contlia(A,B,C,S,Phi);

        x0=zeros(1,n); % [xr, x11, x12]
        x0(1)=0;
        z0=zeros(1,n);

        rho=1;
        delta=1.0;%0.001;

        SimStopTime=500;

    case 'dcmotor'
        % reduced dc motor model with speed and current only
        A(2,2)=-1e-2; % adding friction -Bm/Jm
        A=A([2:end],[2:end]);
        B=B([2:end],:);
        C=C([2:end],[2:end]);
        [n,m]=size(B);
        [p,n]=size(C);
        D=zeros(p,m);
        C=C(:,:); % current and speed measurement
        [p,n]=size(C);
        [pp,nn]=size(C); % n-state vector length, p-number of outputs

        Q=eye(n);
        [~,E]=lqcfCA(A,B,Q);

        %%%%%%%%%%%%%%%%%%%%%%%%%% Observer Design %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % move the observer sliding mode poles further to the left
        % multiplied by factor
        psm=[];
        po_factor=5;
        Eid=1;
        if ~isempty(E)
            for i=1:(nn-pp) % ToDo check nn-pp
                % check if the poles are conjugated complex
                if ~isreal(E(Eid))
                    psm(i)=po_factor*real(E(Eid));
                else
                    psm(i)=po_factor*E(Eid);
                end

                if length(E)>(nn-pp)
                    Eid=Eid+2;
                else
                    Eid=Eid+1;
                end
            end
        end

        % check the length of psm
        [~,~,~,~,r]=outfor(A,B,C);
        if (nn-pp-r)>0
            assert(length(psm)==(nn-pp-r));
        else
            assert(length(psm)==(nn-pp)); % ToDo check
        end

        % set the estimation error poles
        per=-10*ones(1,pp);
        assert(length(per)==pp);

        [G,F]=wzobs(A,B,C,psm,per);

        % Design diagonal weighting matrix for the state vector
        Q=diag([5 5 1]);

        Gamma=diag([-10.0]);

        Phi=-0.1;

        Cc=[1 0]; % speed control
        Tc=Cc*inv(eye(n));
        Ramp=1000;

        [Ai,Bi]=intac(A,B,Cc);

        Qi = eye(length(Ai));

        [S,E]=lqcfCA(Ai,Bi,Qi);
        Eopt=E;
        %Eopt(Eopt>-1)=-1;
        S=rpp(Ai,Bi,Eopt);

        Phi=-eye(size(S,1));

        [L,Lr,Lrdot,Sr,Lam,P]=contliaCA(A,B,Cc,S,Phi);
        Lam_inv=inv(Lam);

        pm=size(Cc,1);

        x0=zeros(1,n); % [xr, x11, x12]
        x0(1)=0;
        z0=zeros(1,n);

        rho=1;
        delta=0.001;

        SimStopTime=20;

end




%% Simulate the model
mdl_name='obsint_mdl';
if ~bdIsLoaded(mdl_name)
    open([mdl_name,'.slx']);
end
set_param(mdl_name,'StopTime',sprintf('%d',SimStopTime));
sim(mdl_name);

%% plot simulation results

subplot(2,3,1)
plot(t.Data,u.Data);
grid on;
legend(get_legend('u'));

subplot(2,3,2)
plot(t.Data,s.Data);
grid on;
legend(get_legend('s'));

subplot(2,3,3)
plot(t.Data,r.Data);
hold on;
plot(t.Data,ym.Data);
grid on;
legend([get_legend('r');get_legend('ym')]);

subplot(2,3,4)
plot(t.Data,x.Data);
grid on;
legend(get_legend('x'));

subplot(2,3,5)
plot(t.Data,y.Data);
grid on;
legend(get_legend('y'));

subplot(2,3,6)
plot(t.Data,ey.Data);
grid on;
legend(get_legend('ey'));

sgtitle([mdl_name,' - ',mat_name],'Interpreter','None');

%%
function sys_info(sys, mat_name)
    disp(['Model: ',mat_name]);
    Co=ctrb(sys.A,sys.B);
    unco = length(sys.A) - rank(Co);
    disp(['Uncontrollable states: ', num2str(unco)]);

    Ob=obsv(sys.A,sys.C);
    unob = length(sys.A) - rank(Ob);
    disp(['Unobservable states: ', num2str(unob)]);

    [~,~,~,~,r]=outfor(sys.A,sys.B,sys.C);

    disp(['Invariant zeros: ', num2str(r)]);

%     fprintf('%s', 'A=[');
%     fmt = [repmat('\t%.3f', 1, size(sys.A,2)-1), '\t%.3f\n'];
%     fprintf(fmt, sys.A.');
%     fprintf('%s', ']');

    %disp(['The result is: [' num2str(sys.A) ']']) ;
    %fprintf(num2str(sys.A));
end