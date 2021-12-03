function [L, Lr, Lrdot, P, Lam_inv, S, Sr, E, CA,...
    Ao, Bo, Co, To_inv, G, F, Gl, Gn, P2]=SMCDesign_initIA(A, B, C, Q)

%%%%%%%%%%%%%%%%%%%%%%%%%%% Canonical Form %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[At, Bv, Tr, CA]=SMCCanForm(A, B);

%%%%%%%%%%%%%%% Controller Design wiht Integration Action %%%%%%%%%%%%%%%%%
Ct=C*inv(Tr);

[Ai,Bvi]=intac(At,Bv,Ct);               

[S,E]=lqcfCA(Ai,Bvi,Q);

Phi=-eye(size(S,1));

[L,Lr,Lrdot,Sr,Lam,P]=contliaCA2(At,Bv,Ct,S,Phi);
Lam_inv=inv(Lam);
% [L,Lr,Lrdot,Sr,Lam,P]=contliaCA(At,Bv,Ct,S,Phi);

[nn,ll]=size(Bv);
[pp,nn]=size(Ct);

TT=[eye(size(S,1)) zeros(ll, nn)
    zeros(nn, ll)  Tr];

L=L*TT;
S=S*TT;

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
% %             Eid=Eid+2;
        else
            psm(i)=po_factor*E(Eid);
        end
        Eid=Eid+1;
    end
end

% check the length of psm
[~,~,~,~,r]=outfor(At,Bv,Ct);
if (nn-pp-r)>0
    assert(length(psm)==(nn-pp-r));
else
    assert(length(psm)==(nn-pp)); % ToDo check
end

% set the estimation error poles
per=-10*ones(1,pp);
assert(length(per)==pp);

Ao=At;
Bo=Bv;
Co=Ct;
To_inv=inv(Tr);

[G,F]=wzobs(At,Bv,Ct,psm,per);

[Gl,Gn,P2]=disobs(At,Bv,Ct,psm,per);

end
