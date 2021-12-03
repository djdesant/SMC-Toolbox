function [L, P, Lam_inv, S, E, CA,...
    Ao, Bo, Co, To_inv, G, F, Gl, Gn, P2]=SMCDesign_init(A, B, C, Q)

%%%%%%%%%%%%%%%%%%%%%%%%%%% Canonical Form %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[At, Bv, Tr, CA]=SMCCanForm(A, B);

%%%%%%%%%%%%%%%%%%%%%%%%%% Controller Design %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ct=C*inv(Tr);

[S,E]=lqcfCA(At,Bv,Q);

Phi=-eye(size(S,1));

[L,P,Lam]=contl(At,Bv,S,Phi);
Lam_inv=inv(Lam);

[nn,ll]=size(Bv);
[pp,nn]=size(Ct);

L=L*Tr;
S=S*Tr;

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
