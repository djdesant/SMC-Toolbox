function [L, Lr, P, Lam, S, E, CA]=SMCDesign_initIA2(A, B, C, Q)

[At, Bv, Tr, CA]=SMCCanForm(A, B);

Ct=C*inv(Tr);

[Ai,Bvi]=intac(At,Bv,Ct); 

[A, Bv, Tr, CA]=SMCCanForm(Ai, Bvi)

[S,E]=lqcfCA(Ai,Bvi,Q);

Phi=-eye(size(S,1));

% [L,Lr,Lrdot,Sr,Lam,P]=contliaCA(At,Bv,Ct,S,Phi);
[L,Lr,Lrdot,Sr,Lam,P]=contliaCA2(At,Bv,Ct,S,Phi);

[nn ll]=size(Bv);
TT=[eye(size(S,1)) zeros(ll, nn)
    zeros(nn, ll)  Tr];

L=L*TT;
S=S*TT;

end
