function [A, Bv, Tr, CA]=SMCCanForm(A, B)

[nn,mm]=size(B); % B = [4x8]


%--------------------------------------------------------------------------%
% Perform QR decomposition on the input distribution matrix
%--------------------------------------------------------------------------%
[Tr,temp]=qr(B); % Tr = [4x4]
Tr=flip(Tr');


%--------------------------------------------------------------------------%
% Obtain (Areg,Breg); regular form description
%--------------------------------------------------------------------------%
A=Tr*A*Tr'; % A = [4x4]
B=Tr*B; % B = [4x8]

BB=round(B/(max(max(abs(B)))), 2); % BB = [4x8]

ll=rank(BB); % ll = 2
%--------------------------------------------------------------------------%
% Obtain matrix sub-blocks for sliding mode controller design
%--------------------------------------------------------------------------%
A11 = A(1:nn-ll,1:nn-ll);  % A11 = [2x2]
A12 = A(1:nn-ll,nn-ll+1:nn); % A12 = [2x2]
B2 = B(nn-ll+1:nn,1:end); % B2 = [2x8]

Bv=[zeros(nn-ll ,ll)
    eye(ll)]; % Bv = [4x2]

CA=pinv(B2); % Ca = [8x2]

end
