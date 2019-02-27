function D=createDer2(ZZ,cl,M,Al);
LL=round(ZZ/2);
% This creates the derivative
K=-pi*1i*[0:LL-1,0,-(LL-1):-1]';
fd=2*(cl+2)/2;
p=[K(1:(M+1));K((ZZ)-M+1:(ZZ))]/fd;
D=(eye(length(p))-Al*diag(p.*p));