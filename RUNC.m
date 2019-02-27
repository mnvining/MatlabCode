function [d1,F1]=RUNC(NG,Al,n,d,cl,NoMo)
% uses NG grid points (int)
% uses Al for alpha - should be a big float
% n for number of coarse grid pts on -1,1 (int)
% d for degree of gram poly (int)
% cl continuation length (extended domain is [-1,1+cl)- cl must be type Int
% opt: 0 for even, 1 for odd

syms x
h=2/(n-1);
xx=linspace(-1,1-h,NG);

f=fmaker(d,n);
ef=double(subs(f,xx));
efforplot=double(subs(f,xx));

M=floor(NoMo/2);% this is x2 because we need double modes

ZZ=(cl*NG/2+NG); % this is x2 because of repeated dom.
A=dftmtx(ZZ);
A=A';
A=A/sqrt(ZZ);
Modes=[1:(M+1),(ZZ-M+1):ZZ];
% cont matrix to eval at all
A=A(:,Modes);
ACont1=A(1:NG,:);

AugMat=ACont1;



AugVec=[ef'];


[U,S,V]=svd(AugMat);
Sdag=pinv(S,1e-13);
Res=V*(Sdag*(U'*AugVec));
Cin=Res;

F1=real(A*Cin);


d1=norm(F1(1:NG)-ef(1:NG)','inf');

end
