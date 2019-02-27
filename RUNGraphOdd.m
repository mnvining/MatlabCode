
function [d1,stab,F1,ef]=RUNGraphOdd(NG,Al,n,d,cl,NoMo,opt)
% uses NG grid points (int)
% uses Al for alpha - should be a big float
% n for number of coarse grid pts on -1,1 (int)
% d for degree of gram poly (int)
% cl continuation length (extended domain is [-1,1+cl)- cl must be type Int
% opt: 0 for even, 1 for odd

syms x

h1=2/(n-1);
xcoarse=-1:h1:1-h1;
h2=h1/((NG/(n-1)+1));
xx=-1:h2:1-h1;
size(xx)
LD=round(NG/(n-1)+1);



f=fmaker(d,n);
ef=double(subs(f,xx));
efff=double(subs(f,xcoarse));


[A,Modes,ZZ,M]=CreateAMat(NoMo,cl,NG);
% cont matrix to eval on only pts where f is defined
ACont1=A([1:NG],:); % first part

D=createDer2(ZZ,cl,M,Al);

[EGGf,Eh1f,Eh2f]=BCoeffCalc2(xx,xcoarse,Al,d);

% this calculates D_Dag for use in operator
[Ud,Sd,Vd]=svd(D);
Sd_dag=pinv(Sd,1e-13);
D_Dag=Vd*Sd_dag*transpose(Ud);

AugMat=[ACont1,zeros(NG,2);ACont1*D_Dag,Eh1f',Eh2f'];


if opt==0
    if mod(NoMo,2)==0
        AugMat(:,[2:2:NoMo,NoMo+1:2:end-2])=[];
    else
        AugMat(:,[2:2:NoMo,NoMo+2:2:end-2])=[];
    end
else
    if mod(NoMo,2)==0
        AugMat(:,[1:2:NoMo,NoMo+2:2:end-2])=[];
    else
        AugMat(:,[1:2:NoMo,NoMo+1:2:end-2])=[];
    end
end



AugVec=[ef';EGGf'];


[U,S,V]=svd(AugMat);
Sdag=pinv(S,1e-13);
Res=V*(Sdag*(transpose(U)*AugVec));
size(Res)
Cin=Res(1:end-2);
if opt == 0
    if mod(NoMo,2)==0
        A(:,[2:2:NoMo,NoMo+1:2:end])=[];
    else
        A(:,[2:2:NoMo,NoMo+2:2:end])=[];
    end
else
    if mod(NoMo,2)==0
        A(:,[1:2:NoMo,NoMo+2:2:end])=[];
    else
        A(:,[1:2:NoMo,NoMo+1:2:end])=[];
        size(A)
    end
end

if opt == 0
    if mod(NoMo,2)==0
        ACont1(:,[2:2:NoMo,NoMo+1:2:end])=[];
    else
        ACont1(:,[2:2:NoMo,NoMo+2:2:end])=[];
    end
else
    if mod(NoMo,2)==0
        ACont1(:,[1:2:NoMo,NoMo+2:2:end])=[];
    else
        ACont1(:,[1:2:NoMo,NoMo+1:2:end])=[];
    end
end
if opt == 0
    if mod(NoMo,2)==0
        D_Dag([2:2:NoMo,NoMo+1:2:end],:)=[];
    else
        D_Dag([2:2:NoMo,NoMo+2:2:end],:)=[];
    end
else
    if mod(NoMo,2)==0
        D_Dag([1:2:NoMo,NoMo+2:2:end],:)=[];
    else
        D_Dag([1:2:NoMo,NoMo+1:2:end],:)=[];
    end
end

if opt == 0
    if mod(NoMo,2)==0
        D_Dag(:,[2:2:NoMo,NoMo+1:2:end])=[];
    else
        D_Dag(:,[2:2:NoMo,NoMo+2:2:end])=[];
    end
else
    if mod(NoMo,2)==0
        D_Dag(:,[1:2:NoMo,NoMo+2:2:end])=[];
    else
        D_Dag(:,[1:2:NoMo,NoMo+1:2:end])=[];
    end
end





F1=real(A*Cin);
%d1=norm(F1(1:NG)-ef(1:NG)','inf')

%  U1=real(A*D_Dag'*Cin);

U1_1=real(ACont1*D_Dag*Cin);%+Res(end-1)*Eh1f'+Res(end)*Eh2f');
%     U1_2=real(ACont2*D_Dag*Cin+(-1)^opt*Res(end-1)*Eh1f'+(-1)^opt*Res(end)*Eh2f');
%

solcoarse=U1_1([1:LD:NG]);

d1=norm(F1(1:NG)-ef(1:NG)','inf');
%     d2=norm(U1_1-EGGf');
%
%
%
stab=norm(solcoarse)/norm(efff);
ef = ef'

end
