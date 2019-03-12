function [d1,stab,F1,ef]=RUNToo(NG,Al,n,d,cl,NoMo,opt)

syms x

h1=2/(n-1);
xcoarse=-1:h1:1-h1;
h2=h1/((NG/(n-1))+1);
xx=-1:h2:1-h1;
LD=round(NG/(n-1)+1);


f=fmaker(d,n);
ef=double(subs(f,xx));
efff=double(subs(f,xcoarse));

M=floor(NoMo/2);% this is x2 because we need double modes

ZZ=round(NG*(cl+2)/2); % this is x2 because of repeated dom.
A=dftmtx(ZZ);
A=A';
A=A/sqrt(ZZ);
if mod(NoMo,2)==0
    Modes=[1:M,ZZ-M+1:ZZ];
else
    Modes=[1:M,(ZZ-M):ZZ];
end

%% Fix: modes are off so you aren't deleting the right modes
%%


% cont matrix to eval at all pts 
A=A(:,Modes);

% cont matrix to eval on only pts where f is defined
ACont1=A([1:NG],:); % first part

LL=round(ZZ/2);

% This creates the derivative

K=-pi*1i*[0:LL-1,0,-(LL-1):-1]';
fd=(cl+2)/2;
if mod(NoMo,2)==0
    p=[K(1:(M));K((ZZ)-M+1:(ZZ))];
else   
    p=[K(1:(M));K((ZZ)-M:(ZZ))];
end
p=p/fd;

D=(eye(length(p))-Al*diag(p.*p));

% calculate B(f)
[EGGf,Eh1f,Eh2f]=BCoeffCalc2(xx,xcoarse,Al,d);

% this calculates D_Dag for use in operator
[Ud,Sd,Vd]=svd(D);
Sd_dag=pinv(Sd,1e-13);
D_Dag=Vd*Sd_dag*Ud';

%AugMat=[ACont1,zeros(NG,2);ACont2,zeros(NG,2);ACont1*D_Dag,Eh1f',Eh2f';ACont2*D_Dag,(-1)^opt*Eh1f',(-1)^opt*Eh2f'];
AugMat=[0*ACont1,zeros(NG,2);ACont1*D_Dag,Eh1f',Eh2f'];

size(AugMat)

if opt==0
    % delete odd modes
    if mod(NoMo,2)==0
        % even # of modes
        AugMat(:,[2:2:M,M+1:2:end-2])=[];
    else
       % odd # of modes
        AugMat(:,[2:2:M+1,M+2:2:end-2])=[];
    end
else
    % delete even modes
    if mod(NoMo,2)==0
        % even # of modes
        AugMat(:,[1:2:M,M+2:2:end-2])=[];
    else
        % odd # of modes
        AugMat(:,[1:2:M,M+3:2:end-2])=[];
    end
end



AugVec=[0*ef';EGGf'];


% [U,S,V]=svd(AugMat);
% Sdag=pinv(S,1e-13);
% Res=V*(Sdag*(U'*AugVec));

Res=AugMat\AugVec;


Cin=Res(1:end-2);
Res(end-1:end)
if opt == 0
    if mod(NoMo,2)==0
        A(:,[2:2:M,M+1:2:end])=[];
    else
        A(:,[2:2:M+1,M+2:2:end])=[];
    end
else
    if mod(NoMo,2)==0
        A(:,[1:2:M,M+2:2:end])=[];
    else
        A(:,[1:2:M,M+3:2:end])=[];
    end
end

if opt == 0
    if mod(NoMo,2)==0
        ACont1(:,[2:2:M,M+1:2:end])=[];
    else
        ACont1(:,[2:2:M+1,M+2:2:end])=[];
    end
else
    if mod(NoMo,2)==0
        ACont1(:,[1:2:M,M+2:2:end])=[];
    else
        ACont1(:,[1:2:M,M+3:2:end])=[];
    end
end
if opt == 0
    if mod(NoMo,2)==0
        D_Dag([2:2:M,M+1:2:end],:)=[];
    else
        D_Dag([2:2:M+1,M+2:2:end],:)=[];
    end
else
    if mod(NoMo,2)==0
        D_Dag([1:2:M,M+2:2:end],:)=[];
    else
        D_Dag([1:2:M,M+3:2:end],:)=[];
    end
end

if opt == 0
    if mod(NoMo,2)==0
        D_Dag(:,[2:2:M,M+1:2:end])=[];
    else
        D_Dag(:,[2:2:M+1,M+2:2:end])=[];
    end
else
    if mod(NoMo,2)==0
        D_Dag(:,[1:2:M,M+2:2:end])=[];
    else
        D_Dag(:,[1:2:M,M+3:2:end])=[];
    end
end






F1=real(A*Cin);
%d1=norm(F1(1:NG)-ef(1:NG)','inf')

%  U1=real(A*D_Dag'*Cin);

U1_1=real(ACont1*D_Dag*Cin);%+Res(end-1)*Eh1f'+Res(end)*Eh2f');
U1a=real(ACont1*D_Dag*Cin+Res(end-1)*Eh1f'+Res(end)*Eh2f');
%     U1_2=real(ACont2*D_Dag*Cin+(-1)^opt*Res(end-1)*Eh1f'+(-1)^opt*Res(end)*Eh2f');
%
LD=round(NG/(n-1)+1);
solcoarse=U1_1([1:LD:NG]);


%norm(F1(1:NG)-real(ACont1*Cin),'inf')
fG=ef';

d1=norm(F1(1:NG)-fG(1:NG),'inf')
d2=norm(U1a-EGGf','inf')
%
%
%
stab=norm(solcoarse)/norm(efff);
ef=ef';

end