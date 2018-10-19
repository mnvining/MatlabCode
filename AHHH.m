function AHHH(NG,cl,n,NoMo)
syms x

xx=linspace(-1,1,NG+1);
xx=xx(1:end-1);
h=xx(2)-xx(1);
x2=-1:h:1+cl-h;
xcoarse1=linspace(-1,1,n);
xcoarse=xcoarse1(1:end-1);
hh=xcoarse(2)-xcoarse(1);
xcoarseext=-1:hh:1+cl-hh;
u=cos(x);
f=cos(x);
ef=double(subs(f,xx));
efff=double(subs(f,xcoarse));
M=floor(NoMo/2);

ZZ=cl*NG/2+NG;
Modes=[1:M,ZZ-M+1:ZZ];
length(Modes)

A=createFWDcont(NG,Modes,cl);

fd=(cl+2)/2;%+cl/(NG-1);
ZZ=round(cl*NG/2+NG);
LL=round(ZZ/2);

K=-pi*1i*[0:LL-1,0,-(LL-1):-1]';

% derivative coefficients pertaining to the FC
p=[K(1:M);K((ZZ)-M+1:(ZZ))]/fd;
size(p)

D=diag(p);
[Ud,Sd,Vd]=svd((D));
Sd_Dag=zeros(size(Sd));
Sd_Dag(Sd>0)=1./Sd(Sd>0);
D_Dag=Vd*(Sd_Dag'*Ud');

% Apply differential operator to cosine. 
Cin=A\ef';
%Diffop=real(A*D_Dag'*Cin);
Diffop1=real(A*D'*Cin);
figure(12)

plot(xx,Diffop1,xx,-sin(xx))

%%%acc=norm(Diffop+sin(xx'))
acc1=norm(Diffop1+sin(xx'))
