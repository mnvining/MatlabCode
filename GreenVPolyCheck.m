
syms x
xx=linspace(-1,1,257);
xx=xx(1:end-1);
xcoarse1=-1:.25:.75;
d=0;
n=9;
Al=.01;
f(x)=x.^d;%=fmaker(d,n);
ef=double(subs(f,xx));

PSOL=PlainPolySol(d,n,xx,Al);

GI(x)= GreensInt(d,Al)/Al;


[~,~,Eh1f,Eh2f]=BCoeffCalc(xx,xcoarse1,Al,d,1);
EGGf=double(subs(GI,xx));



Mat=[Eh1f',Eh2f'];
B=EGGf'-PSOL;
co=Mat\B;

plot(xx,abs(EGGf'-co(1)*Eh1f'-co(2)*Eh2f'-PSOL))

