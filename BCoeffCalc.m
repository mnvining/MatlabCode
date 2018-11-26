function [EGGf,Eh1f,Eh2f]=BCoeffCalc(Fine_Grid,Coarse_Grid,Al,d)
syms x;
n=length(Coarse_Grid);
GG(x)=totalGreen(d,n+1,Al);
[h1,h2]=hsol(Coarse_Grid,Al);

EGGf=double(subs(GG,Fine_Grid));
Eh1f=double(subs(h1,Fine_Grid));
Eh2f=double(subs(h2,Fine_Grid));

end