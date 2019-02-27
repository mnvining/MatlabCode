function [EGGf,Eh1f,Eh2f]=BCoeffCalc2(Fine_Grid,Coarse_Grid,Al,d)
syms x;
n=length(Coarse_Grid);
GG(x)=totalGreen2(d,n+1,Al,Fine_Grid);
[h1,h2]=hsol(Fine_Grid,Al);

EGGf=double(subs(GG,Fine_Grid));
Eh1f=double(subs(h1,Fine_Grid));
Eh2f=double(subs(h2,Fine_Grid));

end