function [cobb,EGGf,Eh1f,Eh2f]=BCoeffCalc(Fine_Grid,Coarse_Grid,Al,d,opt)
syms x;
n=length(Fine_Grid);
GG(x)=totalGreen(d,n,Al);
[h1,h2]=hsol(Coarse_Grid,Al);
d2h1=diff(h1,2);
d2h2=diff(h2,2);
d2GG=diff(GG,2);


EGG=double(subs(GG,Coarse_Grid));
Eh1=double(subs(h1,Coarse_Grid));
Eh2=double(subs(h2,Coarse_Grid));

EG2=double(subs(d2GG,Coarse_Grid));
EH1=double(subs(d2h1,Coarse_Grid));
EH2=double(subs(d2h2,Coarse_Grid));


EGGf=double(subs(GG,Fine_Grid));
Eh1f=double(subs(h1,Fine_Grid));
Eh2f=double(subs(h2,Fine_Grid));

if opt==1
    Z=[Eh1*Eh1' Eh2*Eh1'; Eh1*Eh2' Eh2*Eh2'];
    zz=[-EGG*Eh1';-EGG*Eh2'];
    cobb=Z\zz;
    
elseif opt==2
    Z=[EH1*EH1' EH2*EH1'; EH1*EH2' EH2*EH2'];
    zz=[-EG2*EH1';-EG2*EH2'];
    cobb=Z\zz;
    
elseif opt==3
    Z=[Eh1*Eh1'+Eh1*EH1'+EH1*Eh1'+EH1*EH1' Eh2*Eh1'+EH2*Eh1'+Eh2*EH1'+EH2*EH1'; Eh1*Eh2'+Eh1*EH2'+EH1*Eh2'+EH1*EH2' Eh2*Eh2'+EH2*Eh2'+Eh2*EH2'+EH2*EH2'];
    zz=[-(EGG*Eh1'+EG2*Eh1'+EGG*EH1'+EG2*EH1');-(EGG*Eh2'+EG2*Eh2'+EGG*EH2'+EG2*EH2')];
    cobb=Z\zz;
end

end