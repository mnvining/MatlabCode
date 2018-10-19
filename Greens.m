function G=Greens(Al)
syms x a


    A1=-(1/2)*exp(-4/sqrt(Al))*sqrt(Al)*(-exp((a-1)/sqrt(Al))+exp(-(a-1)/sqrt(Al)))/(exp(-4/sqrt(Al))-1);
    A2=(1/2)*sqrt(Al)*(-exp((a-1)/sqrt(Al))+exp(-(a-1)/sqrt(Al)))/(exp(-4/sqrt(Al))-1);
    B1=-(1/2)*sqrt(Al)*(exp(-(a-1)/sqrt(Al))*exp(-4/sqrt(Al))-exp((a-1)/sqrt(Al)))/(exp(-4/sqrt(Al))-1);
    B2=(1/2)*sqrt(Al)*(exp(-(a-1)/sqrt(Al))*exp(-4/sqrt(Al))-exp((a-1)/sqrt(Al)))/(exp(-4/sqrt(Al))-1);

G=piecewise(x<a,A1*exp((1-x)/sqrt(Al))+A2*exp((x-1)/sqrt(Al)),x>=a,B1*exp((1-x)/sqrt(Al))+B2*exp((x-1)/sqrt(Al)));

end