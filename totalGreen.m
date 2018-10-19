function Gf=totalGreen(d,n,Al)
    syms x;
    expr=0*x;
    zz=coeffcalc(d,n);
    for j=1:length(zz)
        m(x)=GreensInt(j-1,Al);
        expr=expr+zz(j)*m(x);
    end
    Gf= -expr/Al;
end