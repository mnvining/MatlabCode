function Gf=totalGreen2(d,n,Al,Fine_Grid)
% Note: GreensInt is already scaled by Alpha
p=Fine_Grid(end);
    syms x;
    expr=0*x;
    zz=coeffcalc(d,n);
    for j=1:length(zz)
        m(x)=GreensIntInc(j-1,Al,p);
        expr=expr+zz(j)*m(x);
    end
    Gf= expr;
end