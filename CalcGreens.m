function Gf=CalcGreens(d,n,Al)
    syms x a;
    assume(x>=-1)
    assume(x<=1)
    G(x,a)=Greens(Al);
    f(x)=fmaker(d,n); 
    Gf=int(f(a)*G(x,a),a,-1,1);
end