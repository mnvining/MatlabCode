function co=coeffcalc(d,n)
    syms x
    exx=fmaker(d,n);
    co=double(coeffs(exx));
end