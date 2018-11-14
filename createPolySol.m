function PS=createPolySol(d,n,xx,Al)
% function PS=createPolySol(d,n,xx)
% inputs: d - degree of poly desired
%           n - number of points on coarse grid +1
%           xx - fine grid [row vector for consistency w/ RUNPolys code]
%           Al - alpha

    z=coeffcalc(d,n);
    %z=fliplr(z);
    xx=xx';
    
    if d==0
        PS=z(1)*xx.^0;
    elseif d==1
        PS=z(1)*xx+z(2);
    elseif d==2
        PS=z(1)*xx.^2+z(2)*xx+(z(3)+2*z(1)*Al);
    elseif d==3
        PS=z(1)*xx.^3+z(2)*xx.^2+(z(3)+6*z(1)*Al)*xx+(z(4)+2*z(2)*Al);
    elseif d==4
        PS=z(1)*xx.^4+z(2)*xx.^3+(z(3)+12*z(1)*Al)*xx.^2+(z(4)+6*z(2)*Al)*xx+(z(5)+2*z(3)*Al+24*z(1)*(Al^2));
    elseif d==5
        PS=z(1)*xx.^5+z(2)*xx.^4+(z(3)+20*z(1)*Al)*xx.^3+(z(4)+12*z(2)*Al)*xx.^2+(z(5)+6*z(3)*Al+120*z(1)*(Al^2))*xx+(z(6)+2*z(4)*Al + 24*z(2)*(Al^2));
    elseif d==6
        PS=z(1)*xx.^6+z(2)*xx.^5+(z(3)+30*z(1)*Al)*xx.^4+(z(4)+20*z(2)*Al)*xx.^3+(360*z(1)*(Al^2)+12*z(3)*Al+z(5))*xx.^2+(120*z(2)*(Al^2)+6*Al*z(4)+z(6))*xx+(720*z(1)*(Al^3)+24*z(3)*(Al^2)+2*Al*z(5)+z(7));
    else
        PS=z(1)*xx.^7+z(2)*xx.^6+(z(3)+42*z(1)*Al)*xx.^5+(z(4)+30*z(2)*Al)*xx.^4+(840*z(1)*(Al^2)+20*z(3)*Al+z(5))*xx.^3+(360*z(2)*(Al^2)+12*Al*z(4)+z(6))*xx.^2+(5040*z(1)*(Al^3)+120*z(3)*(Al^2)+6*z(5)*Al+z(7))*xx+(720*z(2)*(Al^3)+24*z(4)*(Al^2)+2*z(6)*Al+z(8));
    end
    
        
        
        