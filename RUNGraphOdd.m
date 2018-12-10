
function [F1,Cin]=RUNGraphOdd(NG,Al,n,d,cl,NoMo,opt)
    % uses NG grid points (int)
    % uses Al for alpha - should be a big float
    % n for number of coarse grid pts on -1,1 (int)
    % d for degree of gram poly (int)
    % cl continuation length (extended domain is [-1,1+cl)- cl must be type Int
    % opt: 0 for even, 1 for odd
    
    syms x

    xx=linspace(-1,1,NG+1);
    xx=xx(1:end-1);
    h=xx(2)-xx(1);
    x2=6:h:1+cl+2-h;
    xall=-1:h:1+2*cl+2-h;
    length(xall)
    xcoarse1=linspace(-1,1,n);
    xcoarse2=linspace(6,8,n);
    xcoarse=xcoarse1(1:end-1);
    xcoarse2=xcoarse2(1:end-1);
    hh=xcoarse(2)-xcoarse(1);
    xcoarseext=-1:hh:1+cl+2-hh;
    f=fmaker(d,n);
    ef=double(subs(f,xx));
    ef2=double(subs((-1)^opt*f,xx));
    ef=[ef ef2];
    efforplot=double(subs(f,xx));
    ef2forplot=double(subs((-1)^opt*f,xx));
    efff=double(subs(f,xcoarse));
    nfc=double(subs((-1)^opt*f,xcoarse)); 
    tots=[efff,nfc];
    
    M=2*floor(NoMo/2);% this is x2 because we need double modes
  
    ZZ=2*(cl*NG/2+NG); % this is x2 because of repeated dom.
    A=dftmtx(ZZ);
    A=A';
    A=A/sqrt(ZZ);
    Modes=[1:(M+1),(ZZ-M+1):ZZ];
    % cont matrix to eval at all
    A=A(:,Modes);
    % cont matrix to eval on only pts where f is defined
    ACont1=A([1:NG],:); % first part
    ACont2=A([ZZ/2+1:ZZ/2+NG],:); % second part
    
    LL=round(ZZ/2);
    % This creates the derivative
    K=-pi*1i*[0:LL-1,0,-(LL-1):-1]';
    fd=2*(cl+2)/2;
    p=[K(1:(M+1));K((ZZ)-M+1:(ZZ))]/fd;
    D=(eye(length(p))-Al*diag(p.*p));
    
    % calculate B(f)
    [EGGf,Eh1f,Eh2f]=BCoeffCalc(xx,xcoarse,Al,d);

    % this calculates D_Dag for use in operator    
    [Ud,Sd,Vd]=svd(D);
    Sd_dag=pinv(Sd,1e-13);
    D_Dag=Vd*Sd_dag*Ud';

    AugMat=[ACont1,zeros(NG,2);ACont2,zeros(NG,2);ACont1*D_Dag',Eh1f',Eh2f';ACont2*D_Dag',(-1)^opt*Eh1f',(-1)^opt*Eh2f'];
    
    AugVec=[ef';EGGf';(-1)^opt*EGGf'];

    
    [U,S,V]=svd(AugMat);
    Sdag=pinv(S,1e-13);
    Res=V*(Sdag*(U'*AugVec));
    
    Cin=Res(1:end-2);
    size(Cin)
    size(A)
    
    F1=real(A*Cin);
    
    U1=real(A*D_Dag'*Cin);
    U1_1=real(ACont1*D_Dag'*Cin+Res(end-1)*Eh1f'+Res(end)*Eh2f');
    U1_2=real(ACont2*D_Dag'*Cin+(-1)^opt*Res(end-1)*Eh1f'+(-1)^opt*Res(end)*Eh2f');
    size(U1_1)
    %plot(xx,U1_1,'b',x2,U1_2,'b')%,xx,EGGf,'k',x2,-EGGf,'r')
    %hold on
    %axis([-1 13 -5 5])
    % coarse grid things
    LD=round(NG/(n-1));
    solcoarse=U1_1([1:LD:NG]);
    %plot(xcoarse,solcoarse,'ro','MarkerSize',5)
    %figure(2)
    %hold on
    %plot(xcoarse,F1(1:LD:NG),'ro','MarkerSize',6)
    
    norm(F1(1:NG)-real(ACont1*Cin))
    
    d1=norm(F1(1:NG)-ef(1:NG)');
    d2=norm(U1_1-EGGf');
    
    

    stab=norm(solcoarse)/norm(efff);
    %stab2=norm([EGGf(1:LD:end),-EGGf(1:LD:end)])/norm(tots);
    %stab=1;

end
