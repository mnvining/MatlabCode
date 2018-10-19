
function [stab]=RUNGraphOdd(NG,Al,n,d,cl,opt,NoMo)
    % uses NG grid points (int)
    % uses Al for alpha - should be a big float
    % n for number of coarse grid pts on -1,1 (int)
    % d for degree of gram poly (int)
    % cl continuation length (extended domain is [-1,1+cl)- cl must be type Int
    % opt: 1 for func, 2 for 2nd der, 3, for sobolev
    
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
    ef2=double(subs(-f,xx));
    ef=[ef ef];
    efforplot=double(subs(f,xx));
    ef2forplot=double(subs(-f,xx));
    efff=double(subs(f,xcoarse));
    nfc=double(subs(-f,xcoarse)); 
    tots=[efff,efff];
    
    M=floor(NoMo/2);
  
    ZZ=(cl*NG/2+NG)*2;
    A=dftmtx(ZZ);
    A=A';
    A=A/sqrt(ZZ);
    Modes=[1:(M),(ZZ-M+1):ZZ];
    length(Modes)
    % cont matrix to eval at all
    A=A(:,Modes);
    % cont matrix to eval on only pts where f is defined
    ACont1=A([1:NG],:); % first part
    ACont2=A([ZZ/2+1:ZZ/2+NG],:); % second part
    
    LL=ZZ/2;
    % This creates the derivative
    K=-pi*1i*[0:LL-1,0,-(LL-1):-1]';
    fd=(cl+cl/(2*(NG-1)))*2;
    p=[K(1:(M));K((ZZ)-M+1:(ZZ))]/fd;
    D=(eye(size(p))-diag(p.*p));

    
    % calculate B(f)
    [cobb,EGGf,Eh1f,Eh2f]=BCoeffCalc(xx,xcoarse,Al,d,opt);

    BB=EGGf+cobb(1)*Eh1f+cobb(2)*Eh2f;

%     D=createId2Derivative(NG,M,Al,cl);

    
    % this calculates D_Dag for use in operator    
    [Ud,Sd,Vd]=svd(D);
    Sd_dag=1./Sd(Sd>=1e-14);
    Sd_dag=diag(Sd_dag);
    D_Dag=Vd*Sd_dag*Ud';
    lambdam=1;%1e15;

    AugMat=[ACont1;ACont2;lambdam*ACont1*D_Dag';lambdam*ACont2*D_Dag'];



    AugVec=[ef';lambdam*BB';lambdam*BB'];

    Cin=AugMat\AugVec;
    %Cin=[ACont1;ACont2]\ef';
    size(xall)
    
    
    F1=real(A*Cin);
    figure(1)
    plot(xall,F1,'b',xx,efforplot,'k',x2,efforplot,'r')
    figure(2)
    plot(xall,F1,'b',xx,efforplot,'k',x2,efforplot,'r')
    axis([-1 1 -5 5])
    figure(3)
    U1=real(A*D_Dag'*Cin);
    plot(xall,U1,'b',xx,BB,'k',x2,-BB,'r')
    hold on
    axis([-1 1 -5 5])
    % coarse grid things
    LD=round(NG/(n-1));
    solcoarse=U1([1:LD:NG]);
    plot(xcoarse,solcoarse,'ro','MarkerSize',5)
    figure(2)
    hold on
    plot(xcoarse,F1(1:LD:NG),'ro','MarkerSize',6)
  
    stab=norm(solcoarse)/norm(efff)
    stab2=norm([BB(1:LD:end),-BB(1:LD:end)])/norm(tots);

end
