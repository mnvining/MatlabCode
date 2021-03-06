
function [AugMat,A,F1,Cin]=RUNGraphTen(NG,Al,n,d,cl,NoMo,opt)
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
    x2=6:h:1+cl+2;
    xall=-1:h:1+2*cl+2;
    xcoarse1=linspace(-1,1,n);
    xcoarse2=linspace(6,8,n);
    hh=xcoarse1(2)-xcoarse1(1);
    xcoarseext=-1:hh:1+cl+2;
    f=fmaker2(d,n);
    ef=double(subs(f,xx));
    plot(xx,ef);
    ef2=double(subs((-1)^opt*f,xx));
    ef=[ef ef2];
    efforplot=double(subs(f,xx));
    ef2forplot=double(subs((-1)^opt*f,xx));
    efff=double(subs(f,xcoarse1));
    nfc=double(subs((-1)^opt*f,xcoarse2)); 
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
    [EGGf,Eh1f,Eh2f]=BCoeffCalc2(xx,xcoarse1,Al,d);

    % this calculates D_Dag for use in operator    
    [Ud,Sd,Vd]=svd(D);
    Sd_dag=pinv(Sd,1e-13);
    D_Dag=Vd*Sd_dag*Ud';
    
    
    size(ACont1)
    size(Eh1f')

    AugMat=[ACont1,zeros(NG,2);ACont2,zeros(NG,2);ACont1*D_Dag,Eh1f',Eh2f';ACont2*D_Dag,(-1)^opt*Eh1f',(-1)^opt*Eh2f'];
    
    
    if opt==0
        if mod(NoMo,2)==0
            AugMat(:,[2:2:NoMo,NoMo+1:2:end-2])=[];
        else
            AugMat(:,[2:2:NoMo,NoMo+2:2:end-2])=[];
        end
    else
        if mod(NoMo,2)==0
            AugMat(:,[1:2:NoMo,NoMo+2:2:end-2])=[];
        else
            AugMat(:,[1:2:NoMo,NoMo+1:2:end-2])=[];
        end
    end
    
    
    
    AugVec=[ef';EGGf';(-1)^opt*EGGf'];

    
    [U,S,V]=svd(AugMat);
    Sdag=pinv(S,1e-13);
    Res=V*(Sdag*(U'*AugVec));
    
    Cin=Res(1:end-2);
    if opt == 0
        if mod(NoMo,2)==0
            A(:,[2:2:NoMo,NoMo+1:2:end])=[];
        else
            A(:,[2:2:NoMo,NoMo+2:2:end])=[];
        end
    else 
        if mod(NoMo,2)==0
            A(:,[1:2:NoMo,NoMo+2:2:end])=[];
        else
            A(:,[1:2:NoMo,NoMo+1:2:end])=[];
        end
    end
    
    
        
    
    F1=real(A*Cin);
    %d1=norm(F1(1:NG)-ef(1:NG)','inf')
    
%     U1=real(A*D_Dag'*Cin);
%     U1_1=real(ACont1*D_Dag*Cin+Res(end-1)*Eh1f'+Res(end)*Eh2f');
%     U1_2=real(ACont2*D_Dag*Cin+(-1)^opt*Res(end-1)*Eh1f'+(-1)^opt*Res(end)*Eh2f');
% 
%     LD=round(NG/(n-1));
%     solcoarse=U1_1([1:LD:NG])

    
    %norm(F1(1:NG)-real(ACont1*Cin),'inf')
    
    %d1=norm(F1(1:NG)-ef(1:NG)','inf')
%     d2=norm(U1_1-EGGf');
%     
%     
% 
%     stab=norm(solcoarse)/norm(efff);

end
