
function [stab,acc,acc2]=RUNGraph(NG,Al,n,d,cl,opt,NoMo,lambda)
    % function [stab,acc,acc2]=RUNGraph(NG,Al,n,d,cl,opt,NoMo,lambda)
    % uses NG grid points (int)
    % uses Al for alpha - should be a big float
    % n for number of coarse grid pts on -1,1 (int)
    % d for degree of gram poly (int)
    % cl continuation length (extended domain is [-1,1+cl)- cl must be type Int
    % opt: 1 for func, 2 for 2nd der, 3, for sobolev
    % NoMo number of modes desired [max for d=8 is 28]
    % lambda = weighting of function vs derivative in matrix 
    
    syms x

    xx=linspace(-1,1,NG+1);
    xx=xx(1:end-1);
    h=xx(2)-xx(1);
    x2=-1:h:1+cl-h;
    xcoarse1=linspace(-1,1,n);
    xcoarse=xcoarse1(1:end-1);
    hh=xcoarse(2)-xcoarse(1);
    xcoarseext=-1:hh:1+cl-hh;
    f=fmaker(d,n);
    ef=double(subs(f,xx));
    efff=double(subs(f,xcoarse));
    figure(10)
    plot(xx,ef)
    M=floor(NoMo/2);
  
    ZZ=cl*NG/2+NG;
    Modes=[1:M,ZZ-M+1:ZZ];

    A=createFWDcont(NG,Modes,cl);
    B=createBackDoubleDomain(NG,Modes,cl);

    [cobb,EGGf,Eh1f,Eh2f]=BCoeffCalc(xx,xcoarse1,Al,d,opt);

    BB=EGGf'+cobb(1)*Eh1f'+cobb(2)*Eh2f';
    
    D=createId2Derivative(NG,M,Al,cl);
    
    % FIX THIS STUFF!!
    [Ud,Sd,Vd]=svd(D);
    Sd_Dag=zeros(size(Sd));
    Sd_Dag(Sd>1e-13)=1./(Sd(Sd>1e-13));
    D_Dag=Vd*Sd_Dag'*Ud';
    
    [Ua,Sa,Va]=svd(A);
    SaDag=zeros(size(Sa));
    SaDag(Sa>1e-13)=1./(Sa(Sa>1e-13));
    
    h1_C=Va*(SaDag'*(Ua'*Eh1f'));
    h2_C=Va*(SaDag'*(Ua'*Eh2f'));
    r1=(real(A*h1_C)-Eh1f');
    r2=(real(A*h2_C)-Eh2f');
    if norm(r1)<1e-16
        K=ef';
    else
        Mat1=[r1'*Eh1f',r1'*Eh2f';r2'*Eh1f',r2'*Eh2f'];
        cond(Mat1)
        RHS1=[-r1'*EGGf';-r2'*EGGf']; 
        Coeff=Mat1\RHS1;
        K=EGGf'+Coeff(1)*Eh1f'+Coeff(2)*Eh2f';
    end
    
    
    
    

    %AugMat=[lambda*A,zeros(NG,2);A*D_Dag',A*h1_C,A*h2_C];
    
    %cond([A;A*D_Dag'])

    AugVec=[lambda*ef';K];
    
    AugMat=[lambda*A;A*D_Dag'];
    %AugVec=[lambda*ef';BB];
    cond(AugMat)
    [U,S,V]=svd(AugMat);
    S_Dag=zeros(size(S));
    S_Dag(S>1e-13)=1./S(S>1e-13);
    S_Dag=S_Dag';
    Res=V*(S_Dag*(U'*AugVec));
    %Cin=Res(1:end-2);
 
    Cin=Res;
    
    figure(90)
    loglog(real(U'*AugVec),real(S_Dag'),'o','MarkerSize',6)

    Uex=real(B*D_Dag*Cin);
    U2=real(A*Cin);
    U1=real(A*D_Dag*Cin);%*Cin+Res(end-1)*Eh1f'+Res(end)*Eh2f');
    %U1=real(A*D_Dag*Cin);
    
    %figure(22)
    %plot(xx,U2,xx,ef,xx,BB)
    %axis([-1 1 -5 5])
    %legend('U','f','B')
    
    U2=real(A*Cin);
    U2ex=real(B*Cin);
    acc=norm(U1-EGGf');
    acc2=norm(U2-ef');
    % coarse grid things
    LD=round(NG/(n-1));
    solcoarse=U1(1:LD:end);
    %plot(xcoarse,U2(1:LD:end),'ro','MarkerSize',6)
    %plot(xcoarse,efff,'ko','MarkerSize',6)
    %figure(43)
    %semilogy(abs(U2-ef'))
    
    %figure(91)
    %semilogy(abs(U1-EGGf'))
  
    stab=norm(solcoarse)/norm(efff);

end
