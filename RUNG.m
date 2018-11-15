function [stab,acc,acc2]=RUNG(NG,Al,n,d,cl,NoMo,lambda,opt)
    % function [stab,acc,acc2]=RUNGraph(NG,Al,n,d,cl,opt,NoMo,lambda)
    % uses NG grid points (int)
    % uses Al for alpha - should be a big float
    % n for number of coarse grid pts on -1,1 (int)
    % d for degree of gram poly (int)
    % cl continuation length (extended domain is [-1,1+cl)- cl must be type Int
    % opt: 1 for func 2 for homogeneous solutions
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
%     figure(10)
%     plot(xx,ef)
%     title('function only')
    M=floor(NoMo/2);
  
    ZZ=cl*NG/2+NG;
    Modes=[1:M,ZZ-M+1:ZZ];

    A=createFWDcont(NG,Modes,cl);
    B=createBackDoubleDomain(NG,Modes,cl);

    [EGGf,Eh1f,Eh2f]=BCoeffCalc(xx,xcoarse1,Al,d);
    
    D=createId2Derivative(NG,M,Al,cl);
    
    % FIX THIS STUFF!!
    [Ud,Sd,Vd]=svd(D);
    Sd_Dag=zeros(size(Sd));
    Sd_Dag(Sd>1e-13)=1./(Sd(Sd>1e-13));
    D_Dag=Vd*Sd_Dag'*Ud';

    %PSOL=createPolySol(d,n,xx,Al);
    
    %c=[Eh1f',Eh2f']\[EGGf'-PSOL];
    %C=EGGf'-c(1)*Eh1f'-c(2)*Eh2f';
    
    

    AugVec=[lambda*ef';EGGf'];
    
    AugMat=[lambda*A,ones(NG,2); A*D_Dag',Eh1f',Eh2f'];
    %AugMat=[lambda*A;A*D_Dag'];
  
    [U,S,V]=svd(AugMat);
    S_Dag=zeros(size(S));
    S_Dag(S>1e-13)=1./S(S>1e-13);
    S_Dag=S_Dag';
    Res=V*(S_Dag*(U'*AugVec));

    Cin=Res(1:end-2);
    %Cin=Res;

 
    
    figure(90)
    loglog(real(U'*AugVec),real(S_Dag'),'o','MarkerSize',6)
    title('decay of singular vals')
    
    % U2 is fc 
    U2=real(A*Cin);
    
    
    % U1 is uc
    %U1=real(A*D_Dag'*Cin+Res(end-1)*Eh1f'+Res(end)*Eh2f');
    U1=real(A*D_Dag'*Cin);

%     
    figure(22)
    plot(xx,ef,xx,U2)
    axis([-1 1 -5 5])
    title('f vs fc')
    figure(24)
    
    plot(xx,EGGf',xx,U1)
    title('computed uc vs greens')

    acc=norm(U1-EGGf');
    acc2=norm(U2-ef');
    % coarse grid things
    LD=round(NG/(n-1));
    solcoarse=U1(1:LD:end);

%     figure(43)
%     semilogy(abs(U2-ef'))
%     title('norm f fc')
%     figure(91)
%     semilogy(abs(U1-EGGf'))
%     title('norm uc greens')
  
    stab=norm(solcoarse)/norm(efff);

end
