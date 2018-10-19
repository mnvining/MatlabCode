
function [Opt1,Opt2]=PlsWork(NG,Al,n,cl,opt,NoMo)
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
    x2=-1:h:1+cl-h;
    xcoarse1=linspace(-1,1,n);
    xcoarse=xcoarse1(1:end-1);
    hh=xcoarse(2)-xcoarse(1);
    xcoarseext=-1:hh:1+cl-hh;
    u=cos(x);
    f=(1+Al)*cos(x);
    ef=double(subs(f,xx));
    efff=double(subs(f,xcoarse));
    M=floor(NoMo/2);
  
    ZZ=cl*NG/2+NG;
    Modes=[1:M,ZZ-M+1:ZZ];
    length(Modes)

    A=createFWDcont(NG,Modes,cl);
    B=createBackDoubleDomain(NG,Modes,cl);

    [cobb,EGGf,Eh1f,Eh2f]=BCoeffCalcCos(xx,xcoarse1,Al,opt);
    %figure(8)
    %plot(xx,cos(xx),xx,EGGf)

    BB=EGGf+cobb(1)*Eh1f+cobb(2)*Eh2f;

    D=createId2Derivative(NG,M,Al,cl);
    
    % FIX THIS STUFF!!
    [Ud,Sd,Vd]=svd(D);
    Sd_dag=1./Sd(Sd>=1e-14);
    Sd_dag=diag(Sd_dag);
    D_Dag=Vd*Sd_dag*Ud';
    lambda=1*10^0;

    %AugMat=[lambda*A,zeros(NG,2);A*D_Dag',Eh1f',Eh2f'];

    %AugVec=[lambda*ef';cos(xx)']; 
    
    AugMat=A;
    AugVec=ef';
    
    [U,S,V]=svd(AugMat);
    S_Dag=zeros(size(S));
    S_Dag(S>0)=1./S(S>0);
    S_Dag=S_Dag';
    size(S_Dag)
    Res=V*(S_Dag*(U'*AugVec));
    %Cin=Res(1:end-2);
    %Cin=A\ef';
    %Cin=[A*D_Dag']\BB';
    Cin=Res;

    G=cos(xx')-real(A*D_Dag*Cin);
    H=[Eh1f',Eh2f'];
    Coef=H\G;
    
    Opt1=norm(real(A*Cin)-(1+Al)*cos(xx'))
    Opt2=norm(cos(xx')-real(A*D_Dag*Cin)-Coef(1)*Eh1f'-Coef(2)*Eh2f')
    figure(11)
    plot(abs(cos(xx')-real(A*D_Dag*Cin)-Coef(1)*Eh1f'-Coef(2)*Eh2f'))
end
