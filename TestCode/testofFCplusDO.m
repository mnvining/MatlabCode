% Test Code;
n=255;          m=27;           alpha=0.00000001;
%h=2/(n-1);      x=[-1:h:1]';    y=0*x+1;    yu=y;
h=2/(n-1);      x=[-1:h:1]';    y=x.^5;     yu=x.^5+20*alpha*x.^3+120*alpha*alpha*x;
[AF,AD,DO]=FCplusDO(n,m,alpha);    AI=[AF;AD];
RHS=[y;yu];
partialsystem=0;
if partialsystem
    [U,S,V]=svd(AF);
    fc=DO*V*(pinv(S,1e-12)*(U'*(yu)));
    size(fc)
else
    [U,S,V]=svd(AI);
    fc=V*(pinv(S,1e-12)*(U'*(RHS)));
end
    
f=AF*fc;
u=AD*fc;

figure(1)
plot(real(f))
hold on;
plot(real(u))
hold off;

figure(5)
plot(real(f-y))
hold on;
plot(real(u-yu))
hold off;

figure(2)

G=diag(S'); GC=abs(U'*RHS);
loglog(G(1:m),GC(1:m),'o')
figure(3)
loglog(real(GC(1:m)./G(1:m)),'o');
max(abs(f-y))
max(abs(u-yu))
