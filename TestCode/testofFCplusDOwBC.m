% Test Code;
n=255;          m=27;           alpha=0.001;
%h=2/(n-1);      x=[-1:h:1]';    y=0*x+1;    yu=y;
h=2/(n-1);      x=[-1:h:1]';  

h1=exp((-x-1)/sqrt(alpha));
h2=exp((x-1)/sqrt(alpha));
y=x.^0;     yu=1+0*x;

Al=alpha;
for j=1:length(x)
GI(j)=(1/2)*Al*(-exp((2*(2*x(j)+1))/sqrt(Al))+exp((2*(x(j)+2))/sqrt(Al))+2*exp((1+3*x(j))/sqrt(Al))-2*exp((x(j)+3)/sqrt(Al))-exp(2*x(j)/sqrt(Al))+exp(2/sqrt(Al)))*exp(-2*x(j)/sqrt(Al))/(exp(4/sqrt(Al))-1)+(1/2)*Al*(exp((2*(2*x(j)+1))/sqrt(Al))-2*exp((3*(x(j)+1))/sqrt(Al))+exp((2*(x(j)+2))/sqrt(Al))-exp(2*x(j)/sqrt(Al))+2*exp((x(j)+1)/sqrt(Al))-exp(2/sqrt(Al)))*exp(-2*x(j)/sqrt(Al))/(exp(4/sqrt(Al))-1);
end

yu=yu-yu(1)*h1-yu(end)*h2;
figure(10)

plot(x,abs(GI'/Al-yu))
yu=GI'/Al;


%[zeros(255,2);[h1,h2]]
[AF,AD,DO]=FCplusDO(n,m,alpha);    AI=[[AF;AD],[zeros(255,2);[h1,h2]]];
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
    
f=AF*fc(1:end-2);
u=AD*fc(1:end-2)+fc(end-1)*h1+fc(end)*h2;

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
