% Test Code;
n=255;          m=27;           alpha=0.001;
%h=2/(n-1);      x=[-1:h:1]';    y=0*x+1;    yu=y;
h=2/(n-1);      x=[-1:h:1]';  

h1=exp((-x-1)/sqrt(alpha));
h2=exp((x-1)/sqrt(alpha));
z=[1,0,0,0];
Al=alpha;
y=x.^3;     yu=z(1)*x.^3+z(2)*x.^2+(z(3)+6*z(1)*Al)*x+(z(4)+2*z(2)*Al);


for j=1:length(x)
GI(j)=6*Al^2*(exp(1/sqrt(Al)))^3/(((exp(1/sqrt(Al)))^4-1)*exp(x(j)/sqrt(Al)))-6*exp(x(j)/sqrt(Al))*Al^2*exp(1/sqrt(Al))/((exp(1/sqrt(Al)))^4-1)+Al*(exp(1/sqrt(Al)))^3/(((exp(1/sqrt(Al)))^4-1)*exp(x(j)/sqrt(Al)))-exp(x(j)/sqrt(Al))*Al*exp(1/sqrt(Al))/((exp(1/sqrt(Al)))^4-1)-6*exp(x(j)/sqrt(Al))*Al^2*(exp(1/sqrt(Al)))^3/((exp(1/sqrt(Al)))^4-1)-exp(x(j)/sqrt(Al))*Al*(exp(1/sqrt(Al)))^3/((exp(1/sqrt(Al)))^4-1)+6*Al^2*exp(1/sqrt(Al))/(((exp(1/sqrt(Al)))^4-1)*exp(x(j)/sqrt(Al)))+Al*exp(1/sqrt(Al))/(((exp(1/sqrt(Al)))^4-1)*exp(x(j)/sqrt(Al)))-Al*x(j)^3/((exp(1/sqrt(Al)))^4-1)-6*Al^2*x(j)/((exp(1/sqrt(Al)))^4-1)+(exp(1/sqrt(Al)))^4*Al*x(j)^3/((exp(1/sqrt(Al)))^4-1)+6*(exp(1/sqrt(Al)))^4*Al^2*x(j)/((exp(1/sqrt(Al)))^4-1);
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
