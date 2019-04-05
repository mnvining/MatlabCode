function c=HTest(d,Al)

clear A
x=[-1:2/9:1]';
xe=[-1:.001:1]';

for i=1:d+1
    A(:,i)=x.^(i-1);
end

[~,R]=qr(A,0);
Rinv=inv(R);
y=zeros(size(x));
ye=zeros(size(xe));
Sol=y;
Sole=ye;
for i=1:d+1
    co(i)=Rinv(i,d+1);
    y=y+Rinv(i,d+1)*x.^(i-1);
    ye=ye+Rinv(i,d+1)*xe.^(i-1);
    s=PSols(i-1,Al);
    Sol=Sol+Rinv(i,d+1)*s(x);
    Sole=Sole+Rinv(i,d+1)*s(xe);
end
h1=exp((-1-x)/sqrt(Al));
h2=exp(-(-x+1)/sqrt(Al));


f=@(c) norm(Sol-(-1)^d*c*h1+-c*h2)-1;
c=fzero(f,Al);




