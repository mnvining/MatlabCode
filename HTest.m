d=input('What is d?: ')
clear A
x=[-1:2/9:1]';
xe=[-1:.001:1]';

for i=1:d+1
    A(:,i)=x.^(i-1);
end

[Q,R]=qr(A,0);
Rinv=inv(R);
y=zeros(size(x));
ye=zeros(size(xe));
for i=1:d+1
    co(i)=Rinv(i,d+1);
    y=y+Rinv(i,d+1)*x.^(i-1);
    ye=ye+Rinv(i,d+1)*xe.^(i-1);
end

h=zeros(size(x));
h(1)=1;h(end)=1;



