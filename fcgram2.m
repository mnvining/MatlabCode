function f=fcgram2(x,y,d)
% function f=fcgram(x,y)
% input:    x, on [-1,1]
%           y, on [-1,1]
x1=x(10:-1:1);
x2=x(end-9:end);
% First 10 y-values
y1=y(10:-1:1);
% Last 10 y-values
y2=y(end-9:end);

% on the right side
C1=ptstocoeffs2(x2,y2,d);
% on the left side
C2=ptstocoeffs2(x1,y1,d);

% loads the even coefficients for degrees 0:8
load('EvenCo10Stab.mat')
% loads the odd coefficients for degrees 0:8
load('OddCo10Stab.mat')


CEven=(C1+C2)/2
COdd=(C1-C2)/2

EvenC=EN3FCEven(1:end-25,1:d+1)*CEven
figure(10)
plot(EvenC)
OddC=EN3FCOdd(1:end-25,1:d+1)*COdd;


f=EvenC+OddC;
f=f(1:end-10);

