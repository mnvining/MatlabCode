function f=fcgram(x,y)
% function f=fcgram(x,y)
% input:    x, on [-1,1]
%           y, on [-1,1]
x1=x(1:8);
x2=x(end-7:8);
% First 8 y-values
y1=y(1:8);
% Last 8 y-values
y2=y(end-7:end);

C1=ptstocoeffs(x2,y2,5);
C2=ptstocoeffs(x1,y1,5);

load('EvenCo.mat')
load('OddCo.mat')
load('ContMatr.mat')

LSE=EvenCo*C1;
LSO=OddCo*C1;
RSE=EvenCo*C2;
RSO=OddCo*C2;

Coeffs=LSE/2+LSO/2+RSE/2-RSO/2;

Cont=ContMatr*Coeffs;
f=real(Cont);

