function f=fcgram(x,y,d)
% function f=fcgram(x,y)
% input:    x, on [-1,1]
%           y, on [-1,1]
x1=x(1:8);
x2=x(end-7:end);
% First 8 y-values
y1=y(1:8);
% Last 8 y-values
y2=y(end-7:end);

% on the right side
C1=ptstocoeffs(x2,y2,d)
% on the left side
C2=ptstocoeffs(x1,y1,d)
C2(2)

% loads the even coefficients for degrees 0:8
load('EvenCo.mat')
% loads the odd coefficients for degrees 0:8
load('OddCo.mat')
% loads the Continuation Matrix "A" in other graphs
load('ContMatr.mat')


LSE=EvenCo(1:end-20,1:d+1)*C1;
LSO=OddCo(1:end-20,1:d+1)*C1;
RSE=EvenCo(1:end-20,1:d+1)*C2;
RSO=OddCo(1:end-20,1:d+1)*C2;

f=RSE/2-RSO/2+(LSE/2+LSO/2);
f=f(1:end-8);

