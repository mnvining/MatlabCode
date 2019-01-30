function f=fcgram2(x,y,d)
% function f=fcgram(x,y)
% input:    x, on [-1,1]
%           y, on [-1,1]
x1=x(1:10);
x2=x(end-9:end);
% First 8 y-values
y1=y(1:10);
% Last 8 y-values
y2=y(end-9:end);

% on the right side
C1=ptstocoeffs2(x2,y2,d);
% on the left side
C2=ptstocoeffs2(x1,y1,d);

% loads the even coefficients for degrees 0:8
load('EvenCo10.mat')
% loads the odd coefficients for degrees 0:8
load('OddCo10.mat')



LSE=EvenCo10(1:end-25,1:d+1)*C1;
LSO=OddCo10(1:end-25,1:d+1)*C1;
RSE=EvenCo10(1:end-25,1:d+1)*C2;
RSO=OddCo10(1:end-25,1:d+1)*C2;

f=RSE/2-RSO/2+(LSE/2+LSO/2);
f=f(1:end-10);

