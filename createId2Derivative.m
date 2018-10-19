function L2=createId2Derivative(N,M,Alpha,cl)
% function L2=createId2Derivative(N,M,Alpha,cl)
% Inputs:   N, fine grid length [on original domain]
%           M, half the # modes desired
%           Alpha - alpha value for operator
%           cl, desired continuation length [length of interval past [-1,1)
%               for 3.5 times periodicity, cl=5.


fd=(cl+2)/2;
ZZ=round(cl*N/2+N);
LL=round(ZZ/2);

K=-pi*1i*[0:LL-1,0,-(LL-1):-1]';

% derivative coefficients pertaining to the FC
n=[K(1:M);K((ZZ)-M+1:(ZZ))]/fd;


L2=eye(length(n))-Alpha*diag(n.*n);