function [Res,acc,stab]=Cont35Odd(points,NModes)
% function Res=Cont35Odd(points)
% inputs:   points to be continued on a 3.5 periodic grid
%           NModes the number of modes desired to be used
% output:   Res: f_c for those points [full, w/ continuation]
%           acc: difference between f and f_c [both f's]
%           stab: norm(f_c)/norm(f)

N=length(points);
Z=2*N*3.5;
M=floor(NModes/2);
modes=[1:M,Z-M-1:Z];

% create continuation matrix
A=dftmtx(Z);
A=A';
A=A(:,modes);
% pull out specific point sections
A1=A(1:N,:);
A2=A(Z/2+1:Z/2+N,:);
% create matrix and vector [odd] for solve
F=[points;-points];
AMat=[A1;A2];
size(AMat)
size(F)
Cin=AMat\F;
% full continuation
f_c=real(A*Cin);
% section 1
f_c1=real(A1*Cin);
% section 2
f_c2=real(A2*Cin);
Res=f_c;
acc=abs([f_c1;f_c2]-F);
stab=norm(f_c1)/norm(points);
end