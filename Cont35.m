function [Res,acc,stab]=Cont35(points,NModes)
% function Res=Cont35(points)
% inputs:   points to be continued on a 3.5 periodic grid
%           NModes the number of modes desired to be used
% output: f_c for those points

N=length(points);
Z=N*3.5;
M=floor(NModes/2);
modes=[1:M,Z-M-1:Z];
A=createFWDcont(N,modes,5);
Cin=A\points;
f_c=real(A*Cin);
Res=f_c;
acc=abs(f_c-points);
stab=norm(f_c)/norm(points);
end

