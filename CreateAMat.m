function [A,Modes,ZZ,M]=CreateAMat(NoMo,cl,NG)
M=2*floor(NoMo/2);% this is x2 because we need double modes

ZZ=2*(cl*NG/2+NG); % this is x2 because of repeated dom.
A=dftmtx(ZZ);
A=A';
A=A/sqrt(ZZ);
Modes=[1:(M+1),(ZZ-M+1):ZZ];
% cont matrix to eval at all
A=A(:,Modes);