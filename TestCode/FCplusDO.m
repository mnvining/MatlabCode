function [AF,AD,DO]=FCplusDO(n,m,alpha)
%n is the number of points on the fine grid. point 1 is assumed to be -1,
%point n is assumed to be 1
%n should be odd.
%m is the number of modes in the Fourier Continuation
%m should be odd.
%-alpha u''+u = f  do the differential opperator is I-alpha D_xx
h=2/(n-1);
ml=floor(m/2);
xf=[-1:h:6-h]';
nf=length(xf)
AF=dftmtx(nf);
AF=AF'/sqrt(nf);
AF=AF(1:n,[1:ml+1,nf-ml+1:nf]);
DD=(2*pi/7)^2 * diag(-([[0:ml],[-ml:-1]]).^2);
DO=(-alpha*DD+eye(m,m));
DP=pinv(DO,1e-12);
AD=AF*DP;





%Code for testing
% [U,S,V]=svd(AF);
% cf=V*(pinv(S,1e-12)*(U'*(xf(1:n).^2)));
% size(cf)
% size(DD)
% size(AF)
% test=AF*(DD*cf);
% plot(xf(1:n),real(test))

end




