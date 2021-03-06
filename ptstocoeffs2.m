function C=ptstocoeffs2(x,y,d)
% function C=ptstocoeffs(x,y)
% calculates the coefficients of the gram polynomial expansion for the
% points (x,y).
% inputs: x - x-coords
%           y - corresponding y-coords [must be column vector]
% outputs: c - coefficients c_0,c_1... corresponding to the Gram Polys.

% Uses 0 - d.


% loads the even coefficients for degrees 0:8
load('EvenCo10Stab.mat')
p=EN3FCEven(1:10,1:d+1);
C=p\y

    