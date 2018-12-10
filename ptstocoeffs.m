function C=ptstocoeffs(x,y,d)
% function C=ptstocoeffs(x,y)
% calculates the coefficients of the gram polynomial expansion for the
% points (x,y).
% inputs: x - x-coords
%           y - corresponding y-coords [must be column vector]
% outputs: c - coefficients c_0,c_1... corresponding to the Gram Polys.

% Uses 0 - d.
p=zeros(length(x),d+1);
for i=1:length(x)
    for j=0:d
        p(i,j+1)=evalgp(x(i),j,9);
    end
end
size(p)
size(y)
C=p\y;

    