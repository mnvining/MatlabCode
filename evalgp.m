function y=evalgp(xx,d,n)
% function y=evalgp(x,d,n)
% evaluate gram poly of deg d created from n-1 pts at point x. 
syms x
f = fmaker(d,n);
y=double(subs(f,xx));

