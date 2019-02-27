n=11;
h=2/(n-1);
FP=2+26*h;
x=[-1:h:1];
x=x';
y=1./(1+25*x.^2);
f=fcgram2(x,y,9);
f=[y(1:end-10);f];
FT=fft(f);
FT=[FT(1:18);FT(19)/2;zeros(9*length(FT)-1,1);FT(19)/2;FT(20:end)];
IFT=real(ifft(10*FT));
xx=[0:h/10:FP-h/10]-1;
xx=xx';
yy=1./(1+25*xx.^2);
figure(2)
EV=20/h+1;
plot(xx(1:EV),yy(1:EV)-IFT(1:EV))
figure(3)
plot(xx(1:EV),yy(1:EV),xx(1:EV),IFT(1:EV))

% focus on error for first 91 values