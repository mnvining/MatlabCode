for i=1:50
c(i)=HTest(5,Al(i));
end
loglog(Al,abs(c))
H=[Al',c'];
semilogx(Al,H(:,2)./(H(:,1).^2))
c'
H(:,2)./(H(:,1).^2)
