Al=logspace(-4,1,30);
for i = 0:7
    for j=1:30
        [a(i+1,j),F1,~,d1(i+1,j),d2(i+1,j)]=RUNGraphOdd(256,Al(j),9,i,5,27,0);
        figure(i+1)
                hold on
        plot(F1)
        title(['f_c ',num2str(i),'th degree (even)'])

        [a2(i+1,j),F2,~,d12(i+1,j),d22(i+1,j)]=RUNGraphOdd(256,Al(j),9,i,5,27,1);
        figure(i+10)
                hold on
        plot(F2)
        title(['f_c ',num2str(i),'th degree (odd)'])

    end
    figure(i+20)
    loglog(Al,d1(i+1,:))
    title(['Accuracy in f even, deg=', num2str(i)])
    figure(i+30)
    loglog(Al,d12(i+1,:))
    title(['Accuracy in f odd, deg=',num2str(i)])
end


        
