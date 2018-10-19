Al=logspace(-4.5,1,150);
S=zeros(150,8);
warning('off','all')
for j=0:7
    for i=1:150
        OP=RUNGraph(32,Al(i),9,j,5,3);
        S(i,j+1)=OP;
    end
end
