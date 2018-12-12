NG=256; Al=1; n=9; cl=5; NoMo=27; 
EvenCo=zeros(53,8);
OddCo=zeros(53,8);
for i = 0:7
    [F1,Cin]=RUNGraphOdd(NG,Al,n,i,cl,NoMo,1);
    [F2,Cin2]=RUNGraphOdd(NG,Al,n,i,cl,NoMo,0);
    EvenCo(:,i+1)=Cin;
    OddCo(:,i+1)=Cin2;
end
