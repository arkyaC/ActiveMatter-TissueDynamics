rho=linspace(0,3,15);
N=10;
correl=zeros(1,length(rho));
err=zeros(1,length(rho));
for k=1:length(rho)
    trials=zeros(1,N);
    for i=1:N
        trials(i)=vicsek(rho(k));
    end
    correl(k)=sum(trials)/N;
    err(k)=sqrt(sum((trials-correl(k)).^2)/N);
end
errorbar(rho,correl,err,'b.');
axis([0,1,0,1]);