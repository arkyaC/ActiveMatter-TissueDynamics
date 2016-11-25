rhoNorm=linspace(0,.6,20);
N=10;
correl=zeros(1,length(rhoNorm));
err=zeros(1,length(rhoNorm));
for k=1:length(rhoNorm)
    trials=zeros(1,N);
    for i=1:N
        trials(i)=noisyNucl(rhoNorm(k));
    end
    correl(k)=sum(trials)/N;
    err(k)=sqrt(sum((trials-correl(k)).^2)/N);
end
errorbar(rhoNorm,correl,err,'b.');
axis([0,1,0,1]);