rho = 10.4;
N = 10;
noiseMax = 2*pi;
noise = linspace(0,noiseMax,18);
correl=zeros(1,length(noise));
err=zeros(1,length(noise));
for k=1:length(noise)
    trials=zeros(1,N);
    for i=1:N
        trials(i)=vicsek(rho,noise(k));
    end
    correl(k)=sum(trials)/N;
    err(k)=sqrt(sum((trials-correl(k)).^2)/N);
end
errorbar(noise,correl,err,'b.');
axis([0,noiseMax,0,1]);