% rhoNorm=linspace(0,.6,11);
% rhoNorm=rhoNorm(2:end);
rhoNorm = 0.6;
noise = linspace(0,1.4,20);
noise = noise(2:end);
N=10;
correl=zeros(1,length(noise));
err=zeros(1,length(noise));
for k=1:length(noise)
    trials=zeros(1,N);
    for i=1:N
        orderN=noisyNucl(rhoNorm,noise(k));
        trials(i) = mean(orderN(end-100:end));
    end
    correl(k)=mean(trials);
    err(k)=sqrt(mean( (trials-correl(k)).^2 ));
end
%writing data to file
A=[noise;correl;err];
fileID = fopen('data/dump.txt','w');
fprintf(fileID,'%10s %8s %8s\n','noise','order','error');
fprintf(fileID,'%6.5f %5.4f %6.5f\n',A);
fclose(fileID);
%plotting data
errorbar(rhoNorm,correl,err,'b.');
axis([0,1,0,1]);