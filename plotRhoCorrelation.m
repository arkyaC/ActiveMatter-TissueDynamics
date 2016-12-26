rhoNorm=linspace(0,.6,11);
rhoNorm=rhoNorm(2:end);
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
%writing data to file
A=[rhoNorm;correl;err];
fileID = fopen('run 1.txt','w');
fprintf(fileID,'%10s %8s %8s\n','rhoNorm','order','error');
fprintf(fileID,'%6.5f %5.4f %6.5f\n',A);
fclose(fileID);
%plotting data
% errorbar(rhoNorm,correl,err,'b.');
% axis([0,1,0,1]);