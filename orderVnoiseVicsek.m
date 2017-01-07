rho = 40/(3.1^2);
N = 10;
noiseMax = 5;
noise = linspace(0,noiseMax,30);
correl=zeros(1,length(noise));
err=zeros(1,length(noise));
for k=1:length(noise)
    trials=zeros(1,N);
    for i=1:N
        orderParamData=vicsek(rho,noise(k));
        trials(i) = sum(orderParamData(end-500:end))/500;
    end
    correl(k)=sum(trials)/N;
    err(k)=sqrt(sum((trials-correl(k)).^2)/N);
end
% writing data
dataToWrite = [noise;correl;err];
fileID = fopen('dump.txt','w');
fprintf(fileID,'%8s \t %8s \t %8s\n','noise','order','error');
fprintf(fileID,'%5.4f \t %5.4f \t %6.5f\n',dataToWrite);
fclose(fileID);
% plotting graph
errorbar(noise,correl,err,'b.');
axis([0,noiseMax,0,1]);
