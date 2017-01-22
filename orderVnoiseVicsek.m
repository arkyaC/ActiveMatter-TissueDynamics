rho = 4;
N = 5;
noiseMax = 5.5;
noise = linspace(0,noiseMax,20);
correl=zeros(1,length(noise));
err=zeros(1,length(noise));
for k=1:length(noise)
    display(k); %to check progress
    trials=zeros(1,N);
    for i=1:N
        orderParamData=vicsek(rho,noise(k));
        trials(i) = mean(orderParamData(end-200:end));
    end
    correl(k)=mean(trials);
    err(k)=sqrt(mean((trials-correl(k)).^2));
end
% writing data
dataToWrite = [noise;correl;err];
fileID = fopen('data/dump.txt','w');
fprintf(fileID,'%8s \t %8s \t %8s\n','noise','order','error');
fprintf(fileID,'%5.4f \t %5.4f \t %6.5f\n',dataToWrite);
fclose(fileID);
% plotting graph
errorbar(noise,correl,err,'b.');
axis([0,noiseMax,0,1]);
xlabel('Noise');ylabel('Order Parameter');
