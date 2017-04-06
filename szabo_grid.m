%max allowed density is 1.39
% rho = (pi * Req^2 * N) / L^2
% delX^2 * N = L^2
% rhoMax = 1.39 when delX = 1.5*Req

function orderN=szabo_grid(rho,noise)
%rhoNorm = 0.1;
%L=14;
numberOfPoints = 100;N=numberOfPoints; % must be perfect square
v0 = 1;
mu = 1;
tau = 1;
Req = 5/6;R0 = 1;
Fadh = 0.75;Frep = 30;%Fwall = 50;
% Fadh = 0.75/6;Frep = 30/6;%Fwall = 50/6;

L = Req * sqrt(pi*numberOfPoints/rho);
delX = L/sqrt(N);

x = delX/2:delX:delX*(sqrt(N)-1)+delX/2;
keeperX = repmat(x,1,sqrt(N));
y = x;
tempY = repmat(y,sqrt(N),1);
keeperY = reshape(tempY,1,N);

Nsteps=20000; %number of timeSteps
% cutoffIter=Nsteps-100;
theta = 2*pi*rand(1,numberOfPoints) - pi; %need it in the range of -pi to pi
%timedelta=0.05*R0/v0;
timedelta=0.005;

y=zeros(Nsteps+1,3*numberOfPoints);
y(1,:) = [keeperX,keeperY,theta];
%vel=zeros(1,2*numberOfPoints);
orderN = zeros(1,Nsteps);
for k=1:Nsteps
	posX=y(k,1:numberOfPoints); % x position matrix
	posY=y(k,numberOfPoints+1:2*numberOfPoints); % y position matrix
	theta=y(k,2*numberOfPoints+1:end);
    
	drift=[v0*cos(theta) v0*sin(theta)];
    rhs = zeros(1,2*numberOfPoints);
    for i=1:numberOfPoints
        interaxn=[0,0];
        for j=1:numberOfPoints
            if j==i
                continue;
            end
            dx = posX(i) - posX(j);
			dy = posY(i) - posY(j);
            if abs(dx) > 0.5*L
				dx = dx - L*sign(dx);
            end
            if abs(dy) > 0.5*L
				dy = dy - L*sign(dy);
            end
            dijSq = dx^2 + dy^2;
            dij = sqrt(dijSq);
            eij = [dx,dy]/dij;
%             dij=sqrt((posX(i)-posX(j)).^2 + (posY(i) - posY(j)).^2);
%             eij=[posX(i)-posX(j),posY(i)-posY(j)]/dij;
            if dijSq<Req^2
                interaxn = interaxn - Frep*(dij-Req)/Req * eij;
            elseif dijSq>=Req^2 && dijSq<=R0^2
                interaxn = interaxn - Fadh*(dij-Req)/(R0-Req) * eij;
            end
        end
        %for wall force
    %         Fx = 0;Fy = 0;
    %         if posX(i)<R0
    %             Fx = Fx + Fwall*exp(-2*posX(i)/R0);
    %         end
    %         if L-posX(i)<R0
    %             Fx = Fx - Fwall*exp(-2*(L-posX(i))/R0);
    %         end
    %         if posY(i)<R0
    %             Fy = Fy + Fwall*exp(-2*posY(i)/R0);
    %         end
    %         if L-posY(i)<R0
    %             Fy = Fy - Fwall*exp(-2*(L-posY(i))/R0);
    %         end
    %         interaxn = interaxn + [Fx Fy];
        interaxn=mu*interaxn;
        rhs(i)=drift(i)+interaxn(1);
        rhs(numberOfPoints+i)=drift(numberOfPoints+i)+interaxn(2);
    end
    s=[0,0];
    res = zeros(1,numberOfPoints);
    for i=1:numberOfPoints
        ni=[drift(i) drift(numberOfPoints+i)]/v0;
        vi=[rhs(i) rhs(numberOfPoints+i)];
        normvi=sqrt(vi(1)^2+vi(2)^2);
        vi=vi/normvi;
        s=s+[vi(1) vi(2)];
        res(i)=asin(ni(1)*vi(2)-ni(2)*vi(1));
    end
    orderN(k)=sqrt(s(1)^2+s(2)^2)/numberOfPoints;
	res = (1/tau)*res;
	posX = posX + timedelta*rhs(1:numberOfPoints);
    posY = posY + timedelta*rhs(numberOfPoints+1:end);
	posX = posX - L * floor(posX/L);
	posY = posY - L * floor(posY/L);
% 	theta = theta + timedelta*res + sqrt(timedelta*(noise^2)/12)*normrnd(0,1,[1 length(theta)]);
    noiseAmp = noise/sqrt(timedelta);
    theta = theta + timedelta * (res + noiseAmp * (rand(1,numberOfPoints) - 0.5));
	y(k+1,:) = [posX,posY,theta];%Euler stepping
    if(mod(k,200) == 0) 
        k %display status
    end
end

% write data to dump
% timeSteps = 1:Nsteps;
% fileID = fopen('data/dump.txt','w');
% fprintf(fileID,'%d \t %6.5f \n',[timeSteps;orderN]);
% fclose(fileID);

% %movie
% for k=1:size(y)
% 	posX=y(k,1:numberOfPoints); % x position matrix
% 	posY=y(k,numberOfPoints+1:2*numberOfPoints); % y position matrix
% 	plot(posX,posY,'b.'); %plot instantaneous position
%     axis([0,L,0,L]);
%     axis square;
%     grid on;
%     %frames(k) = getframe(gcf);
%  	pause(.0005);
% end

%plotting the order parameter value against step number
figure
plot(linspace(0,Nsteps,Nsteps),orderN);
axis([0,Nsteps,0,1]);
xlabel('Time step');ylabel('Order Parameter');
% drawnow

%writing to a file
% A=[rho;correl;err];
% fileID = fopen('run 1.txt','w');
% fprintf(fileID,'%10s %8s %8s\n','rho','order','error');
% fprintf(fileID,'%6.5f %5.4f %6.5f\n',A);
% fclose(fileID);
