function correl=noisyNucl(rhoNorm)
%rhoNorm = 0.1;
%L=14;
driftspeed = 1;
mu = 1;
tau = 1;
Req = 5/6;R0 = 1;
Fadh = 0.75/6;Frep = 5;Fwall = 50/6;
minAllowableDistance = R0;
%numberOfPoints = floor(2*rhoNorm*L^2/3.14);
numberOfPoints = 1000;
L = sqrt(numberOfPoints*pi/(2*rhoNorm));
x = L* rand(1, 10000);
y = L* rand(1, 10000);
noise = 0.6; %stochastic noise parameter

% Initialize first point.
keeperX = x(1);
keeperY = y(1);
% Try dropping down more points.
counter = 1;
k=1;
while counter<numberOfPoints
	% Get a trial point.
	thisX = x(k);
	thisY = y(k);
	% See how far is is away from existing keeper points.
	distances = sqrt((thisX-keeperX).^2 + (thisY - keeperY).^2);
	minDistance = min(distances);
    if minDistance >= minAllowableDistance
        counter = counter + 1;
        keeperX(counter) = thisX;
        keeperY(counter) = thisY;
    end
    k=k+1;
end

numberOfPoints = length(keeperX);
Nsteps=6000;cutoffIter=Nsteps-100;
theta = 2*pi*rand(1,numberOfPoints) - pi; %need it in the range of -pi to pi
timedelta=0.05*R0/driftspeed;

y=zeros(Nsteps+1,3*numberOfPoints);
y(1,:) = [keeperX,keeperY,theta];
%vel=zeros(1,2*numberOfPoints);
correl=0;
orderN = zeros(Nsteps);
for k=1:Nsteps
	posX=y(k,1:numberOfPoints); % x position matrix
	posY=y(k,numberOfPoints+1:2*numberOfPoints); % y position matrix
	theta=y(k,2*numberOfPoints+1:end);
	drift=[driftspeed*cos(theta) driftspeed*sin(theta)];
    rhs = zeros(1,2*numberOfPoints);
	for i=1:length(posX)
		interaxn=[0,0];
		for j=1:length(posX)
			if j==i
				continue;
			end
			dij=sqrt((posX(i)-posX(j)).^2 + (posY(i) - posY(j)).^2);
			eij=[posX(i)-posX(j),posY(i)-posY(j)]/dij;
			if dij<Req
			    interaxn=interaxn-Frep*(dij-Req)/Req * eij;
			elseif dij>=Req && dij<=R0
			    interaxn=interaxn-Fadh*(dij-Req)/(R0-Req) * eij;
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
    res = zeros(1,length(theta));
	for i=1:length(theta)
		ni=[drift(i) drift(numberOfPoints+i)]/driftspeed;
		vi=[rhs(i) rhs(numberOfPoints+i)];
		normvi=sqrt(vi(1)^2+vi(2)^2);
		vi=vi/normvi;
        s=s+[vi(1) vi(2)];
        res(i)=asin(ni(1)*vi(2)-ni(2)*vi(1));
		%res(i)=asin(dot(ez,cross(ni,vi)));
    end
    if k>cutoffIter
        correl=correl+sqrt(s(1)^2+s(2)^2)/numberOfPoints; %correlation factor
    end
    orderN(k)=sqrt(s(1)^2+s(2)^2)/numberOfPoints;
	res=(1/tau)*res;
	posX=posX+timedelta*rhs(1:numberOfPoints);
    posY=posY+timedelta*rhs(numberOfPoints+1:end);
	posX=L*(posX/L-floor(posX/L));
	posY=L*(posY/L-floor(posY/L));
	theta=theta+timedelta*res+sqrt(timedelta*(noise^2)/12)*normrnd(0,1,[1 length(theta)]);
	y(k+1,:)=[posX,posY,theta];%Euler stepping
end
%plotting the order parameter value against step number
plot(linspace(0,Nsteps,Nsteps),orderN);
axis([0,Nsteps,0,1]);
xlabel('Time step');ylabel('Order Parameter');
drawnow
%plotting the data
% for k=1:size(y)
% 	posX=y(k,1:numberOfPoints); % x position matrix
% 	posY=y(k,numberOfPoints+1:2*numberOfPoints); % y position matrix
% 	plot(posX,posY,'b.'); %plot instantaneous position
%     axis([0,L,0,L]);
%     axis square;
%     grid on;
%     %frames(k) = getframe(gcf);
%  	pause(.001);
% end

%writing to a file
% A=[rhoNorm;correl;err];
% fileID = fopen('run 1.txt','w');
% fprintf(fileID,'%10s %8s %8s\n','rhoNorm','order','error');
% fprintf(fileID,'%6.5f %5.4f %6.5f\n',A);
% fclose(fileID);
correl=correl/(Nsteps-cutoffIter);
