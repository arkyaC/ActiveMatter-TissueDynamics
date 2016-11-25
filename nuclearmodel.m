clear all;close all;clc;
L=14;
driftspeed = 1;
mu=1;
tau=1;
Req=5/6;R0=1;
Fadh=0.75;Frep=30;
rhoNorm=0.12;
x =L* rand(1, 10000);
y =L* rand(1, 10000);
minAllowableDistance = Req;
numberOfPoints = floor(2*rhoNorm*L^2);
% Initialize first point.
keeperX = x(1);
keeperY = y(1);
% Try dropping down more points.
counter = 1;
k=1;
while counter<=numberOfPoints
	% Get a trial point.
	thisX = x(k);
	thisY = y(k);
	% See how far is is away from existing keeper points.
	distances = sqrt((thisX-keeperX).^2 + (thisY - keeperY).^2);
	minDistance = min(distances);
    if minDistance >= minAllowableDistance
        keeperX(counter) = thisX;
        keeperY(counter) = thisY;
        counter = counter + 1;
    end
    k=k+1;
end
% plot(keeperX, keeperY, 'b.');
% axis([0,L,0,L]);
% axis square; %make a 1:1 plot
% grid on;
% pause(0.2);

numberOfPoints = length(keeperX);
Nsteps=4000;
theta = 2*pi*rand(1,numberOfPoints);
timedelta=0.05*R0/driftspeed;

y=zeros(Nsteps+1,3*numberOfPoints);
y(1,:) = [keeperX,keeperY,theta];
correl=0;
for k=1:Nsteps
	posX=y(k,1:numberOfPoints); % x position matrix
	posY=y(k,numberOfPoints+1:2*numberOfPoints); % y position matrix
	theta=y(k,2*numberOfPoints+1:end);
	drift=[driftspeed*cos(theta) driftspeed*sin(theta)];
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
		interaxn=mu*interaxn;
		rhs(i)=drift(i)+interaxn(1);
		rhs(numberOfPoints+i)=drift(numberOfPoints+i)+interaxn(2);
    end
    s=[0,0];
	for i=1:length(posX)
		ni=[drift(i) drift(numberOfPoints+i) 0]/driftspeed;
		vi=[rhs(i) rhs(numberOfPoints+i) 0];
		normvi=sqrt(vi(1)^2+vi(2)^2);
		vi=vi/normvi;
        s=s+[vi(1) vi(2)];
		ez=[0,0,1];
		res(i)=asin(dot(ez,cross(ni,vi)));
    end
    correl=correl+sqrt(s(1)^2+s(2)^2)/numberOfPoints; %correlation factor
	res=(1/tau)*res;
	posX=posX+timedelta*rhs(1:numberOfPoints);
    posY=posY+timedelta*rhs(numberOfPoints+1:end);
	posX=L*(posX/L-floor(posX/L));
	posY=L*(posY/L-floor(posY/L));
	theta=theta+timedelta*res;
	y(k+1,:)=[posX,posY,theta];%Euler stepping
end
correl=correl/Nsteps
for k=1500:size(y)
	posX=y(k,1:numberOfPoints); % x position matrix
	posY=y(k,numberOfPoints+1:2*numberOfPoints); % y position matrix
	plot(posX,posY,'b.'); %plot instantaneous position
    axis([0,L,0,L]);
    axis square;
    grid on;
	pause(.005);
end
