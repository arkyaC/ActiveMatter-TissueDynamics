#incomplete
import numpy as np
import matplotlib

L=14;
driftspeed = 1;
mu=1;
tau=1;
Req=5.0/6;R0=1;
Fadh=0.75;Frep=30;
x =np.random.uniform(0,L,10000);
y =np.random.uniform(0,L,10000);
minAllowableDistance = 0.001;
numberOfPoints = floor(2*rhoNorm* L**2 /(3.14* Req**2));
noise=.6; #stochastic noise parameter

# Initialize first point.
keeperX = keeperY = np.zeros(1)
keeperY[0] = y[0];
keeperX[0] = x[0];
# Try dropping down more points.
counter = 1;
k=1;
while counter<=numberOfPoints
	# Get a trial point.
	thisX = np.multiply(x[k],np.ones_like(keeperX));
	thisY = np.multiply(y[k],np.ones_like(keeperY));
	# See how far is is away from existing keeper points.
	distances = sqrt((thisX - keeperX)**2 + (thisY - keeperY)**2); #need to change this!!
	minDistance = min(distances);
    if minDistance >= minAllowableDistance
        keeperX.append(thisX);
        keeperY.append(thisY);
        counter = counter + 1;

    k=k+1;

# plot(keeperX, keeperY, 'b.');
# axis([0,L,0,L]);
# axis square; #make a 1:1 plot
# grid on;
# pause(0.2);

numberOfPoints = length(keeperX);
Nsteps=4000;cutoffIter=1500;
theta = 2*pi*rand(1,numberOfPoints);
timedelta=0.05*R0/driftspeed;

y=zeros(Nsteps+1,3*numberOfPoints);
y(1,:) = [keeperX,keeperY,theta];
correl=0;
for k=1:Nsteps
	posX=y(k,1:numberOfPoints); # x position matrix
	posY=y(k,numberOfPoints+1:2*numberOfPoints); # y position matrix
	theta=y(k,2*numberOfPoints+1:);
	drift=[driftspeed*cos(theta) driftspeed*sin(theta)];
	for i=1:length(posX)
		interaxn=[0,0];
		for j=1:length(posX)
			if j==i
				continue;

			dij=sqrt((posX(i)-posX(j)).^2 + (posY(i) - posY(j)).^2);
			eij=[posX(i)-posX(j),posY(i)-posY(j)]/dij;
			if dij<Req
			    interaxn=interaxn-Frep*(dij-Req)/Req * eij;
			elseif dij>=Req && dij<=R0
			    interaxn=interaxn-Fadh*(dij-Req)/(R0-Req) * eij;


		interaxn=mu*interaxn;
		rhs(i)=drift(i)+interaxn(1);
		rhs(numberOfPoints+i)=drift(numberOfPoints+i)+interaxn(2);

    s=[0,0];
	for i=1:length(posX)
		ni=[drift(i) drift(numberOfPoints+i) 0]/driftspeed;
		vi=[rhs(i) rhs(numberOfPoints+i) 0];
		normvi=sqrt(vi(1)^2+vi(2)^2);
		vi=vi/normvi;
        s=s+[vi(1) vi(2)];
		ez=[0,0,1];
		res(i)=asin(dot(ez,cross(ni,vi)));

    if k>cutoffIter
        correl=correl+sqrt(s(1)^2+s(2)^2)/numberOfPoints; #correlation factor

	res=(1/tau)*res;
	posX=posX+timedelta*rhs(1:numberOfPoints);
    posY=posY+timedelta*rhs(numberOfPoints+1:);
	posX=L*(posX/L-floor(posX/L));
	posY=L*(posY/L-floor(posY/L));
	theta=theta+timedelta*res+sqrt(timedelta*(noise^2)/12)*normrnd(0,1,[1 length(theta)]);
	y(k+1,:)=[posX,posY,theta];#Euler stepping

correl=correl/(Nsteps-cutoffIter);
