function orderN = vicsek(rho,noise)
numberOfPoints = 40;
%rho = 10.4;
L = sqrt(numberOfPoints/rho);
%noise = .4;
v = 0.03; %as stated in the paper (for optimum results)
r = 1; %definition of neighbourhood for averaging
x =L* rand(1, 1000000);
y =L* rand(1, 1000000);
minAllowableDistance = 0; %arbitrary (just for better visualisation) %CHANGE to 0 for testing
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
        keeperX(counter) = thisX;
        keeperY(counter) = thisY;
        counter = counter + 1;
    end
    k=k+1;
end

Nsteps=10000;cutOffIter=Nsteps-200;
theta = 2*pi*rand(1,numberOfPoints) - pi; %-pi to pi
%theta = pi/2 * ones(1,numberOfPoints);
timedelta = 1; %as mentioned in paper
y = zeros(Nsteps+1,3*numberOfPoints);
y(1,:) = [keeperX,keeperY,theta];
correl = 0; %to be evaluated later
counter = 0;
orderN = zeros(1,Nsteps);
for k=1:Nsteps
	posX=y(k,1:numberOfPoints); % x position matrix
	posY=y(k,numberOfPoints+1:2*numberOfPoints); % y position matrix
	theta=y(k,2*numberOfPoints+1:end);
    newTheta=theta; %just arbit initialisation
	vel=[v*cos(theta) v*sin(theta)];
    s = [0,0]; %avg vel init
    for i=1:numberOfPoints %evaluating updated angle parameter
        avg = 0;
        ctr = 0;
        ni = [vel(i) vel(numberOfPoints+i)]/v;
        s = s + ni;
        for j=1:numberOfPoints
            distIJ = sqrt((posX(j)-posX(i))^2 + (posY(j) - posY(i))^2);
            if distIJ<r
                avg = avg + theta(j);
                ctr = ctr + 1;
            end
        end
        avg = avg/ctr;
        temp = avg + (noise*rand - noise/2);
        %keep angle between -pi and pi
        newTheta(i) = temp - 2*pi*floor( (temp+pi)/(2*pi) );
        %newTheta(i) = wrapToPi(temp);
        %newTheta(i) = 2*pi * (temp/(2*pi) - floor(temp/(2*pi)));
    end
    
    orderN(k)=sqrt(s(1)^2+s(2)^2)/numberOfPoints;
    if k>cutOffIter
        correl = correl + orderN(k);
        counter = counter + 1;
    end
%     histogram(theta,100);
%     axis([-pi,pi,0,100]);
%     pause(0.05);
    theta = newTheta;
    %vel = [v*cos(theta) v*sin(theta)]; %updated vel
    posX = posX + vel(1:numberOfPoints)*timedelta;
    posY = posY + vel(numberOfPoints+1:end)*timedelta;
    posX = L*(posX/L-floor(posX/L));
	posY = L*(posY/L-floor(posY/L));
    y(k+1,:)=[posX,posY,theta];%time stepping
end
% histogram(theta,100);
% axis([-pi,pi,0,100]);

% write data to dump
% timeSteps = 1:Nsteps;
% fileID = fopen('dump.txt','w');
% fprintf(fileID,'%d \t %6.5f \n',[timeSteps;orderN]);
% fclose(fileID);

%plot order parameter against time
plot(linspace(0,Nsteps,Nsteps),orderN);
axis([0,Nsteps,0,1]);
xlabel('Time step');ylabel('Order Parameter');

% pause(10);

%movie
% plot2 = figure;
% for k=1:size(y)
% 	posX=y(k,1:numberOfPoints); % x position matrix
% 	posY=y(k,numberOfPoints+1:2*numberOfPoints); % y position matrix
% 	plot(posX,posY,'b.'); %plot instantaneous position
%     axis([0,L,0,L]);
%     axis square;
%     grid on;
% 	pause(.1);
% end