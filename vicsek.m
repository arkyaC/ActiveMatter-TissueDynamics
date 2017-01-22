function orderN = vicsek(rho,noise)
numberOfPoints = 100;
L = sqrt(numberOfPoints/rho);
v = 0.03; %as stated in the paper (for optimum results)
r = 1; %definition of neighbourhood for averaging

pointsX =L* rand(1, numberOfPoints);
pointsY =L* rand(1, numberOfPoints);
Nsteps=1000;
theta = 2*pi*rand(1,numberOfPoints) - pi; %-pi to pi
%theta = pi/2 * ones(1,numberOfPoints);
timedelta = 1; %as mentioned in paper
y = zeros(Nsteps+1,3*numberOfPoints); % initialisation
y(1,:) = [pointsX,pointsY,theta];
orderN = zeros(1,Nsteps);
for k=1:Nsteps
	posX=y(k,1:numberOfPoints); % x position matrix
	posY=y(k,numberOfPoints+1:2*numberOfPoints); % y position matrix
	theta=y(k,2*numberOfPoints+1:end); % theta matrix
    
    newTheta=theta; %just arbit initialisation, supposed to be faster this way
	vel=[v*cos(theta) v*sin(theta)];
    s = [0,0]; %avg vel init
    for i=1:numberOfPoints %evaluating updated angle parameter
        avgcos = 0;
        avgsin = 0;
        ctr = 0;
        ni = [vel(i) vel(numberOfPoints+i)]/v;
        s = s + ni;
        for j=1:numberOfPoints
        	dx = posX(i) - posX(j);
			dy = posY(i) - posY(j);
            if abs(dx) > 0.5*L
				dx = dx - L*sign(dx);
            end
            if abs(dy) > 0.5*L
				dy = dy - L*sign(dy);
            end
            if (dx^2+dy^2)<(r^2)
                avgcos = avgcos + cos(theta(j));
                avgsin = avgsin + sin(theta(j));
                ctr = ctr + 1;
            end
        end
        avgcos = avgcos/ctr;
        avgsin = avgsin/ctr;
        angleMean = atan2(avgsin,avgcos);
        temp = angleMean + (noise*rand - noise/2);
        %keep angle between -pi and pi
        newTheta(i) = temp-2*pi*floor((temp+pi)/(2*pi));
    end
    orderN(k)=sqrt(s(1)^2+s(2)^2)/numberOfPoints;
    theta = newTheta; % update angle values
    posX = posX + vel(1:numberOfPoints)*timedelta;
    posY = posY + vel(numberOfPoints+1:end)*timedelta;
    posX = posX-L*floor(posX/L); % periodic BC
	posY = posY-L*floor(posY/L);
    y(k+1,:)=[posX,posY,theta];
end
% histogram(theta,100);
% axis([-pi,pi,0,100]);

% write data to dump
% timeSteps = 1:Nsteps;
% fileID = fopen('data/dump.txt','w');
% fprintf(fileID,'%d \t %6.5f \n',[timeSteps;orderN]);
% fclose(fileID);

%plot order parameter against time
% plot(linspace(0,Nsteps,Nsteps),orderN);
% axis([0,Nsteps,0,1]);
% xlabel('Time step');ylabel('Order Parameter');

%movie
% pause(10);
% for k=1:size(y)
% 	posX=y(k,1:numberOfPoints); % x position matrix
% 	posY=y(k,numberOfPoints+1:2*numberOfPoints); % y position matrix
% 	plot(posX,posY,'b.'); %plot instantaneous position
%     axis([0,L,0,L]);
%     axis square;
%     grid on;
% 	pause(.1);
% end