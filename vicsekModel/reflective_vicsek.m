%reflective boundaries

noise = 0.6;
numberOfPoints = 100;
% L = sqrt(numberOfPoints/rho);
L = 5;
rho = numberOfPoints/(L^2);
v = 0.05;
r = 0.2; %definition of neighbourhood for averaging
r0 = 0.3; %border offset for reflection

pointsX = L * rand(1, numberOfPoints);
pointsY = L * rand(1, numberOfPoints);
Nsteps=500;
theta = 2*pi*rand(1,numberOfPoints) - pi; %-pi to pi
timedelta = 1;

y = zeros(Nsteps+1,3*numberOfPoints); % initialisation of solution vector
y(1,:) = [pointsX,pointsY,theta];
for k=1:Nsteps
    posX=y(k,1:numberOfPoints); % x position matrix
    posY=y(k,numberOfPoints+1:2*numberOfPoints); % y position matrix
    theta=y(k,2*numberOfPoints+1:end); % theta matrix

    newTheta=theta; %just arbit initialisation, supposed to be faster this way
    vel=[v*cos(theta) v*sin(theta)];
    for i=1:numberOfPoints %evaluating updated angle parameter
        if L-posX(i)<r0 ||posX(i)<r0||L-posY(i)<r0||posY(i)<r0
            continue
        end
        avgcos = 0;
        avgsin = 0;
        ctr = 0;
        for j=1:numberOfPoints
            if L-posX(j)<r0 ||posX(j)<r0||L-posY(j)<r0||posY(j)<r0
                continue
            end
            dx = posX(i) - posX(j);
            dy = posY(i) - posY(j);
            if (dx^2+dy^2)<(r^2)
                avgcos = avgcos + cos(theta(j));
                avgsin = avgsin + sin(theta(j));
                ctr = ctr + 1;
            end
        end
        avgcos = avgcos/ctr;
        avgsin = avgsin/ctr;
        angleMean = atan2(avgsin,avgcos);
        newTheta(i) = angleMean + (noise*rand - noise/2);
    end
    theta = newTheta; % update angle values
    posX = posX + vel(1:numberOfPoints)*timedelta;
    posY = posY + vel(numberOfPoints+1:end)*timedelta;
    for i = 1:numberOfPoints
        if (L-posX(i)<r0 && cos(theta(i))>0)||(posX(i)<r0 && cos(theta(i))<0)
            %theta(i) = sign( sin(theta(i)) )*pi-theta(i);
            theta(i) = atan2(sin(theta(i)),-cos(theta(i)));
        end
        if (L-posY(i)<r0 && sin(theta(i))>0)||(posY(i)<r0 && sin(theta(i))<0)
            theta(i) = atan2(-sin(theta(i)),cos(theta(i)));
            %theta(i) = -theta(i)
        end
    end
    y(k+1,:)=[posX,posY,theta];
end

%movie
for k=1:size(y)
	posX=y(k,1:numberOfPoints); % x position matrix
	posY=y(k,numberOfPoints+1:2*numberOfPoints); % y position matrix
	plot(posX,posY,'b.'); %plot instantaneous position
    axis([0,L,0,L]);
    axis square;
    grid on;
	pause(.1);
end