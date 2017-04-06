%max allowed density is 1.39
% rho = (pi * Req^2 * N) / L^2
% delX^2 * N = L^2
% rhoMax = 1.39 when delX = 1.5*Req

function orderN=szabo_grid(rho,noise)

numberOfPoints = 100;N=numberOfPoints; % must be perfect square
v0 = 1; %self-propelling speed
tau = 1;mu = 1; %parameters
Req = 5/6; %equilibrium radius
R0 = 1; %cut-off radius

Fadh = 0.75;Frep = 30; %force parameters
% Fwall = 50;

L = Req * sqrt(pi*numberOfPoints/rho); %domain size (side length)
delX = L/sqrt(N); %grid size (length)

x = delX/2:delX:delX*(sqrt(N)-1)+delX/2; %x positions
keeperX = repmat(x,1,sqrt(N)); %sqrt(N) rows of particles parallel to x axis, listed as a 1d array; sqrt(N)*sqrt(N) = N entries
y = x; %y positions, identical distribution as for x
tempY = repmat(y,sqrt(N),1);
keeperY = reshape(tempY,1,N); %reshape sqrt(N)xsqrt(N) array into a 1d array, with N entries, each representing one particle

Nsteps=20000; %number of timeSteps
% cutoffIter=Nsteps-100; %not being used currently
%timedelta=0.05*R0/v0; %not being used currently
theta = 2*pi*rand(1,numberOfPoints) - pi; %need random direction of motion, in the range of -pi to pi
timedelta=0.005; %time stepping interval

y=zeros(Nsteps+1,3*numberOfPoints); % solution matrix, initialised with a bunch of zeros
y(1,:) = [keeperX,keeperY,theta]; % initial position
%----------------Explanation of the structure of the solution matrix    y: ---------------
%each row represents a single time instant
%first N entries in a row represent the N x-coordinates at that time instant
%next N entries are the corresponding N y-coordinates 
%last N entries are the corresponding angle that the velocity vector makes with the x axis
%-----------------------------------------------------------------------------------------

orderN = zeros(1,Nsteps); %stores order parameter at each time step
for k=1:Nsteps
	posX=y(k,1:numberOfPoints); % x position matrix
	posY=y(k,numberOfPoints+1:2*numberOfPoints); % y position matrix
	theta=y(k,2*numberOfPoints+1:end); % direction matrix
    
	drift=[v0*cos(theta) v0*sin(theta)]; 
    % first N entries store x components of velocity, next N store the corresponding y components
    
    rhs = zeros(1,2*numberOfPoints); % rhs of the position updation Diff. eqn
    % structure of rhs exactly identical to that of drift
    % first N entries store x components of force, next N store the corresponding y components
    
    for i=1:numberOfPoints %iterating over each point
        interaxn=[0,0]; %initialise force of interaction, i.e. [Fx Fy]
        for j=1:numberOfPoints % to calculate the interaction force on i due to j
            if j==i
                continue;
            end
            dx = posX(i) - posX(j);
			dy = posY(i) - posY(j);
            %periodicity conditions
            if abs(dx) > 0.5*L
				dx = dx - L*sign(dx);
            end
            if abs(dy) > 0.5*L
				dy = dy - L*sign(dy);
            end
            
            dijSq = dx^2 + dy^2; % square of distance between i and j
            dij = sqrt(dijSq);
            eij = [dx,dy]/dij; %unit vector e_{ij} joining i and j
            if dijSq<Req^2 %attractive interaction force
                interaxn = interaxn - Frep*(dij-Req)/Req * eij;
            elseif dijSq>=Req^2 && dijSq<=R0^2 %repulsive interaction force
                interaxn = interaxn - Fadh*(dij-Req)/(R0-Req) * eij;
            end
        end
    %    %for wall force (in case wall is included)
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
        interaxn=mu*interaxn; % multiply by mobility parameter
        rhs(i)=drift(i)+interaxn(1); %x components, i.e. first N entries of rhs
        rhs(numberOfPoints+i)=drift(numberOfPoints+i)+interaxn(2); % y components
    end
    s=[0,0]; %avg velocity initialise
    res = zeros(1,numberOfPoints); %rhs for the angle updation diff eqn
    for i=1:numberOfPoints
        ni=[drift(i) drift(numberOfPoints+i)]/v0; %unit vector in direction of self-propulsion of i-th particle
        vi=[rhs(i) rhs(numberOfPoints+i)]; %velocity vector of i-th particle
        normvi=sqrt(vi(1)^2+vi(2)^2);
        vi=vi/normvi; %unit vector in the direction of velocity of i-th particle
        s=s+[vi(1) vi(2)]; %calculation of avg velocity (for calculation of order parameter)
        res(i)=asin(ni(1)*vi(2)-ni(2)*vi(1)); %calculation of the rhs of i-th particel theta diff. eqn
    end
    orderN(k)=sqrt(s(1)^2+s(2)^2)/numberOfPoints; %order parameter at the k-th step
	res = (1/tau)*res; %parameter tau being factored in
    %updation of position as per Euler stepping scheme
	posX = posX + timedelta*rhs(1:numberOfPoints);
    posY = posY + timedelta*rhs(numberOfPoints+1:end);
    %periodic boundary conditions
	posX = posX - L * floor(posX/L);
	posY = posY - L * floor(posY/L);
% 	theta = theta + timedelta*res + sqrt(timedelta*(noise^2)/12)*normrnd(0,1,[1 length(theta)]);
    noiseAmp = noise/sqrt(timedelta); %amplitude of noise
    theta = theta + timedelta * (res + noiseAmp * (rand(1,numberOfPoints) - 0.5)); %updation of theta, i.e. direction
	y(k+1,:) = [posX,posY,theta]; %Euler stepping in the solution matrix
    
%     if(mod(k,200) == 0) 
%         k %display status of execution
%     end
end

% plotting the order parameter value against step number
plot(linspace(0,Nsteps,Nsteps),orderN);
axis([0,Nsteps,0,1]);
xlabel('Time step');ylabel('Order Parameter');
% drawnow

% % write data to ./dump.txt
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
