import numpy as np
import matplotlib.pyplot as plt
import time

numberOfPoints = 4000
rho = 10.4;noise = .4 #parameters
L = np.sqrt(numberOfPoints/rho)
#noise = .4
v = 0.03 #as stated in the paper (for optimum results)
r = 1 #definition of neighbourhood for averaging
timedelta = 1 #as mentioned in paper

pointsX =np.random.uniform(0,L,numberOfPoints)
pointsY =np.random.uniform(0,L,numberOfPoints)
Nsteps=1500
theta = np.random.uniform(-np.pi,np.pi,numberOfPoints) #-pi to pi

#y = zeros(Nsteps+1,3*numberOfPoints)
#y(1,:) = [pointsX,pointsY,theta]
solX = np.empty([Nsteps+1,numberOfPoints]);solX[0][:] = pointsX
solY = np.empty([Nsteps+1,numberOfPoints]);solY[0][:] = pointsY
solTheta = np.empty([Nsteps+1,numberOfPoints]);solTheta[0][:] = theta
orderN = np.empty(Nsteps)

for k in range(Nsteps):
	#posX=solX[k][:] # x position matrix
	#posY=solY[k][:] # y position matrix
	#theta=solTheta[k][:] # theta matrix

    #newTheta=theta #just arbit initialisation, supposed to be faster this way
	vel=[v*cos(theta) v*sin(theta)]
    s = [0,0] #avg vel init
    for i in range(numberOfPoints): #evaluating updated angle parameter
        avg = 0
        ctr = 0
        ni = [vel(i) vel(numberOfPoints+i)]/v
        s = s + ni
        for j in range(numberOfPoints):
            distIJsq = (posX(j)-posX(i))**2 + (posY(j) - posY(i))**2
            if distIJsq<=(r**2):
                avg = avg + theta(j)
                ctr = ctr + 1


        avg = avg/ctr
        temp = avg + (noise*rand - noise/2)
        #keep angle between -pi and pi
        newTheta(i) = temp-2*pi*floor((temp+pi)/(2*pi))

    orderN(k)=sqrt(s(1)**2+s(2)**2)/numberOfPoints

    theta = newTheta # update angle values
    posX = posX + vel(1:numberOfPoints)*timedelta
    posY = posY + vel(numberOfPoints+1:end)*timedelta
    posX = posX-L*floor(posX/L) # periodic BC
	posY = posY-L*floor(posY/L)
    y(k+1,:)=[posX,posY,theta]

# write data to dump
f = open('dump.txt','w')
f.write('Time \tOrder \n')
for k in range(Nsteps):
	f.write('%d \t%6.5f \n'%((k+1) , orderN[k]))
f.close()

#plot order parameter against time
plot(linspace(0,Nsteps,Nsteps),orderN)
axis([0,Nsteps,0,1])
xlabel('Time step');ylabel('Order Parameter')

# pause(10)

#movie
# for k=1:size(y)
# 	posX=y(k,1:numberOfPoints) # x position matrix
# 	posY=y(k,numberOfPoints+1:2*numberOfPoints) # y position matrix
# 	plot(posX,posY,'b.') #plot instantaneous position
#     axis([0,L,0,L])
#     axis square
#     grid on
# 	pause(.1)
#
