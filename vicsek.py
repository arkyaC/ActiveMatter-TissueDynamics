import numpy as np
import time
start_time = time.time()

numberOfPoints = 40
rho = 4.0;noise = 5.0 #parameters
L = np.sqrt(numberOfPoints/rho)
#noise = .4
v = 0.03 #as stated in the paper (for optimum results)
r = 1.0 #definition of neighbourhood for averaging
timedelta = 1.0 #as mentioned in paper

pointsX =np.random.uniform(0,L,numberOfPoints)
pointsY =np.random.uniform(0,L,numberOfPoints)
Nsteps=1500
theta = np.random.uniform(-np.pi,np.pi,numberOfPoints) #-pi to pi

#y = zeros(Nsteps+1,3*numberOfPoints)
#y(1,:) = [pointsX,pointsY,theta]
solX = np.empty([Nsteps+1,numberOfPoints]);solX[0] = pointsX
solY = np.empty([Nsteps+1,numberOfPoints]);solY[0] = pointsY
solTheta = np.empty([Nsteps+1,numberOfPoints]);solTheta[0] = theta
orderN = np.empty(Nsteps)

for k in xrange(Nsteps):
	newTheta = np.empty(numberOfPoints)
	xvel = v*np.cos(solTheta[k])
	yvel = v*np.sin(solTheta[k])
	avgVel = np.array([0.0,0.0]) #avg vel init
	for i in xrange(numberOfPoints): #evaluating updated angle parameter
	    avg = 0.0;ctr = 0
	    ni = np.array([xvel[i]/v,yvel[i]/v])
	    avgVel = avgVel + ni
	    for j in xrange(numberOfPoints):
	        distIJsq = (solX[k][j]-solX[k][i])**2 + (solY[k][j] - solY[k][i])**2
	        if distIJsq<(r**2):
	        	avg = avg + solTheta[k][j]
	        	ctr = ctr + 1
		if ctr>0:
			avg = avg/ctr
		temp = avg + np.random.uniform(-noise/2,noise/2)
		#keep angle between -pi and pi
		newTheta[i] = temp-2*np.pi*np.floor((temp+np.pi)/(2*np.pi))

	orderN[k] = np.sqrt(avgVel[0]**2+avgVel[1]**2)/numberOfPoints
	posX = solX[k] + xvel*timedelta #updating position
	posY = solY[k] + yvel*timedelta
	solX[k+1] = posX-L*np.floor(posX/L) # periodic BC
	solY[k+1] = posY-L*np.floor(posY/L)
	solTheta[k+1] = newTheta #update direction values

# write data to dump
f = open('dump.txt','w')
f.write('Time \tOrder \n')
for k in xrange(Nsteps):
	f.write('%d \t%8.7f \n'%((k+1) , orderN[k]))
f.close()
print(np.mean(orderN[Nsteps-400:]))

#plot order parameter against time
#plot(linspace(0,Nsteps,Nsteps),orderN)
#axis([0,Nsteps,0,1])
#xlabel('Time step');ylabel('Order Parameter')
print("--- %s seconds ---" % (time.time() - start_time))
