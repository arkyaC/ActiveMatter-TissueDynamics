#incomplete
import numpy as np
import matplotlib.pyplot as plt
import time

L=14;driftspeed = 1;mu=1;tau=1;Req=5.0/6;R0=1;Fadh=0.75;Frep=30 #parameters
noise=.6 #stochastic noise parameter
rhoNorm = 0.1 #parameter for normalized density
x =np.random.uniform(0,L,10000)
y =np.random.uniform(0,L,10000)
minAllowableDistance = Req
numberOfPoints = int(np.floor(2*rhoNorm* L**2 /(3.14* Req**2)))

# Initialize first point.
keeperX = keeperY = np.empty(1)
keeperY[0] = y[0]
keeperX[0] = x[0]
# Try dropping down more points.
counter = 1
k=1
while counter<numberOfPoints:
	thisX = x[k]*np.ones_like(keeperX)
	thisY = y[k]*np.ones_like(keeperY)
	# See how far is this away from existing keeper points.
	distances=np.sqrt((thisX - keeperX)**2 + (thisY - keeperY)**2)
	minDistance = np.min(distances)
	if minDistance >= minAllowableDistance:
	   X = thisX[0];Y = thisY[0]
	   keeperX = np.append(keeperX,X)
	   keeperY = np.append(keeperY,Y)
	   counter = counter + 1
	k=k+1
Nsteps=4000;cutoffIter=1500
theta = 2 * np.pi * np.random.uniform(0,1,numberOfPoints)
timedelta=0.05*R0/driftspeed

y=np.zeros((Nsteps+1,3*numberOfPoints))
y[0,:] = np.concatenate((keeperX,keeperY,theta)) #ICs
correl=0 #init order param
for k in xrange(Nsteps):
	posX=y[k,0 : numberOfPoints] # x position matrix
	posY=y[k,numberOfPoints: 2*numberOfPoints] # y position matrix
	theta=y[k,2*numberOfPoints:]
	drift=np.concatenate((driftspeed*np.cos(theta),driftspeed*np.sin(theta)))
	rhs = np.empty_like(drift)
	for i in xrange(len(posX)):
		interaxn=np.array([0,0])
		for j in xrange(len(posX)): #force of interaction on ith particle
			if j==i:
				continue
			dij=np.sqrt((posX[i]-posX[j])**2 + (posY[i] - posY[j])**2)
			eij=np.array([ posX[i]-posX[j],posY[i]-posY[j] ])/dij
			if dij<Req:
				interaxn = interaxn - Frep*(dij-Req)/Req * eij
			elif dij>=Req and dij<=R0:
				interaxn = interaxn - Fadh*(dij-Req)/(R0-Req) * eij
		interaxn=mu*interaxn
		rhs[i] = drift[i] + interaxn[0]
		rhs[numberOfPoints+i] = drift[numberOfPoints+i] + interaxn[1]
	s=np.array([0,0])
	res = np.empty_like(theta)
	for i in xrange(len(theta)):
		ni=np.array([drift[i], drift[numberOfPoints+i], 0])/driftspeed
		vi=np.array([rhs[i], rhs[numberOfPoints+i], 0])
		normvi=np.sqrt(vi[0]**2 + vi[1]**2)
		vi=vi/normvi
		s=s+np.array([vi[0],vi[1]])
		ez=np.array([0,0,1])
		res[i]=np.arcsin(np.inner(ez,np.cross(ni,vi)))
	if k>cutoffIter:
		correl=correl+np.sqrt(s[0]**2+s[1]**2)/numberOfPoints #order parameter
	res=res/tau
	posX=posX+timedelta*rhs[0:numberOfPoints]
	posY=posY+timedelta*rhs[numberOfPoints:]
	posX=L*(posX/L-np.floor(posX/L)) # periodic BC
	posY=L*(posY/L-np.floor(posY/L)) # periodic BC
	theta=theta + timedelta*res + np.sqrt(timedelta*(noise**2)/12) * np.random.normal(0,1,len(theta))
	y[k+1,:]=np.concatenate((posX,posY,theta))#Euler stepping

correl=correl/(Nsteps-cutoffIter)
print correl

plt.axis([0, L, 0, L])
plt.ion()
for i in xrange(len(y)): #problem with animation
	plt.axis([0, L, 0, L])
	plt.scatter(y[i,0 : numberOfPoints],y[i,numberOfPoints: 2*numberOfPoints])#,'b.')
	plt.pause(0.00001)
	plt.gcf().clear()
