import numpy as np
from scipy.integrate import odeint
from scipy.optimize import curve_fit
from numpy import genfromtxt

def f(x,a,b,c): # rhs for ODE solver
    return a*np.power(x,c)+b
	#return a*np.power(x,.5)+b

file = open('./data/R_gyr-mon.txt') 
data = np.genfromtxt(file, names=['monomer','radius'])

R_gyr = data['radius']
mon_units = data['monomer']

#popt, cov = curve_fit(f, mon_units, R_gyr, [.5,2.5])
popt, cov = curve_fit(f, mon_units, R_gyr, [.5,2.5,0.5]) # extract fit results
a_opt, b_opt, c_opt = popt
#a_opt, b_opt = popt

print("a = %g" % a_opt)
print("b = %g" % b_opt)
print("c = %g" % c_opt)

import matplotlib.pyplot as plt
xdata = np.linspace(mon_units[0],mon_units[len(mon_units)-1], 2000)
axes = plt.gca()

plt.xlabel("No. of monomers")
plt.ylabel("Radius of gyration")
plt.title("c_opt="+str(c_opt));
plt.plot(mon_units, R_gyr, 'g.', label='Simulation data')
plt.plot(xdata, f(xdata, a_opt, b_opt, c_opt), 'r-', label='Fitted curve a(x^c) + b')
#plt.plot(xdata, f(xdata, a_opt, b_opt), 'r-', label='Fitted curve a(x^0.5) + b')
plt.legend(loc='lower right')
plt.gcf().set_size_inches(6, 4)
#plt.savefig('out.png', dpi=96) #to save the fit result
plt.show()
