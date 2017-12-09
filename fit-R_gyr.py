import numpy as np
from scipy.integrate import odeint
from scipy.optimize import curve_fit
from numpy import genfromtxt

def f(x,b,c): # rhs for ODE solver
    return c*x+b

file = open('./data/e2e_dist-mon_500.txt') #enter correct filename HERE
data = np.genfromtxt(file, names=['monomer','radius'])

R_gyr = data['radius']
mon_units = data['monomer']

#popt, cov = curve_fit(f, mon_units, R_gyr, [.5,2.5])
popt, cov = curve_fit(f, mon_units, R_gyr, [2.5,0.5]) # extract fit results
b_opt, c_opt = popt

print("b = %g" % b_opt)
print("c = %g" % c_opt)

import matplotlib.pyplot as plt
xdata = np.linspace(mon_units[0],mon_units[len(mon_units)-1], 2000)
axes = plt.gca()

plt.xlabel("No. of monomers")
# plt.ylabel("Radius of gyration sq") #for radius of gyration
plt.ylabel("End to end distance sq") #for end to end distance
plt.title("slope="+str(c_opt));
plt.plot(mon_units, R_gyr, 'g.', label='Simulation data')
plt.plot(xdata, f(xdata, b_opt, c_opt), 'r-', label='Fitted curve cx + b')
plt.legend(loc='lower right')
plt.gcf().set_size_inches(6, 4)
#plt.savefig('out.png', dpi=96) #to save the fit result
plt.show()
