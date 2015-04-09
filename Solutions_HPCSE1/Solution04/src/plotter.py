from math import pi
import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt('integral.dat', unpack=True)

x = data[0]
y = data[1]
y_errors = data[2]

# === approximations

plt.xscale('log')
plt.yscale('linear')

current_plot ,= plt.plot(x,y,'-o',markersize=4.0)
plt.errorbar(x,y,yerr=y_errors,fmt='o',markersize=4.0,color=current_plot.get_color())

plt.axhline(y = 0.3, linestyle='dashed', color='k')

plt.ylabel(r'$I$',fontsize=16)

plt.xlabel(r'$N$',fontsize=16)

plt.savefig('integral.pdf',dpi=600)

plt.show()

# === error convergence

plt.xscale('log')
plt.yscale('log')

current_plot ,= plt.plot(x,y_errors,'-o',markersize=4.0)

plt.ylabel(r'$\Delta I$',fontsize=16)

plt.xlabel(r'$N$',fontsize=16)

plt.savefig('integral_errors.pdf',dpi=600)

plt.show()

# === scaling

data = np.loadtxt('scaling.dat', unpack=True)

x = data[0]
y = data[1]
y_errors = data[2]

y = y[0] / y

plt.xscale('linear')
plt.yscale('linear')

current_plot ,= plt.plot(x,y,'-o',label='N = 2^17 (static)',markersize=4.0)
plt.legend (loc = 'best')
xmin, xmax = plt.xlim()
plt.ylim( (xmin, xmax) )

plt.ylabel(r'$T_1/T_N$',fontsize=16)

plt.xlabel(r'#threads',fontsize=16)

plt.savefig('scaling.pdf',dpi=600)

plt.show()

