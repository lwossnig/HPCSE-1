#
# Configuration plot script for the HPCSE I course
# (c) 2014 Andreas Hehn <hehn@phys.ethz.ch>, ETH Zurich
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import re
import sys

if len(sys.argv) != 2:
    print "usage:", sys.argv[0], "<sim_output_file>"
    sys.exit(1)

filename = sys.argv[1]

def get_particlenumber(filename):
    f = open(filename,'r')
    f.readline()
    scndline = f.readline()
    f.close()
    match = re.search('nparticles = (\d+)',scndline)
    if match:
        return int(match.group(1))
    else:
        raise RuntimeError("nparticles not found in second line of input file!")

fig = plt.figure()
nparticles= get_particlenumber(filename)
fcontent = np.loadtxt(filename, float)
num_per_row = len(fcontent[0])
fcontent = np.reshape(fcontent, (-1,nparticles,num_per_row))

print "Plotting", len(fcontent), "frames with", nparticles, "particles..."
data = []
maxforce = 0.0
for ap in fcontent:
    a = ap.transpose()
    force = np.sqrt(a[4]*a[4]+a[5]*a[5])
    f = force.max()
    if f > maxforce:
        maxforce = f
    data += [np.array([a[0], a[1], a[2], a[3], force])]
 
ims = []
pd = data[0]
ppd = data[0]
for d in data:
    ims.append( (plt.scatter(ppd[0],ppd[1], c=ppd[4]/maxforce, cmap=plt.cm.cool, alpha=0.10), plt.scatter(pd[0],pd[1], c=pd[4]/maxforce, cmap=plt.cm.cool, alpha=0.25), plt.scatter(d[0],d[1], c=d[4]/maxforce, cmap=plt.cm.cool),))
    ppd = pd
    pd = d

plt.xlim(0,1)
plt.ylim(0,1)

ani = animation.ArtistAnimation(fig, ims, interval=10, repeat=True)
plt.show()
