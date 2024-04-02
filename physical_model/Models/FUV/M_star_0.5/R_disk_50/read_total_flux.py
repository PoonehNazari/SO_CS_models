import numpy as np
from radmc3dPy.analyze import *
import matplotlib.pylab as plt
from radmc3dPy.natconst import *

nx =1005
ny =400
nz = 1
#analyze.radmc3dData(grid)

plt.figure()
a    = readData(ddens=True)
r    = a.grid.x[:]
theta = a.grid.y[:]
theta[0]=0.
theta[-1]=np.pi
r2,theta2 = np.meshgrid(r/au,theta)
dens_2d = a.rhodust[:,:,0,0].T
cmap = plt.get_cmap('tab10')
c = plt.pcolormesh(r2 * np.sin(theta2), r2 * np.cos(theta2),
    np.log10(dens_2d), cmap=cmap)
plt.xlim(0,1000)
plt.ylim(0,1000)


data = np.loadtxt('userdef_total_flux.out',skiprows=2)
hdr = np.loadtxt('userdef_total_flux.out')[0:2]
if len(data) != hdr[1]:
    print('Your data shape is different from what the header says')

data = np.reshape(data, [1,nz, ny, nx])
data = np.swapaxes(data, 0, 3)
data = np.swapaxes(data, 1, 2)

plt.figure()
c = plt.pcolormesh(r2 * np.sin(theta2), r2 * np.cos(theta2),
    np.log10(data[:,:,0,0].T), cmap=cmap, vmax = 6,vmin=-4)
plt.colorbar()
plt.xlim(0,1000)
plt.ylim(0,1000)
plt.show()
