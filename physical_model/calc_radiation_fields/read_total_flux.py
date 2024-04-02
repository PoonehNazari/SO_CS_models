import numpy as np
from radmc3dPy.analyze import *
import matplotlib.pylab as plt
from radmc3dPy.natconst import *

nx =1005
ny =400
nz = 1
p_mas = 1.67e-24 #gr
#analyze.radmc3dData(grid)
cmap = plt.get_cmap('jet')
cmap.set_bad(color = 'k')
plt.figure()
a    = readData(ddens=True,dtemp=True)
r    = a.grid.x[:]
theta = a.grid.y[:]
theta[0]=0.
theta[-1]=np.pi
r2,theta2 = np.meshgrid(r/au,theta)
dens_2d = a.rhodust[:,:,0,0].T * 100./(1.4* p_mas)
temp2d = a.dusttemp[:,:,0,0].T
dens_2d_masked = np.ma.masked_where(dens_2d<1.05e3,dens_2d)
temp2d_masked = np.ma.masked_where(dens_2d<1.05e3,temp2d)
c = plt.pcolormesh(r2 * np.sin(theta2), r2 * np.cos(theta2),
    np.log10(temp2d_masked), cmap=cmap)
plt.xlim(0,80)
plt.ylim(0,40)


data = np.loadtxt('userdef_total_flux.out',skiprows=2)
hdr = np.loadtxt('userdef_total_flux.out')[0:2]
if len(data) != hdr[1]:
    print('Your data shape is different from what the header says')

data = np.reshape(data, [1,nz, ny, nx])
data = np.swapaxes(data, 0, 3)
data = np.swapaxes(data, 1, 2)
data = data[:,:,0,0].T
data[np.where(data==0.0000000000000000)] = 1e-100
plt.figure()
cmap = plt.get_cmap('tab20')
cmap.set_bad(color = 'k')
data_masked = np.ma.masked_where(dens_2d<1.05e3,data)
outflow_elements = np.where(dens_2d<1.05e3)
data_elemnts = np.where(data==0.0000000000000000)
print('lalalalallalalallalla',data_elemnts)
c = plt.pcolormesh(r2 * np.sin(theta2), r2 * np.cos(theta2),
    np.log10(data_masked), cmap=cmap, vmax = 6,vmin=-4)
plt.colorbar()
plt.xlim(0,1000)
plt.ylim(0,1000)
plt.show()
