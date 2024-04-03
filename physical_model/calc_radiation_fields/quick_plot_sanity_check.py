import numpy as np
import matplotlib.pylab as plt

data = np.loadtxt('../../inputs/input_points_M_star_0.5.txt')
r = data[:,0]
theta=data[:,1]
dens = data[:,2]

plt.scatter(r * np.sin(theta), r * np.cos(theta), c=np.log10(dens), cmap='jet')
plt.ylim(0,200)
plt.show()
