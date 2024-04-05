import numpy as np
import matplotlib.pylab as plt

for i in range(5):
    data = np.loadtxt('../../inputs/M_star_1.0/input_points_'+str(i)+'.txt')
    r = data[:,0]
    theta=data[:,1]
    dens = data[:,2]

    plt.scatter(r * np.sin(theta), r * np.cos(theta), c=np.log10(dens), cmap='jet',vmin=6.5,vmax=10)
plt.ylim(0,200)
plt.show()
