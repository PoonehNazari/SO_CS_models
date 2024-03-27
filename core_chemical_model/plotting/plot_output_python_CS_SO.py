import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib
import matplotlib.cm as cm
import matplotlib.font_manager
from matplotlib.ticker import LogFormatter 
from matplotlib.pylab import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.interpolate import griddata
from matplotlib.colors import ListedColormap
#import cmasher as cmr

cmap = cm.jet

matplotlib.rc('font', family='sans-serif') 
matplotlib.rc('font', serif='Helvetica Neue') 
matplotlib.rc('text', usetex='false') 
matplotlib.rcParams.update({'font.size': 12})

fig = plt.figure(figsize=[6,5],constrained_layout=True)
ax = plt.subplot(111)

point1 = np.loadtxt('test_point_lowC_G01.dat')

time = point1[:,0]
SO = point1[:,3]
CS = point1[:,4]
C2H = point1[:,5]
CH3OH = point1[:,6]
CH3CN = point1[:,7]

ax.plot(time,SO,linestyle='solid',color=cmap(0),lw='2',label='SO')
ax.plot(time,CS,linestyle='solid',color=cmap(0.5),lw='2',label='CS')
ax.plot(time,C2H,linestyle='solid',color=cmap(0.7),lw='2',label='C2H')
ax.plot(time,CH3OH,linestyle='solid',color=cmap(0.8),lw='2',label='CH3OH')
ax.plot(time,CH3CN,linestyle='solid',color=cmap(0.9),lw='2',label='CH3CN')

plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=4,frameon=False,handlelength=1)



#ax.set_xlim(10**3,1*10**7)
#ax.set_ylim(1*10**-12,2*10**-3)
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(True)

# ax.set_title('Test',y=1.12)
ax.set_xlabel('Time (yr)',color='k')
ax.set_ylabel('Abundance wrt $\mathregular{n_H}$ ($\mathregular{cm-3}$)',color='k')
#ax.set_title('Solid: UV =O, Dashed: UV = 4Go',y=1.15)

plt.savefig('test_point_lowC_G01.pdf',dpi=300)
plt.show()
##########


