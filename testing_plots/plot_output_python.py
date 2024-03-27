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

point1 = np.loadtxt('test_point_oxygen.dat')

time = point1[:,0]
gco = point1[:,3]
gh2co = point1[:,4]
gch3oh = point1[:,5]
gch4 = point1[:,6]
gh2o = point1[:,7]
gch3och3 = point1[:,8]
gch3cho = point1[:,9]
ghcooh = point1[:,10]

ax.plot(time,gco,linestyle='solid',color=cmap(0),lw='2',label='CO')
ax.plot(time,gh2co,linestyle='solid',color=cmap(0.1),lw='2',label='$\mathregular{H_2CO}$')
ax.plot(time,gch3oh,linestyle='solid',color=cmap(0.2),lw='2',label='$\mathregular{CH_3OH}$')
ax.plot(time,gch4,linestyle='solid',color=cmap(0.3),lw='2',label='$\mathregular{CH_4}$')
ax.plot(time,gh2o,linestyle='solid',color=cmap(0.4),lw='2',label='$\mathregular{H_2O}$')
ax.plot(time,gch3och3,linestyle='solid',color=cmap(0.5),lw='2',label='$\mathregular{CH_3OCH_3}$')
ax.plot(time,gch3cho,linestyle='solid',color=cmap(0.6),lw='2',label='$\mathregular{CH_3CHO}$')
ax.plot(time,ghcooh ,linestyle='solid',color=cmap(0.7),lw='2',label='$\mathregular{HCOOH}$')


plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=4,frameon=False,handlelength=1)


point1 = np.loadtxt('test_point_oxygen_uv.dat')

time = point1[:,0]
gco = point1[:,3]
gh2co = point1[:,4]
gch3oh = point1[:,5]
gch4 = point1[:,6]
gh2o = point1[:,7]
gch3och3 = point1[:,8]
gch3cho = point1[:,9]
ghcooh = point1[:,10]

ax.plot(time,gco,linestyle='dashed',color=cmap(0),lw='2',label='CO')
ax.plot(time,gh2co,linestyle='dashed',color=cmap(0.1),lw='2',label='$\mathregular{H_2CO}$')
ax.plot(time,gch3oh,linestyle='dashed',color=cmap(0.2),lw='2',label='$\mathregular{CH_3OH}$')
ax.plot(time,gch4,linestyle='dashed',color=cmap(0.3),lw='2',label='$\mathregular{CH_4}$')
ax.plot(time,gh2o,linestyle='dashed',color=cmap(0.4),lw='2',label='$\mathregular{H_2O}$')
ax.plot(time,gch3och3,linestyle='dashed',color=cmap(0.5),lw='2',label='$\mathregular{CH_3OCH_3}$')
ax.plot(time,gch3cho,linestyle='dashed',color=cmap(0.6),lw='2',label='$\mathregular{CH_3CHO}$')
ax.plot(time,ghcooh ,linestyle='dashed',color=cmap(0.7),lw='2',label='$\mathregular{HCOOH}$')



ax.set_xlim(10**3,1*10**7)
ax.set_ylim(1*10**-12,2*10**-3)
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(False)

# ax.set_title('Test',y=1.12)
ax.set_xlabel('Time (yr)',color='k')
ax.set_ylabel('Abundance wrt $\mathregular{n_H}$ ($\mathregular{cm-3}$)',color='k')
ax.set_title('Solid: UV =O, Dashed: UV = 4Go',y=1.15)

plt.savefig('test_point_oxygen.pdf',dpi=300)

##########


fig = plt.figure(figsize=[6,5],constrained_layout=True)
ax = plt.subplot(111)

point1 = np.loadtxt('test_point_nitrogen.dat')

time = point1[:,0]
gco = point1[:,3]
gh2co = point1[:,4]
gch3oh = point1[:,5]
gch4 = point1[:,6]
gh2o = point1[:,7]
gch3cn = point1[:,8]
ghnco = point1[:,9]
gnh2cho = point1[:,10]

ax.plot(time,gco,linestyle='solid',color=cmap(0),lw='2',label='CO')
ax.plot(time,gh2co,linestyle='solid',color=cmap(0.1),lw='2',label='$\mathregular{H_2CO}$')
ax.plot(time,gch3oh,linestyle='solid',color=cmap(0.2),lw='2',label='$\mathregular{CH_3OH}$')
ax.plot(time,gch4,linestyle='solid',color=cmap(0.3),lw='2',label='$\mathregular{CH_4}$')
ax.plot(time,gh2o,linestyle='solid',color=cmap(0.4),lw='2',label='$\mathregular{H_2O}$')
ax.plot(time,gch3cn,linestyle='solid',color=cmap(0.5),lw='2',label='$\mathregular{CH_3CN}$')
ax.plot(time,ghnco,linestyle='solid',color=cmap(0.6),lw='2',label='$\mathregular{HNCO}$')
ax.plot(time,gnh2cho ,linestyle='solid',color=cmap(0.7),lw='2',label='$\mathregular{NH_2CHO}$')


plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=4,frameon=False,handlelength=1)


point1 = np.loadtxt('test_point_nitrogen_uv.dat')

time = point1[:,0]
gco = point1[:,3]
gh2co = point1[:,4]
gch3oh = point1[:,5]
gch4 = point1[:,6]
gh2o = point1[:,7]
gch3cn = point1[:,8]
ghnco = point1[:,9]
gnh2cho = point1[:,10]

ax.plot(time,gco,linestyle='dashed',color=cmap(0),lw='2',label='CO')
ax.plot(time,gh2co,linestyle='dashed',color=cmap(0.1),lw='2',label='$\mathregular{H_2CO}$')
ax.plot(time,gch3oh,linestyle='dashed',color=cmap(0.2),lw='2',label='$\mathregular{CH_3OH}$')
ax.plot(time,gch4,linestyle='dashed',color=cmap(0.3),lw='2',label='$\mathregular{CH_4}$')
ax.plot(time,gh2o,linestyle='dashed',color=cmap(0.4),lw='2',label='$\mathregular{H_2O}$')
ax.plot(time,gch3cn,linestyle='dashed',color=cmap(0.5),lw='2',label='$\mathregular{CH_3CN}$')
ax.plot(time,ghnco,linestyle='dashed',color=cmap(0.6),lw='2',label='$\mathregular{HNCO}$')
ax.plot(time,gnh2cho ,linestyle='dashed',color=cmap(0.7),lw='2',label='$\mathregular{NH_2CHO}$')



ax.set_xlim(10**3,1*10**7)
ax.set_ylim(1*10**-12,2*10**-3)
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(False)

ax.set_title('Solid: UV =O, Dashed: UV = 4Go',y=1.15)
ax.set_xlabel('Time (yr)',color='k')
ax.set_ylabel('Abundance wrt $\mathregular{n_H}$ ($\mathregular{cm-3}$)',color='k')

plt.savefig('test_point_nitrogen.pdf',dpi=300)
plt.show()
