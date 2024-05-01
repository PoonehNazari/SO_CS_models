import matplotlib
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 20
matplotlib.rcParams['lines.linewidth'] = 3.
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['legend.numpoints'] = 1
matplotlib.rcParams['legend.handlelength'] = 2
matplotlib.rcParams['legend.fontsize'] = 23
matplotlib.rcParams['lines.markeredgewidth'] = 2
matplotlib.rcParams['lines.markersize'] = 8
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.major.width'] = 2
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['xtick.minor.width'] = 2.
matplotlib.rcParams['xtick.minor.size'] = 4
matplotlib.rcParams['ytick.major.width'] = 2.
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['ytick.minor.width'] = 2.
matplotlib.rcParams['ytick.minor.size'] = 4
import numpy as np
import matplotlib.pylab as plt
import os
import sys
from scipy.interpolate import griddata
###########################################################################
#################################FUNCS/CLASSES#####################################
class My_source:
    def __init__(self, C_O, CS, SO, C2H, HNCO, r, theta, Mstar):
        self.C_O = C_O #input C/O
        self.CS = CS #final CS abundances
        self.SO = SO #final SO abundances
        self.C2H = C2H
        self.HNCO = HNCO
        self.r = r #array of radii
        self.theta = theta #array of theta
        self.Mstar = Mstar #the stellar mass

    def CS_SO(self):
        ratio = self.CS/self.SO
        return ratio

    def HNCO_C2H(self):
        ratio =  self.HNCO/self.C2H
        return ratio

    def convert_to_cylindrical(self):
        new_r = self.r * np.sin(self.theta)
        new_z = self.r * np.cos(self.theta)
        return new_r, new_z
###########################################################################
###########################################################################
def fill_bad_pixels(data,grid_r, grid_z):
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if np.isnan(data[i, j]) and (grid_r[i, j] < 45) and (grid_z[i, j] < 30):
                # Check left neighbor
                if j > 0 and not np.isnan(data[i, j-1]):
                    data[i, j] = data[i, j-1]
                # Check right neighbor if left is not available
                elif j < data.shape[1] - 1 and not np.isnan(data[i, j+1]):
                    data[i, j] = data[i, j+1]
                # elif i > 0 and not np.isnan(data[i-1, j]):
                #     data[i, j] = data[i-1, j]
                # # Try bottom neighbor
                # elif i < data.shape[0] - 1 and not np.isnan(data[i+1, j]):
                #     data[i, j] = data[i+1, j]
    return data
##########################################
##########################################
def read_output_data(M_star, C_O_ratio):
    outputs = ['output_points_0.dat', 'output_points_1.dat', 'output_points_2.dat', 'output_points_3.dat']#, 'output_points_4.dat']
    home_dir = '/Users/pooneh/Library/Mobile Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/'
    r = []
    theta = []
    CS = []
    SO = []
    C2H  = []
    HNCO = []
    for i_output in range(len(outputs)):
        output = open(home_dir+'outputs/'+M_star+'/C_O_'+str(C_O_ratio)+'/'+outputs[i_output])
        lines_output = output.readlines()
        output.close()
        for line in lines_output:
            if line.startswith(' RADIUS ='):
                current_radius = float(line.split('=')[1].strip().split()[0])
                r.append(current_radius)
            elif line.startswith(' HEIGHT ='):
                current_theta = float(line.split('=')[1].strip().split()[0])
                theta.append(current_theta)
            elif line.startswith(' CS '):
                parts_CS = line.split()
                parts_CS = ['1' if ('+1' in element) and ('E+1' not in element) else element for element in parts_CS]
                value_CS = [float(value_CS) for value_CS in parts_CS[1:-1]][-1]
                CS.append(value_CS)
            elif line.startswith(' SO '):
                parts_SO = line.split()
                parts_SO = ['1' if ('+1' in element) and ('E+1' not in element) else element for element in parts_SO]
                value_SO = [float(value_SO) for value_SO in parts_SO[1:-1]][-1]
                SO.append(value_SO)
            elif line.startswith(' C2H '):
                parts_C2H = line.split()
                parts_C2H = ['1' if ('+1' in element) and ('E+1' not in element) else element for element in parts_C2H]
                value_C2H = [float(value_C2H) for value_C2H in parts_C2H[1:-1]][-1]
                C2H.append(value_C2H)
            elif line.startswith(' HNCO '):
                parts_HNCO = line.split()
                parts_HNCO = ['1' if ('+1' in element) and ('E+1' not in element) else element for element in parts_HNCO]
                value_HNCO = [float(value_HNCO) for value_HNCO in parts_HNCO[1:-1]][-1]
                HNCO.append(value_HNCO)

    SO = np.array(SO)
    CS = np.array(CS)
    C2H = np.array(C2H)
    HNCO = np.array(HNCO)
    r = np.array(r)
    theta = np.array(theta)
    if len(theta) != len(r):
        print('ERROR the lengths are not the same something is not being found')
        sys.exit()
    if len(SO) != len(CS):
        print('ERROR the lengths are not the same something is not being found')
        sys.exit()
    if len(r) != len(CS):
        print('ERROR the lengths are not the same something is not being found')
        sys.exit()
    source = My_source(C_O_ratio, CS, SO, C2H,HNCO, r, theta, M_star)
    return source
###########################################################################
###########################################################################
def plot_maps_write_CS_SO(M_star, C_O_ratio):
    source = read_output_data(M_star, C_O_ratio)
    new_r, new_z = source.convert_to_cylindrical()
    CS_SO = source.CS_SO()
    HNCO_C2H = source.HNCO_C2H()
    CS = source.CS
    SO =  source.SO
    C2H = source.C2H
    HNCO = source.HNCO

    cmap = plt.get_cmap('magma')
    #cmap.set_bad(color = 'k')
    vmin = -14
    vmax = -7
    fig,ax = plt.subplots(1,3,figsize=(15,4))
    grid_r, grid_z = np.meshgrid(np.logspace(np.log10(min(new_r)),np.log10(max(new_r)),150),
                                np.linspace(min(new_z),max(new_z),100))


    grid_C2H = griddata((new_r, new_z), np.log10(C2H), (grid_r, grid_z), method='linear')
    grid_SO = griddata((new_r, new_z), np.log10(SO), (grid_r, grid_z), method='linear')
    grid_CS = griddata((new_r, new_z), np.log10(CS), (grid_r, grid_z), method='linear')
    cax1 = ax[0].scatter(new_r, new_z, c=np.log10(C2H), cmap= cmap, vmin=vmin, vmax = vmax,marker='s',s=80)
    cax2 = ax[1].scatter(new_r, new_z, c=np.log10(SO), cmap= cmap, vmin=vmin, vmax = vmax,marker='s',s=80)
    cax3 = ax[2].scatter(new_r, new_z, c=np.log10(CS), cmap= cmap, vmin=vmin, vmax = vmax,marker='s',s=80)
    #cax1 = ax[0].pcolormesh(grid_r, grid_z,grid_C2H , cmap=cmap,vmin=vmin, vmax = vmax)
    #cax2 = ax[1].pcolormesh(grid_r, grid_z,grid_CS , cmap=cmap,vmin=vmin, vmax = vmax)
    #cax3 = ax[2].pcolormesh(grid_r, grid_z,grid_SO , cmap=cmap,vmin=vmin, vmax = vmax)
    #ax[0,0].annotate('C/O = '+str(C_O_ratio), xycoords='axes fraction', xy=(0.02,0.9),
    #              weight='bold',fontsize=18)
    ax[0].set_xlabel('R [au]')
    ax[1].set_xlabel('R [au]')
    ax[2].set_xlabel('R [au]')
    ax[0].set_ylabel('z [au]')

    cb1 = fig.colorbar(cax1, ax=ax[0])
    cb1.set_label('log(C$_2$H)')
    cb2 = fig.colorbar(cax2, ax=ax[1])
    cb2.set_label('log(CS)')
    cb3 = fig.colorbar(cax3, ax=ax[2])
    cb3.set_label('log(SO)')
    plt.tight_layout()
    plt.savefig('ALMA_proposal/maps'+str(C_O_ratio)+'_'+str(M_star)+'.png')

    matplotlib.rcParams['font.size'] = 25
    plt.figure()
    #grid_CS_SO = griddata((new_r, new_z), np.log10(CS_SO), (grid_r, grid_z), method='linear')
    # Fill bad pixels
    #filled_data = fill_bad_pixels(grid_CS_SO,grid_r, grid_z)
    #plt.pcolormesh(grid_r, grid_z,filled_data , cmap=cmap,vmin=-3, vmax = 2)
    plt.scatter(new_r, new_z, c=np.log10(CS_SO), cmap= cmap, vmin=-3, vmax = 3,marker='s',s=80)
    plt.annotate('C/O = '+str(C_O_ratio), xycoords='axes fraction', xy=(0.02,0.9),
                  weight='bold',fontsize=18)
    plt.xlabel('R [au]')
    plt.ylabel('z [au]')
    plt.colorbar(label='log(CS/SO)')
    plt.tight_layout()
    plt.savefig('ALMA_proposal/map_C_O_'+str(C_O_ratio)+'_'+str(M_star)+'_ALMA.png')

###########################################################################
###########################################################################
#########################MAIN##############################################
C_O_ratios = [0.2,0.44, 0.9,1.2]
M_stars = ['M_star_0.5']
for i_star in range(len(M_stars)):
    for i_ratios in range(len(C_O_ratios)):
        plot_maps_write_CS_SO(M_stars[i_star], C_O_ratios[i_ratios])
