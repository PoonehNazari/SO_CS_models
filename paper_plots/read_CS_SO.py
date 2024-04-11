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
    def __init__(self, C_O, CS, SO, r, theta, Mstar):
        self.C_O = C_O #input C/O
        self.CS = CS #final CS abundances
        self.SO = SO #final SO abundances
        self.r = r #array of radii
        self.theta = theta #array of theta
        self.Mstar = Mstar #the stellar mass

    def CS_SO(self):
        ratio = self.CS/self.SO
        return ratio

    def convert_to_cylindrical(self):
        new_r = self.r * np.sin(self.theta)
        new_z = self.r * np.cos(self.theta)
        return new_r, new_z
###########################################################################
###########################################################################
def read_output_data(M_star, C_O_ratio):
    outputs = ['output_points_0.dat', 'output_points_1.dat', 'output_points_2.dat', 'output_points_3.dat', 'output_points_4.dat']
    home_dir = '/Users/pooneh/Library/Mobile Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/'
    r = []
    theta = []
    CS = []
    SO = []
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
                value_CS = [float(value_CS) for value_CS in parts_CS[1:-2]][-1]
                CS.append(value_CS)
            elif line.startswith(' SO '):
                parts_SO = line.split()
                value_SO = [float(value_SO) for value_SO in parts_SO[1:-2]][-1]
                SO.append(value_SO)

    SO = np.array(SO)
    CS = np.array(CS)
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
    source = My_source(C_O_ratio, CS, SO, r, theta, M_star)
    return source
###########################################################################
###########################################################################
def plot_maps_write_CS_SO(M_star, C_O_ratio):
    source = read_output_data(M_star, C_O_ratio)
    new_r, new_z = source.convert_to_cylindrical()
    CS_SO = source.CS_SO()
    CS = source.CS
    SO =  source.SO
    myfile = open('C_O_'+str(C_O_ratio)+'_'+M_star+'.dat','w')
    myfile.write('r_cyl[au]\t z_cyl[au]\t CS/SO \t CS \t SO \n')
    for i in range(len(new_r)):
        myfile.write(str(new_r[i])+'\t'+str(new_z[i])+'\t'+str(CS_SO[i])+'\t'+str(CS[i])+'\t'+str(SO[i])+'\n')
    myfile.close()

    plt.figure()
    cmap = plt.get_cmap('magma')
    plt.scatter(new_r, new_z, c=np.log10(CS_SO), cmap= cmap, vmin=-2, vmax = 1,marker='s',s=80)
    plt.annotate('C/O = '+str(C_O_ratio), xycoords='axes fraction', xy=(0.02,0.9),
                  weight='bold',fontsize=18)
    plt.xlabel('R [au]')
    plt.ylabel('z [au]')
    plt.colorbar(label='CS/SO')
    plt.tight_layout()
    plt.savefig('map_C_O_'+str(C_O_ratio)+'_'+str(M_star)+'.png')

    plt.figure()
    grid_r, grid_z = np.meshgrid(np.logspace(np.log10(min(new_r)),np.log10(max(new_r)),150),
                                np.linspace(min(new_z),max(new_z),100))
    grid_CS_SO = griddata((new_r, new_z), np.log10(CS_SO), (grid_r, grid_z), method='linear')
    plt.pcolormesh(grid_r, grid_z,grid_CS_SO , cmap=cmap,vmin=-2, vmax = 1)
    plt.annotate('C/O = '+str(C_O_ratio), xycoords='axes fraction', xy=(0.02,0.9),
                  weight='bold',fontsize=18)
    plt.xlabel('R [au]')
    plt.ylabel('z [au]')
    plt.colorbar(label='CS/SO')
    plt.tight_layout()
    plt.savefig('map_C_O_'+str(C_O_ratio)+'_'+str(M_star)+'_interpolator.png')
###########################################################################
###########################################################################
#########################MAIN##############################################
C_O_ratios = [0.2]
M_stars = ['M_star_0.5']
for i_star in range(len(M_stars)):
    for i_ratios in range(len(C_O_ratios)):
        plot_maps_write_CS_SO(M_stars[i_star], C_O_ratios[i_ratios])

# data_test = np.loadtxt(home_dir+'inputs/M_star_0.5/input_points_0.txt')
# RAD = data_test[:,0]
# THET = data_test[:,1]
# RAD = [format(num, ".3e") for num in RAD]
# THET = [format(num, ".3e") for num in THET]
# test1 = np.column_stack((RAD, THET))
# set_test1 = {tuple(row) for row in test1}
#print(set_test1)

# r = [format(num, ".3e") for num in r]
# theta = [format(num, ".3e") for num in theta]
# test2 = np.column_stack((r, theta))
# set_test2 = {tuple(row) for row in test2}
# #print(set_test2)
# difference = np.array(list(set_test1 - set_test2))
# #print(difference)
