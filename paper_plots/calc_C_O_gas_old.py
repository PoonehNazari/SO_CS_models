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
    def __init__(self, C_O, CO, C, CO2, CH4, CH3, O, O2, OH, NO, H2O, CS, SO, C2H2, H2CO, HCN, NH2CHO, C2H6, HNC, CH3NH2, C2H4, CN, C3, C3H, C3H2, CH3OH, CH3CCH , r, theta, Mstar):
        self.C_O = C_O #input C/O
        self.CO = CO
        self.C   = C
        self.CO2 = CO2
        self.CH4   = CH4
        self.CH3   = CH3
        self.O   = O
        self.O2   = O2
        self.OH   = OH
        self.NO   = NO
        self.H2O   = H2O
        self.CS = CS
        self.SO = SO
        self.C2H2 = C2H2
        self.H2CO = H2CO
        self.HCN = HCN
        self.NH2CHO = NH2CHO
        self.C2H6 = C2H6
        self.HNC = HNC
        self.CH3NH2  = CH3NH2
        self.C2H4 = C2H4
        self.CN  = CN
        self.C3  = C3
        self.C3H  = C3H
        self.C3H2 = C3H2
        self.CH3OH  = CH3OH
        self.CH3CCH = CH3CCH
        self.r = r #array of radii
        self.theta = theta #array of theta
        self.Mstar = Mstar #the stellar mass

    def calc_C_O_ratio(self):
        ratio = (self.CO + self.C + self.CO2 + self.CH4 + self.CH3 + self.C2H2 * 2 + self.H2CO + self.HCN + self.NH2CHO + self.C2H6*2 + self.HNC + self.CH3NH2 + self.C2H4*2 + self.CN + 3 * self.C3 +
                3 * self.C3H + 3 * self.C3H2 + self.CH3OH + 3 * self.CH3CCH)/(self.CO + self.O + 2 * self.O2 + 2*self.CO2 + self.NO + self.H2O + self.H2CO + self.NH2CHO + self.CH3OH)
        return ratio

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
    outputs = ['output_points_0.dat', 'output_points_1.dat', 'output_points_2.dat', 'output_points_3.dat']#, 'output_points_4.dat']
    home_dir = '/Users/pooneh/Library/Mobile Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/'
    r = []
    theta = []
    CO = []
    C = []
    CO2 = []
    CH4 = []
    CH3 = []
    O = []
    O2 = []
    OH = []
    NO = []
    H2O = []
    CS = []
    SO = []
    C2H2= []
    H2CO= []
    HCN= []
    NH2CHO= []
    C2H6= []
    HNC= []
    CH3NH2= []
    C2H4= []
    CN= []
    C3= []
    C3H= []
    C3H2= []
    CH3OH= []
    CH3CCH= []
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
            elif line.startswith(' CO '):
                parts_CO = line.split()
                parts_CO = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_CO]
                value_CO = [float(value_CO) for value_CO in parts_CO[1:-1]][-1]
                CO.append(value_CO)
            elif line.startswith(' C '):
                parts_C = line.split()
                parts_C = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_C]
                value_C = [float(value_C) for value_C in parts_C[1:-1]][-1]
                C.append(value_C)
            elif line.startswith(' CO2 '):
                parts_CO2 = line.split()
                parts_CO2 = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_CO2]
                value_CO2 = [float(value_CO2) for value_CO2 in parts_CO2[1:-1]][-1]
                CO2.append(value_CO2)
            elif line.startswith(' CH4 '):
                parts_CH4 = line.split()
                parts_CH4 = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_CH4]
                value_CH4 = [float(value_CH4) for value_CH4 in parts_CH4[1:-1]][-1]
                CH4.append(value_CH4)
            elif line.startswith(' CH3 '):
                parts_CH3 = line.split()
                parts_CH3 = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_CH3]
                value_CH3 = [float(value_CH3) for value_CH3 in parts_CH3[1:-1]][-1]
                CH3.append(value_CH3)
            elif line.startswith(' O '):
                parts_O = line.split()
                parts_O = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_O]
                value_O = [float(value_O) for value_O in parts_O[1:-1]][-1]
                O.append(value_O)
            elif line.startswith(' O2 '):
                parts_O2 = line.split()
                parts_O2 = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_O2]
                value_O2 = [float(value_O2) for value_O2 in parts_O2[1:-1]][-1]
                O2.append(value_O2)
            elif line.startswith(' OH '):
                parts_OH = line.split()
                parts_OH = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_OH]
                value_OH = [float(value_OH) for value_OH in parts_OH[1:-1]][-1]
                OH.append(value_OH)
            elif line.startswith(' NO '):
                parts_NO = line.split()
                parts_NO = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_NO]
                value_NO = [float(value_NO) for value_NO in parts_NO[1:-1]][-1]
                NO.append(value_NO)
            elif line.startswith(' H2O '):
                parts_H2O = line.split()
                parts_H2O = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_H2O]
                value_H2O = [float(value_H2O) for value_H2O in parts_H2O[1:-1]][-1]
                H2O.append(value_H2O)
            elif line.startswith(' CS '):
                parts_CS = line.split()
                parts_CS = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_CS]
                value_CS = [float(value_CS) for value_CS in parts_CS[1:-1]][-1]
                CS.append(value_CS)
            elif line.startswith(' SO '):
                parts_SO = line.split()
                parts_SO = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_SO]
                value_SO = [float(value_SO) for value_SO in parts_SO[1:-1]][-1]
                SO.append(value_SO)


            elif line.startswith(' C2H2 '):
                parts_C2H2 = line.split()
                parts_C2H2 = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_C2H2]
                value_C2H2 = [float(value_C2H2) for value_C2H2 in parts_C2H2[1:-1]][-1]
                C2H2.append(value_C2H2)
            elif line.startswith(' H2CO '):
                parts_H2CO = line.split()
                parts_H2CO = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_H2CO]
                value_H2CO = [float(value_H2CO) for value_H2CO in parts_H2CO[1:-1]][-1]
                H2CO.append(value_H2CO)
            elif line.startswith(' HCN '):
                parts_HCN = line.split()
                parts_HCN= ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_HCN]
                value_HCN = [float(value_HCN) for value_HCN in parts_HCN[1:-1]][-1]
                HCN.append(value_HCN)
            elif line.startswith(' NH2CHO '):
                parts_NH2CHO = line.split()
                parts_NH2CHO = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_NH2CHO]
                value_NH2CHO = [float(value_NH2CHO) for value_NH2CHO in parts_NH2CHO[1:-1]][-1]
                NH2CHO.append(value_NH2CHO)
            elif line.startswith(' C2H6 '):
                parts_C2H6 = line.split()
                parts_C2H6 = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_C2H6]
                value_C2H6 = [float(value_C2H6) for value_C2H6 in parts_C2H6[1:-1]][-1]
                C2H6.append(value_C2H6)
            elif line.startswith(' HNC '):
                parts_HNC = line.split()
                parts_HNC = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_HNC]
                value_HNC = [float(value_HNC) for value_HNC in parts_HNC[1:-1]][-1]
                HNC.append(value_HNC)
            elif line.startswith(' CH3NH2 '):
                parts_CH3NH2 = line.split()
                parts_CH3NH2 = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_CH3NH2]
                value_CH3NH2 = [float(value_CH3NH2) for value_CH3NH2 in parts_CH3NH2[1:-1]][-1]
                CH3NH2.append(value_CH3NH2)
            elif line.startswith(' C2H4 '):
                parts_C2H4 = line.split()
                parts_C2H4 = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_C2H4]
                value_C2H4 = [float(value_C2H4) for value_C2H4 in parts_C2H4[1:-1]][-1]
                C2H4.append(value_C2H4)
            elif line.startswith(' CN '):
                parts_CN = line.split()
                parts_CN = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_CN]
                value_CN = [float(value_CN) for value_CN in parts_CN[1:-1]][-1]
                CN.append(value_CN)
            elif line.startswith(' C3 '):
                parts_C3 = line.split()
                parts_C3 = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_C3]
                value_C3 = [float(value_C3) for value_C3 in parts_C3[1:-1]][-1]
                C3.append(value_C3)
            elif line.startswith(' C3H '):
                parts_C3H = line.split()
                parts_C3H = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_C3H]
                value_C3H = [float(value_C3H) for value_C3H in parts_C3H[1:-1]][-1]
                C3H.append(value_C3H)
            elif line.startswith(' C3H2 '):
                parts_C3H2 = line.split()
                parts_C3H2 = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_C3H2]
                value_C3H2 = [float(value_C3H2) for value_C3H2 in parts_C3H2[1:-1]][-1]
                C3H2.append(value_C3H2)
            elif line.startswith(' CH3OH '):
                parts_CH3OH = line.split()
                parts_CH3OH = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_CH3OH]
                value_CH3OH = [float(value_CH3OH) for value_CH3OH in parts_CH3OH[1:-1]][-1]
                CH3OH.append(value_CH3OH)
            elif line.startswith(' CH3CCH '):
                parts_CH3CCH = line.split()
                parts_CH3CCH = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_CH3CCH]
                value_CH3CCH = [float(value_CH3CCH) for value_CH3CCH in parts_CH3CCH[1:-1]][-1]
                CH3CCH.append(value_CH3CCH)

    CO = np.array(CO)
    C = np.array(C)
    CO2 = np.array(CO2)
    CH4 = np.array(CH4)
    CH3 = np.array(CH3)
    O = np.array(O)
    O2 = np.array(O2)
    OH = np.array(OH)
    NO = np.array(NO)
    H2O = np.array(H2O)
    CS = np.array(CS)
    SO = np.array(SO)
    C2H2= np.array(C2H2)
    H2CO= np.array(H2CO)
    HCN= np.array(HCN)
    NH2CHO= np.array(NH2CHO)
    C2H6= np.array(C2H6)
    HNC= np.array(HNC)
    CH3NH2= np.array(CH3NH2)
    C2H4= np.array(C2H4)
    CN= np.array(CN)
    C3= np.array(C3)
    C3H= np.array(C3H)
    C3H2= np.array(C3H2)
    CH3OH= np.array(CH3OH)
    CH3CCH= np.array(CH3CCH)
    r = np.array(r)
    theta = np.array(theta)
    if len(theta) != len(r):
        print('ERROR the lengths are not the same something is not being found')
        sys.exit()
    if len(OH) != len(CO2):
        print('ERROR the lengths are not the same something is not being found')
        sys.exit()
    if len(r) != len(C):
        print('ERROR the lengths are not the same something is not being found')
        sys.exit()
    source = My_source(C_O_ratio,CO, C, CO2, CH4, CH3, O, O2, OH, NO, H2O, CS, SO, C2H2, H2CO, HCN, NH2CHO, C2H6, HNC, CH3NH2, C2H4, CN, C3, C3H, C3H2, CH3OH, CH3CCH, r, theta, M_star)
    return source
###########################################################################
###########################################################################
def plot_maps_write_CS_SO(M_star, C_O_ratio):
    source = read_output_data(M_star, C_O_ratio)
    new_r, new_z = source.convert_to_cylindrical()
    calc_C_O = source.calc_C_O_ratio()
    CS_SO = source.CS_SO()
    myfile = open('calc_C_O_'+str(C_O_ratio)+'_'+M_star+'.dat','w')
    myfile.write('r_cyl[au]\t z_cyl[au]\t C/O \t CS/SO \n')
    for i in range(len(new_r)):
        myfile.write(str(new_r[i])+'\t'+str(new_z[i])+'\t'+str(calc_C_O[i])+'\t'+str(CS_SO[i])+'\n')
    myfile.close()

    cmap = plt.get_cmap('magma')
    vmin = -3
    vmax = 3

    plt.figure()
    plt.scatter(new_r, new_z, c=np.log10(CS_SO), cmap= cmap, vmin=vmin, vmax = vmax,marker='s',s=80)
    plt.annotate('Initial C/O = '+str(C_O_ratio), xycoords='axes fraction', xy=(0.02,0.9),
                  weight='bold',fontsize=18)
    plt.xlabel('R [au]')
    plt.ylabel('z [au]')
    plt.colorbar(label='CS/SO')
    plt.tight_layout()
    plt.savefig('map_CS_SO_'+str(C_O_ratio)+'_'+str(M_star)+'.png')

    plt.figure()
    plt.scatter(new_r, new_z, c=calc_C_O, cmap= cmap, vmin=0.5, vmax = 1.7,marker='s',s=80)
    plt.annotate('Initial C/O = '+str(C_O_ratio), xycoords='axes fraction', xy=(0.02,0.9),
                  weight='bold',fontsize=18)
    plt.xlabel('R [au]')
    plt.ylabel('z [au]')
    plt.colorbar(label='C/O')
    plt.tight_layout()
    plt.savefig('map_C_O_'+str(C_O_ratio)+'_'+str(M_star)+'.png')

    plt.figure()
    plt.scatter(calc_C_O,CS_SO, c='k', s=20)
    plt.ylabel('CS/SO')
    plt.xlabel('C/O')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(0.1,2)
    plt.tight_layout()
    plt.savefig('C_O_CS_SO_'+str(C_O_ratio)+'_'+str(M_star)+'.png')
###########################################################################
###########################################################################
#########################MAIN##############################################
C_O_ratios = [0.2,0.44, 0.9,1.2]
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
