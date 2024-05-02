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
import re
###########################################################################
#################################FUNCS/CLASSES#####################################
def process_word_C(word):
    # Regular expression to find 'C' optionally followed by a digit
    pattern = r'C(\d?)'
    matches = re.findall(pattern, word)

    # Initialize sum
    total = 0

    # Process each match
    for match in matches:
        if match.isdigit():  # If there's a digit, add it as an integer
            total += int(match)
        else:  # If there's no digit, add 1 (for 'C' with no following digit)
            total += 1

    return total
###########################################################################
###########################################################################
def process_word_O(word):
    # Regular expression to find 'C' optionally followed by a digit
    pattern = r'O(\d?)'
    matches = re.findall(pattern, word)

    # Initialize sum
    total = 0

    # Process each match
    for match in matches:
        if match.isdigit():  # If there's a digit, add it as an integer
            total += int(match)
        else:  # If there's no digit, add 1 (for 'C' with no following digit)
            total += 1

    return total
###########################################################################
###########################################################################
class My_source:
    def __init__(self, C_O, final_C, final_O , CS, SO, HNCO, CH3CN, r, theta, Mstar):
        self.C_O = C_O #input C/O
        self.final_C = final_C
        self.final_O = final_O
        self.CS = CS
        self.SO = SO
        self.HNCO = HNCO
        self.CH3CN = CH3CN
        self.r = r #array of radii
        self.theta = theta #array of theta
        self.Mstar = Mstar #the stellar mass

    def calc_C_O_ratio(self):
        ratio = self.final_C/self.final_O
        return ratio

    def CS_SO(self):
        ratio = self.CS/self.SO
        return ratio

    def CH3CN_HNCO(self):
        ratio = self.CH3CN/self.HNCO
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
    sections = []
    CS = []
    SO = []
    HNCO = []
    CH3CN = []
    for i_output in range(len(outputs)):
        output = open(home_dir+'outputs/'+M_star+'/C_O_'+str(C_O_ratio)+'/'+outputs[i_output])
        lines_output = output.readlines()
        output.close()
        current_section = []
        collect = False
        for line in lines_output:
            if line.startswith(' RADIUS ='):
                current_radius = float(line.split('=')[1].strip().split()[0])
                r.append(current_radius)
                collect = False
            elif line.startswith(' HEIGHT ='):
                current_theta = float(line.split('=')[1].strip().split()[0])
                theta.append(current_theta)
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
            elif line.startswith(' HNCO '):
                parts_HNCO = line.split()
                parts_HNCO = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_HNCO]
                value_HNCO = [float(value_HNCO) for value_HNCO in parts_HNCO[1:-1]][-1]
                HNCO.append(value_HNCO)
            elif line.startswith(' CH3CN '):
                parts_CH3CN = line.split()
                parts_CH3CN = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in parts_CH3CN]
                value_CH3CN = [float(value_CH3CN) for value_CH3CN in parts_CH3CN[1:-1]][-1]
                CH3CN.append(value_CH3CN)
            elif line.startswith(' *'):
                # If current_section is not empty, add it to sections
                if current_section:
                    sections.append(current_section)
                    current_section = []
                # Set flag to true to start collecting lines after this
                collect = True
                continue  # Skip adding the 'RADIUS' line to the section
            # Add line to the current section
            if collect:
                current_section.append(line.strip())
            # Add the last section if not empty
        if current_section:
            sections.append(current_section)

    final_C = []
    final_O = []
    for sect in sections:
        carbons =[]
        oxygens = []
        for new_line in sect:
            new_line = new_line.split()
            new_line = ['1E-20' if ('+1' in element) and ('E+1' not in element) else element for element in new_line]
            new_line = ['1E-20' if ('+2' in element) and ('E+2' not in element) else element for element in new_line]
            if new_line != []:
                value = float(new_line[-1])
                if not new_line[0].startswith('G'):
                    total_C = process_word_C(new_line[0])
                    carbons.append(float(total_C) * value)
                    total_O = process_word_O(new_line[0])
                    oxygens.append(float(total_O) * value)
        final_O.append(sum(oxygens))
        #print(sum(oxygens))
        final_C.append(sum(carbons))

    final_C = np.array(final_C)
    final_O = np.array(final_O)
    r = np.array(r)
    theta = np.array(theta)
    CS = np.array(CS)
    SO = np.array(SO)
    HNCO = np.array(HNCO)
    CH3CN = np.array(CH3CN)
    if len(theta) != len(r):
        print('ERROR the lengths are not the same something is not being found')
        sys.exit()
    if len(theta) != len(final_O):
        print('ERROR the lengths are not the same something is not being found')
        sys.exit()
    source = My_source(C_O_ratio, final_C, final_O , CS, SO,HNCO, CH3CN, r, theta, M_star)
    return source
###########################################################################
###########################################################################
def plot_maps_write_CS_SO(M_star, C_O_ratio):
    source = read_output_data(M_star, C_O_ratio)
    new_r, new_z = source.convert_to_cylindrical()
    calc_C_O = source.calc_C_O_ratio()
    CS_SO = source.CS_SO()
    CH3CN_HNCO = source.CH3CN_HNCO()
    myfile = open('calc_C_O_'+str(C_O_ratio)+'_'+M_star+'.dat','w')
    myfile.write('r_cyl[au]\t z_cyl[au]\t C/O \t CS/SO \t CH3CN/HNCO \n')
    for i in range(len(new_r)):
        myfile.write(str(new_r[i])+'\t'+str(new_z[i])+'\t'+str(calc_C_O[i])+'\t'+str(CS_SO[i])+'\t'+str(CH3CN_HNCO[i])+'\n')
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
    plt.close()

    plt.figure()
    plt.scatter(new_r, new_z, c=calc_C_O, cmap= cmap, vmin=0.01, vmax = 1.3,marker='s',s=80)
    plt.annotate('Initial C/O = '+str(C_O_ratio), xycoords='axes fraction', xy=(0.02,0.9),
                  weight='bold',fontsize=18)
    plt.xlabel('R [au]')
    plt.ylabel('z [au]')
    plt.colorbar(label='C/O')
    plt.tight_layout()
    plt.savefig('map_C_O_'+str(C_O_ratio)+'_'+str(M_star)+'.png')
    plt.close()

    plt.figure()
    plt.scatter(calc_C_O,CS_SO, c='k', s=20)
    plt.ylabel('CS/SO')
    plt.xlabel('C/O')
    plt.yscale('log')
    #plt.xscale('log')
    # plt.ylim(0.01,10)
    plt.xlim(0.1,2)
    plt.tight_layout()
    plt.savefig('C_O_CS_SO_'+str(C_O_ratio)+'_'+str(M_star)+'.png')
    plt.close()
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
