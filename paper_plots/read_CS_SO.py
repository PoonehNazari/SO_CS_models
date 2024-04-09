import numpy as np
import matplotlib.pylab as plt
import os

class Source:
    def _init_(self, C_O, CS, SO, r, theta, Mstar):
        self.C_O = C_O #input C/O
        self.CS = CS #final CS abundances
        self.SO = SO #final SO abundances
        self.r = r #array of radii
        self.theta = theta #array of theta
        self.Mstar = Mstar #the stellar mass

    def CS_SO(CS, SO):
        ratio = CS/SO
        return ratio

    def convert_to_cylindrical(r, theta):
        new_r = r * np.sin(theta)
        new_z = r * np.cos(theta)
        return new_r, new_z


inputs = ['input_points_0.txt','input_points_1.txt','input_points_2.txt','input_points_3.txt','input_points_4.txt']
home_dir = '/Users/pooneh/Library/Mobile Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models'
C_O_ratios = [0.2,0.44]
M_stars = ['M_star_0.5']
all_sources = []
n = 0
for i_star in range(len(M_stars)):
    for i_C_O in range(len(C_O_ratios)):
        r = []
        theta = []
        CS = []
        SO = []
        for i_input in range(len(inputs)):
            input_data = np.loadtxt(home_dir+'inputs/M_star_0.5/'+inputs[i_input])
            r_input = input_data[:,0]
            theta_input = input_data[:,1]
            for i_r in range(len(r)):
                os.system('./time_evolution.pl ../outputs/M_star_0.5/C_O_'+str(C_O_ratios[i_C_O])+'/output_points_'+str(0)+'.dat overwriting_file '\
                        +str(format(r_input[i_r], '.3E'))+' '+str(format(theta_input[i_r], '.3E'))+' CS SO')
                data  = np.loadtxt('overwriting_file.dat', skiprows=8)
                CS_i = data[-1,3]
                SO_i = data[-1,4]
                r.append(r_input[i_r])
                theta.append(theta_input[i_r])
                CS.append(CS_i)
                SO.append(SO_i)

        my_source = Source(C_O_ratios[i_C_O], CS, SO, r, theta, M_stars[i_star])
        all_sources.append(my_source)
        n = n+ 1
