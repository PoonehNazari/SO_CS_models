import numpy as np
import sys
import read_plot_params as rd_mods

r, theta, r2, theta2,dens_1d, dens_1d_rand, dens_2d, temp_1d, temp_1d_rand, temp_2d, rad_1d, rad_1d_rand,rad_2d, dens_2d_write, temp_2d_write, rad_2d_write = give_grid_dens_temp_rad(M_star_dir)
myfile = open('/Users/pooneh/Library/Mobile Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/inputs/'+filename, 'w')

for i_r in range(len(r)):
    for i_theta in range(len(theta)):
        myfile.write(str(r[i_r])+'\t'+str(theta[i_theta])+'\t'+str(dens_2d_write[i_r,i_theta,0,0])+'\t'+str(temp_2d_write[i_r,i_theta,0,0])+'\t'
                     +str(temp_2d_write[i_r,i_theta,0,0])+'\t'+str(rad_2d_write[i_r,i_theta,0,0])+str(0.00000e+00)+'\t'+str(5.00000E-17)+'\t'+str(0.00000e+00))
