import numpy as np
import sys
import read_plot_params as rd_mods


def write_input(M_star_dir):
    r, theta, r2, theta2,dens_1d, dens_1d_rand, dens_2d, temp_1d, temp_1d_rand, temp_2d, rad_1d, rad_1d_rand,rad_2d, dens_2d_write, temp_2d_write, rad_2d_write = rd_mods.give_grid_dens_temp_rad(M_star_dir)
    bad_r_array = np.where(dens_2d_write<1.05e3)[0]
    bad_theta_array = np.where(dens_2d_write<1.05e3)[1]
    f_list = []
    for i_bad in range(len(bad_r_array)):
        f_list.append([bad_r_array[i_bad],bad_theta_array[i_bad]])
    myfile = open('/Users/pooneh/Library/Mobile Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/inputs/'+M_star_dir+'/input_points_'+str(n_i)+'.txt', 'w')
    for i_r in range(len(r)):
        for i_theta in range(len(theta)):
            if [i_r,i_theta] in f_list:
                print('bad index '+str(i_r)+'_'+str(i_theta))
            else:
                myfile.write(f'{r[i_r]:.6f}'+'\t'+f'{theta[i_theta]:.6f}'+'\t'+f'{dens_2d_write[i_r,i_theta]:.6f}'+'\t'+f'{temp_2d_write[i_r,i_theta]:.6f}'+'\t'
                            +f'{temp_2d_write[i_r,i_theta]:.6f}'+'\t'+f'{rad_2d_write[i_r,i_theta]:.6f}'+'\t'+f'{0.0:.6f}'+'\t'+f'{5E-17:.6f}'+'\t'+f'{0.0:.6f}'+'\n')

    myfile.close()

dirs = ['M_star_0.1','M_star_0.5','M_star_1.0']
for i_star in range(len(dirs)):
    write_input(dirs[i_star])
