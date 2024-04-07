import numpy as np
import os
import shutil
###########################################################################
###########################################################################
######################################FUNCS################################
def initiate_grid(M_star, C_O, sub_grids):
    """
    A function to setup the directories and the file_parameters.txt file before compiling and running the chemical code
    inputs:
        M_star: a directory with the various star masses
        C_O: the C/O ratios that we want to run teh model for
        sub_grids: the division of each input file into 3000 grid points
    """
    home_dir = '/Users/pooneh/Library/Mobile Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/'
    for i_M in range(len(M_star)):
        for i_C_O in range(len(C_O)):
            if not os.path.exists(home_dir+M_star[i_M]+'/'+C_O[i_C_O]+'/'):
                os.makedirs(home_dir+M_star[i_M]+'/'+C_O[i_C_O]+'/')
            for i_grid in range(len(sub_grids)):
                if os.path.exists(home_dir+M_star[i_M]+'/'+C_O[i_C_O]+'/'+sub_grids[i_grid]+'/'):
                    shutil.rmtree(home_dir+M_star[i_M]+'/'+C_O[i_C_O]+'/'+sub_grids[i_grid]+'/')
                    os.makedirs(home_dir+M_star[i_M]+'/'+C_O[i_C_O]+'/'+sub_grids[i_grid]+'/')
                else:
                    os.makedirs(home_dir+M_star[i_M]+'/'+C_O[i_C_O]+'/'+sub_grids[i_grid]+'/')
                shutil.copytree(home_dir+'core_chemical_model/',
                                home_dir+M_star[i_M]+'/'+C_O[i_C_O]+'/'+sub_grids[i_grid]+'/core_chemical_model/')
                os.chdir(home_dir+M_star[i_M]+'/'+C_O[i_C_O]+'/'+sub_grids[i_grid]+'/core_chemical_model/')
                f1 = open('file_parameters.txt','w')
                f1.write('! FILE PARAMETERS FOR PROTOPLANETARY DISK MODEL \n')
                f1.write('>>>>>>> Physical conditions input file \n')
                f1.write('./../../../../inputs/'+M_star[i_M]+'/input_points_'+str(i_grid)+'.txt \n')
                f1.write('>>>>>>> UV radiation field input files \n')
                f1.write("'DEFAULT' \n")
                f1.write('>>>>>>> X-ray radiation field input \n')
                f1.write("'DEFAULT' \n")
                f1.write('>>>>>>> Reaction file \n')
                f1.write('./../../../../rates_abundances/rate12_full_v2.rates \n')
                f1.write('>>>>>>> Species file \n')
                f1.write('./../../../../rates_abundances/rate12_init_abunds_'+C_O[i_C_O]+'.specs \n')
                f1.write('>>>>>>> Binding energies file \n')
                f1.write('./../../../../rates_abundances/rate12_binding_v2.dat \n')
                f1.write('>>>>>>> Grain parameters file \n')
                f1.write('./grain_parameters.txt \n')
                f1.write('>>>>>>> Radiation field parameters file \n')
                f1.write('./../../../radiation_parameters.txt \n')
                f1.write('>>>>>>> Reaction switch file \n')
                f1.write('./reaction_parameters.txt \n')
                f1.write('>>>>>>> Molecular abundances output file \n')
                f1.write('./../../../../outputs/'+M_star[i_M]+'/'+C_O[i_C_O]+'/output_points_'+str(i_grid)+'.dat \n')
                f1.write('>>>>>>> Reaction rates output file \n')
                f1.write('./../../../../outputs/'+M_star[i_M]+'/'+C_O[i_C_O]+'/output_points_'+str(i_grid)+'_rates.dat \n')
                f1.write('>>>>>>> Analyse output file \n')
                f1.write('./../../../../outputs/'+M_star[i_M]+'/'+C_O[i_C_O]+'/output_points_'+str(i_grid)+'_analyse.dat \n')
                f1.write('>>>>>>> First output time \n')
                f1.write('5.0E+00 \n')
                f1.write('>>>>>>> Final output time \n')
                f1.write('5.0E+04 \n')
                f1.write('>>>>>>> LOG_10(Ti+1) - LOG_10(Ti): time steps \n')
                f1.write('2.0 \n')
                f1.write('>>>>>>> Time step for analyse subroutine (set to 0 if not called) \n')
                f1.write('62 \n')
                f1.write('>>>>>>>')
                f1.close()
###########################################################################
###########################################################################
##############################Main#########################################
M_star = ['M_star_0.1','M_star_0.5','M_star_1.0']
C_O = ['C_O_0.2','C_O_0.44','C_O_0.9','C_O_1.2','C_O_1.5']
sub_grids = ['grid_0','grid_1','grid_2','grid_3','grid_4']
initiate_grid(M_star, C_O, sub_grids)
