import numpy as np
import os
############################
##Constants
au =  1.49598e13 #cm
p_mas = 1.67e-24 #gr
ls = 3.8525e33 #erg/s
rs = 6.96e10 # cm solar radius
ss  = 5.6703e-5      # Stefan-Boltzmann const  [erg/cm^2/K^4/s]
ms = 1.98892e33 # mass of sun

#######All functions########
def grid_maker(M_star, L_bol, M_env, T_star, Phot, Nx, Ny, Mdisk, Rdisk, Rc, R_in, R_out,R_full, P_factor):
    """
    A function that creates a dictionary to store all parameters related to each M_star for the grid
    Args:
        M_star: array of stellar mass with X masses
        L_bol: array of luminosities with shape of X by Y where Y is the number of luminosities for each mass (the luminosity can change for each stellar mass)
        Epsilon: To define H = epsilon*R. It is assumed to vary for each luminosity and is an array with shape of Y
        M_env: The envelope mass. Has the same dimension as luminosity
        T_star: Temperature of star (assumed as single value due to small variation)
    Returns:
        A dictionary where the stellar mass is linked to other parameters
    """
    grid_dict = {}
    for i_M in range(len(M_star)):
        rstar = np.sqrt(L_bol[i_M] * ls/(4. * np.pi * T_star**4 * ss)) * 1./rs #This is an array
        grid_dict[M_star[i_M]] = {'L_bol':L_bol[i_M], 'M_env': M_env[i_M], 'r_in':R_in, 'r_out':R_out, 'r_full': R_full, 'T_star':T_star,
                                      'p':P_factor, 'mdisk':Mdisk, 'rdisk':Rdisk, 'r_c':Rc, 'r_star':rstar, 'nx':Nx, 'ny':Ny, 'phot':Phot}
    return grid_dict
########
########
########
def directory_maker(M_star):
    """
    Makes a list of directory names
    Args:
        M_star: array of stellar mass with X masses
    Returns:
        Three lists with directories of stellar mass.
    """
    m_star_dir = []
    for i_M in range(len(M_star)):
        m_star_dir.append('M_star_'+str(M_star[i_M]))
    return m_star_dir
########
########
########
def write_dustsetup(M_star, L_bol, M_env, T_star, Phot, Nx, Ny, Mdisk, Rdisk, Rc, R_in, R_out,R_full, P_factor):
    """
    A function that writes the setup_model.py file to set up the dust in RADMC3D
    Returns:
        None
    """
    m_star_dir = directory_maker(M_star)
    grid_dict = grid_maker(M_star, L_bol, M_env, T_star, Phot, Nx, Ny, Mdisk, Rdisk, Rc, R_in, R_out, R_full, P_factor)
    for i_star in range(len(m_star_dir)):
        if not os.path.exists(m_star_dir[i_star]+'/'):
            os.makedirs(m_star_dir[i_star]+'/')
        f1 = open(m_star_dir[i_star]+'/setup_model.py','w')
        f1.write("import radmc3dPy \n")
        f1.write("import os \n")
        f1.write("radmc3dPy.analyze.writeDefaultParfile('disk_env_model_Pooneh_final') \n")
        f1.write("radmc3dPy.setup.problemSetupDust('disk_env_model_Pooneh_final',nphot='"+str(grid_dict[M_star[i_star]]['phot'])+"',nx='"+grid_dict[M_star[i_star]]['nx']+"',"\
                    "ny='"+str(grid_dict[M_star[i_star]]['ny'])+"',xbound='["+str(grid_dict[M_star[i_star]]['r_in'])+"*au,0.5*au,"+str(grid_dict[M_star[i_star]]['r_out'])+"*au]',"\
                    "r_fullenv='"+str(grid_dict[M_star[i_star]]['r_full'])+"*au',mgas='"+str(grid_dict[M_star[i_star]]['M_env'][0])+"*ms',r_in='"+str(grid_dict[M_star[i_star]]['r_in'])+"*au',"\
                    "prho='"+str(grid_dict[M_star[i_star]]['p'])+"', mstar='"+str(M_star[i_star])+"*ms', rstar='"+str(grid_dict[M_star[i_star]]['r_star'][0])+"*rs',"\
                    "tstar='"+str(grid_dict[M_star[i_star]]['T_star'])+"', binary=False,"\
                    "Rc='"+str(grid_dict[M_star[i_star]]['r_c'])+"*au', rdisk='"+str(grid_dict[M_star[i_star]]['rdisk'])+"*au', mdisk='"+str(grid_dict[M_star[i_star]]['mdisk'])+"*ms',"\
                    "scattering_mode_max='1', Tau='3.0', wbound='[0.0912, 7., 25., 1e4]', nw='[20,100,30]',epsilon = '0.2') \n")
        f1.write("os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/datafiles/dustkappa_silicate.inp .') \n")
        f1.write("os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/datafiles/mcmono_wavelength_micron.inp .')\n")
        f1.write("os.system('radmc3d mctherm')\n")
        f1.write("os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/datafiles/main.f90 .')\n")
        f1.write("os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/datafiles/userdef_module.f90 .')\n")
        f1.write("os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/datafiles/Makefile .')\n")
        f1.write("os.system('make')\n")
        f1.write("os.system('./radmc3d mcmono nphot_mono 1000000')")
########
########
########
################
#Main code
################
##Inputs
t_star = 4000 #K, star effective temperature
l_bol = np.array([[0.5], #M_ClassII  ###Luminosity in L_sun
         [10.], #M_ClassI
         [50.]]) #M_Class0
m_env = np.array([[0.13], #M_ClassII  ##Envelope mass in M_sun
         [1.], #M_ClassI
         [4.]]) #M_Class0
m_star = np.array([0.1,0.5,1.]) #M_sun
mdisk = 0.01 #M_sun
rdisk = 50. #au
r_c = 50. #au.
r_in = 0.4 #au
r_out = 200 #au
r_full = 1e4 #au
p = -1.5
phot = 1000000 #Number of photons
nx = '[100,200]' #grid points in x direction
ny = 200 #grid points in y direction
#########
#Write your setup scripts
write_dustsetup(m_star, l_bol, m_env, t_star, phot, nx, ny, mdisk, rdisk, r_c, r_in, r_out,r_full, p)
