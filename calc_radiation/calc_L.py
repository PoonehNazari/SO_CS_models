import numpy as np
######################################################################
au =  1.49598e13 #cm
ls = 3.8525e33 #erg/s
ss  = 5.6703e-5      # Stefan-Boltzmann const  [erg/cm^2/K^4/s]
h_planck = 6.62607015e-27 #cgs
k_B = 1.380649e-16 #cgs
c_light = 2.99792458e10 #cm/s
rs = 6.96e10 #cm
######################################################################
######################################################################
def find_L_UV(T_star):
    lambda_microns = np.linspace(9.120000000e-02, 2.066000000e-01, 500) #lambda in microns
    lambda_cm = lambda_microns[::-1] * 1e-4 # in cm
    freq = c_light/lambda_cm #Hz from low to high
    #Integrate to get the flux and get rid of 1/Hz, multiply by 4pi to get rid of 1/Str
    B = 0.
    for ind_i in range(len(freq)-1):
        B_i = 2. * freq[ind_i]**3 * h_planck/c_light**2 * 1./(np.exp(h_planck*freq[ind_i]/(k_B * T_star)) - 1.) * (freq[ind_i+1] - freq[ind_i])
        B = B + B_i

    B_i = 2. * freq[ind_i+1]**3 * h_planck/c_light**2 * 1./(np.exp(h_planck*freq[ind_i+1]/(k_B * T_star)) - 1.) * (freq[ind_i+1] - freq[ind_i])
    B = B + B_i #in units of erg/s/cm^2/Sr
    B = 4. * np.pi * B #in units of erg/s/cm^2
    #Now convert to luminosity for your star with various R_star (i.e., different luminosity)
    L_star = np.array([0.5,10,50.]) #in solar L
    R_star = np.sqrt(L_star * ls/(4. * np.pi * T_star**4 * ss)) #cm
    L_UV = 4. * np.pi * R_star**2 * B #in erg/s
    print('L_UV (erg/s)= ',L_UV)
    lambda_UV = 1.5e-5 #cm
    E_UV = c_light/lambda_UV * h_planck
    L_UV = L_UV/E_UV #in ph/s
    print('L_UV (ph/s)=',L_UV)
######################################################################
######################################################################
############################Main######################################
T_star = 4000. #K
find_L_UV(T_star)
