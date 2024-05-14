import numpy as np
import matplotlib.pylab as plt
from radmc3dPy.analyze import *
import sys
##################################################################
##################################################################
def rerun(first_file,second_file):
    """
    A function that determines how different the two tempertaure files are
    inputs:
        two numbers that are supposed to be in the names of dust_temperature_**number**.dat
    returns:
        True if the difference is too much and teh code needs to be re-run
        False if the difference is small enough to take the current temperature
    """
    a_T    = readData(dtemp=True)
    r_T    = a_T.grid.x[:]
    theta_T = a_T.grid.y[:]
    ###########Read first mid-plane T
    my_T = np.loadtxt('dust_temperature'+str(first_file)+'.dat',skiprows=3)
    hdr = np.loadtxt('dust_temperature'+str(first_file)+'.dat')[0:3]
    if len(my_T) != hdr[1]:
        print('Your data shape is different from what the header says!!!!')
        sys.exit()
    my_T = np.reshape(my_T, [1,1, len(theta_T), len(r_T)])
    my_T = np.swapaxes(my_T, 0, 3)
    my_T = np.swapaxes(my_T, 1, 2)
    my_T_1d_first = my_T[:,-1,0,0]
    ###########
    ###########Read second mid-plane T
    my_T = np.loadtxt('dust_temperature'+str(second_file)+'.dat',skiprows=3)
    hdr = np.loadtxt('dust_temperature'+str(second_file)+'.dat')[0:3]
    if len(my_T) != hdr[1]:
        print('Your data shape is different from what the header says!!!!')
        sys.exit()
    my_T = np.reshape(my_T, [1,1, len(theta_T), len(r_T)])
    my_T = np.swapaxes(my_T, 0, 3)
    my_T = np.swapaxes(my_T, 1, 2)
    my_T_1d_second = my_T[:,-1,0,0]
    ###Find the diff between first and second
    ###
    plt.figure()
    au = 1.496e13
    plt.loglog(r_T/au,my_T_1d_first,label='first')
    plt.loglog(r_T/au,my_T_1d_second, label='sec')
    plt.legend()
    plt.show()
    ###
    diff = abs(my_T_1d_first - my_T_1d_second)/my_T_1d_second * 100
    is_good = all(x < 20 for x in diff)
    if is_good == True:
        re_run = False
    else:
        re_run = True
    return re_run
##################################################################
##################################################################
print(rerun(1,2))
