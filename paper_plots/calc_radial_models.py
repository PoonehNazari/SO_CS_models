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
#################################FUNCS#####################################
def plot_radial(M_star, C_O_ratios):
    fig, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    for i_ratios in range(len(C_O_ratios)):
        data = np.loadtxt('calc_C_O_'+str(C_O_ratios[i_ratios])+'_'+M_star+'.dat',skiprows=1)
        r = data[:,0]
        z = data[:,1]
        C_O = data[:,2]
        CS_SO = data[:,3]
        CH3CN_HNCO = data[:,4]

        bound = np.linspace(0,200, 100) #in au
        binned_r = []
        binned_C_O =[]
        binned_CS_SO = []
        binned_CH3CN_HNCO = []
        for i_bound in range(len(bound)-1):
            indicies = np.where((r>bound[i_bound]) & (r<bound[i_bound+1]))[0]
            binned_CS_SO.append(np.nanmean(CS_SO[indicies]))
            binned_CH3CN_HNCO.append(np.nanmean(CH3CN_HNCO[indicies]))
            binned_C_O.append(np.nanmean(C_O[indicies]))
            binned_r.append(np.nanmean(r[indicies]))


        for i_bins in range(len(binned_r)):
            if (np.isnan(binned_CS_SO[i_bins])) & (not np.isnan(binned_CS_SO[i_bins-1])):
                binned_CS_SO[i_bins] = binned_CS_SO[i_bins-1]
            if (np.isnan(binned_CS_SO[i_bins])) & (np.isnan(binned_CS_SO[i_bins-1])):
                binned_CS_SO[i_bins] = binned_CS_SO[i_bins+1]

        for i_bins in range(len(binned_r)):
            if (np.isnan(binned_CH3CN_HNCO[i_bins])) & (not np.isnan(binned_CH3CN_HNCO[i_bins-1])):
                binned_CH3CN_HNCO[i_bins] = binned_CH3CN_HNCO[i_bins-1]
            if (np.isnan(binned_CH3CN_HNCO[i_bins])) & (np.isnan(binned_CH3CN_HNCO[i_bins-1])):
                binned_CH3CN_HNCO[i_bins] = binned_CH3CN_HNCO[i_bins+1]

        ax1.plot(binned_r,binned_CS_SO, label=str(C_O_ratios[i_ratios]))
        ax2.plot(binned_r,binned_CH3CN_HNCO, label=str(C_O_ratios[i_ratios]))

    ax1.set_ylabel('CS/SO', color='k')
    ax1.tick_params(axis='y', labelcolor='k')
    ax2.set_ylabel('CH3CN/HNCO', color='k')
    ax1.legend()
    ax2.legend()


    ax1.set_xlabel('R [au]')
    ax2.set_xlabel('R [au]')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax2.set_ylim(0.1,1e7)
    ax1.set_xlim(10,200)
    ax2.set_xlim(10,200)
    fig.tight_layout()
    fig.savefig('Radial_profiles_CS_SO.png')
    plt.close(fig)

    fig2.tight_layout()
    fig2.savefig('Radial_profiles_CH3CN_HNCO.png')
    plt.close(fig2)
###########################################################################
###########################################################################
#########################MAIN##############################################
C_O_ratios = [0.2,0.44, 0.9]#,1.2]
M_stars = ['M_star_0.5']
for i_star in range(len(M_stars)):
    plot_radial(M_stars[i_star], C_O_ratios)
