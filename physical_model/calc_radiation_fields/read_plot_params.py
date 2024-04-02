import matplotlib
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 26
matplotlib.rcParams['lines.linewidth'] = 3.5
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['legend.numpoints'] = 1
matplotlib.rcParams['legend.handlelength'] = 2
matplotlib.rcParams['legend.fontsize'] = 23
matplotlib.rcParams['lines.markeredgewidth'] = 2
matplotlib.rcParams['lines.markersize'] = 8
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.major.width'] = 2.2
matplotlib.rcParams['xtick.major.size'] = 10
matplotlib.rcParams['xtick.minor.width'] = 2.2
matplotlib.rcParams['xtick.minor.size'] = 6
matplotlib.rcParams['ytick.major.width'] = 2.2
matplotlib.rcParams['ytick.major.size'] = 10
matplotlib.rcParams['ytick.minor.width'] = 2.2
matplotlib.rcParams['ytick.minor.size'] = 6
import matplotlib.pylab as plt
import numpy as np
from radmc3dPy.analyze import *
from radmc3dPy.natconst import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
###########Constants###########
p_mas = 1.67e-24 #gr
###########
###########
###########FUNCS###########
def give_grid_dens_temp_rad(M_star_dir):
    os.chdir('/Users/pooneh/Library/Mobile Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/'+M_star_dir)
    a    = readData(ddens=True,dtemp=True)
    r    = a.grid.x[:]
    theta = a.grid.y[:]
    theta[0]=0.
    theta[-1]=np.pi
    r2,theta2 = np.meshgrid(r/au,theta)
    dens_2d = a.rhodust[:,:,0,0].T * 100./(1.4* p_mas)
    temp_2d = a.dusttemp[:,:,0,0].T
    dens_2d_write = a.rhodust[:,:,0,0] * 100./(1.4* p_mas)
    temp_2d_write = a.dusttemp[:,:,0,0]
    dens_1d = a.rhodust[:,int(len(theta)/2),0,0] * 100./(1.4* p_mas)
    temp_1d = a.dusttemp[:,int(len(theta)/2),0,0]
    dens_1d_rand = a.rhodust[:,int(len(theta)/4),0,0] * 100./(1.4* p_mas)
    temp_1d_rand = a.dusttemp[:,int(len(theta)/4),0,0]

    rad = np.loadtxt('userdef_total_flux.out',skiprows=2)
    hdr = np.loadtxt('userdef_total_flux.out')[0:2]
    if len(rad) != hdr[1]:
        print('Your data shape is different from what the header says!!!!')
        sys.exit()
    rad = np.reshape(rad, [1,1, len(theta), len(r)])
    rad = np.swapaxes(rad, 0, 3)
    rad = np.swapaxes(rad, 1, 2)
    rad_1d = rad[:,int(len(theta)/2),0,0]
    rad_1d_rand = rad[:,int(len(theta)/4),0,0]
    rad_2d = rad[:,:,0,0].T
    rad_2d_write = rad[:,:,0,0]
    rad_2d[np.where(rad_2d==0.0000000000000000)] = 1e-100

    return r/au, theta, r2, theta2,dens_1d, dens_1d_rand, dens_2d, temp_1d, temp_1d_rand, temp_2d, rad_1d, rad_1d_rand,rad_2d, dens_2d_write, temp_2d_write, rad_2d_write
###########
###########
###########
def plot_dens_temp_rad(M_star_dir,xlim_2d, ylim_2d,xlim_1d,ylim_1d_rad):
    cmap_temp = plt.get_cmap('hot')
    cmap_temp.set_bad(color = 'k')
    cmap_dens = plt.get_cmap('jet')
    cmap_dens.set_bad(color = 'k')
    cmap_rad = plt.get_cmap('tab20b')
    cmap_rad.set_bad(color = 'k')
    row = 2
    col = 3
    fig,ax = plt.subplots(row,col,figsize=(24,12))
    r, theta, r2, theta2,dens_1d, dens_1d_rand, dens_2d, temp_1d, temp_1d_rand, temp_2d, rad_1d, rad_1d_rand,rad_2d, dens_2d_write, temp_2d_write, rad_2d_write = give_grid_dens_temp_rad(M_star_dir)

    dens_2d_masked = np.ma.masked_where(dens_2d<1.05e3,dens_2d)
    temp_2d_masked = np.ma.masked_where(dens_2d<1.05e3,temp_2d)
    rad_2d_masked = np.ma.masked_where(dens_2d<1.05e3,rad_2d)

    c = ax[0,0].pcolormesh(r2 * np.sin(theta2), r2 * np.cos(theta2),
                           np.log10(dens_2d_masked), cmap=cmap_dens,vmin=6.5,vmax=10)
    ax[1,0].plot(r, dens_1d, c='k')
    ax[1,0].plot(r, dens_1d_rand, c='tab:red')
    #ax[1,0].set_xscale('log')
    ax[1,0].set_yscale('log')
    divider = make_axes_locatable(ax[0,0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = fig.colorbar(c,cax=cax)
    cb.set_label('$\log{n_{\mathdefault{H}}}$ [cm$^{-3}$]', rotation=270.,labelpad = 20.)

    c = ax[0,1].pcolormesh(r2 * np.sin(theta2), r2 * np.cos(theta2),
                       temp_2d_masked, cmap=cmap_temp,vmin=20,vmax=200)
    ax[1,1].plot(r, temp_1d, c='k')
    ax[1,1].plot(r, temp_1d_rand, c='tab:red')
    #ax[1,1].set_xscale('log')
    ax[1,1].set_yscale('log')
    divider = make_axes_locatable(ax[0,1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = fig.colorbar(c,cax=cax)
    cb.set_label('$T$ [K]', rotation=270.,labelpad = 20.)

    c = ax[0,2].pcolormesh(r2 * np.sin(theta2), r2 * np.cos(theta2),
                           np.log10(rad_2d_masked), cmap=cmap_rad,vmin=-4,vmax=6)
    ax[1,2].plot(r, rad_1d, c='k')
    ax[1,2].plot(r, rad_1d_rand, c='tab:red')
    # ax[1,2].set_xscale('log')
    ax[1,2].set_yscale('log')
    divider = make_axes_locatable(ax[0,2])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = fig.colorbar(c,cax=cax)
    cb.set_label('$F_{\mathdefault{UV}}$ [G0]', rotation=270.,labelpad = 20.)

    ax[0,0].set_xlim(0,xlim_2d)
    ax[0,0].set_ylim(0,ylim_2d)
    ax[0,1].set_xlim(0,xlim_2d)
    ax[0,1].set_ylim(0,ylim_2d)
    ax[0,2].set_xlim(0,xlim_2d)
    ax[0,2].set_ylim(0,ylim_2d)

    ax[1,0].set_xlim(0,xlim_1d)
    ax[1,1].set_xlim(0,xlim_1d)
    ax[1,2].set_xlim(0,xlim_1d)

    ax[1,2].set_ylim(1e-4,ylim_1d_rad)


    ax[0,0].set_xlabel('R [au]')
    ax[0,1].set_xlabel('R [au]')
    ax[0,2].set_xlabel('R [au]')

    ax[0,0].set_ylabel('z [au]')
    ax[0,1].set_ylabel('z [au]')
    ax[0,2].set_ylabel('z [au]')

    ax[1,0].set_xlabel('R [au]')
    ax[1,1].set_xlabel('R [au]')
    ax[1,2].set_xlabel('R [au]')

    ax[1,0].set_ylabel('$\log{n_{\mathdefault{H}}}$ [cm$^{-3}$]')
    ax[1,1].set_ylabel('$T$ [K]')
    ax[1,2].set_ylabel('$F_{\mathdefault{UV}}$ [G0]')
    plt.tight_layout()
    plt.savefig(M_star_dir+'_all.png')
###########
###########
###########
######################Main######################
################################################
dirs = ['M_star_0.1','M_star_0.5','M_star_1.0']
for i in range(len(dirs)):
    plot_dens_temp_rad(dirs[i],200, 200,200,1e6)
