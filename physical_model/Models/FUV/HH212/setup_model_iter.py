import radmc3dPy
import os
import sys
radmc3dPy.analyze.writeDefaultParfile('disk_env_model_Pooneh_final_iter_scale_height')
radmc3dPy.setup.problemSetupDust('disk_env_model_Pooneh_final_iter_scale_height',nphot='1000000',nx='[50,150]',ny='100',xbound='[0.4*au,0.5*au,100*au]',r_fullenv='2000.0*au',mgas='0.1*ms',r_in='0.4*au',prho='-1.5', mstar='0.22*ms', rstar='6.26405274185309*rs',tstar='4000', binary=False, ybound='[0.0, pi/2.]',Rc='88.0*au', rdisk='44.0*au', mdisk='0.03*ms',scattering_mode_max='1', Tau='3.0', wbound='[0.0912, 7., 25., 1e4]', nw='[20,100,30]', final_index = '1')
os.system('radmc3d mctherm')

