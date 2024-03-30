import radmc3dPy 
import os 
radmc3dPy.analyze.writeDefaultParfile('disk_env_model_Pooneh_final') 
radmc3dPy.setup.problemSetupDust('disk_env_model_Pooneh_final',nphot='1000000',nx='[300,700]',ny='400',xbound='[0.4*au,0.5*au,10000.0*au]',r_fullenv='10000.0*au',mgas='0.13*ms',r_in='0.4*au',prho='-1.5', mstar='0.1*ms', rstar='1.476451390491502*rs',tstar='4000', binary=False,Rc='50.0*au', rdisk='50.0*au', mdisk='0.01*ms',scattering_mode_max='1', Tau='3.0') 
os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/gap_opening/datafiles/dustkappa_silicate.inp .') 
os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/gap_opening/datafiles/mcmono_wavelength_micron.inp .')
#os.system('radmc3d mctherm')
