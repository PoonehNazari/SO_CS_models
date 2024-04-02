import radmc3dPy 
import os 
radmc3dPy.analyze.writeDefaultParfile('disk_env_model_Pooneh_final') 
radmc3dPy.setup.problemSetupDust('disk_env_model_Pooneh_final',nphot='1000000',nx='[100,200]',ny='200',xbound='[0.4*au,0.5*au,200*au]',r_fullenv='10000.0*au',mgas='1.0*ms',r_in='0.4*au',prho='-1.5', mstar='0.5*ms', rstar='6.60289134922617*rs',tstar='4000', binary=False,Rc='50.0*au', rdisk='50.0*au', mdisk='0.01*ms',scattering_mode_max='1', Tau='3.0', wbound='[0.0912, 7., 25., 1e4]', nw='[20,100,30]',epsilon = '0.2') 
os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/datafiles/dustkappa_silicate.inp .') 
os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/datafiles/mcmono_wavelength_micron.inp .')
os.system('radmc3d mctherm')os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/datafiles/main.f90 .')
os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/datafiles/userdef_module.f90 .')
os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/datafiles/Makefile .')
os.system('make')
os.system('./radmc3d mcmono nphot_mono 1000000')