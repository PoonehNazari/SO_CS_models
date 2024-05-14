import radmc3dPy 
import os 
#radmc3dPy.analyze.writeDefaultParfile('disk_env_model_Pooneh_final') 
#radmc3dPy.setup.problemSetupDust('disk_env_model_Pooneh_final',nphot='1000000',nx='[50,150]',ny='100',xbound='[0.4*au,0.5*au,100*au]',r_fullenv='2000.0*au',mgas='0.1*ms',r_in='0.4*au',prho='-1.5', mstar='0.22*ms', rstar='6.26405274185309*rs',tstar='4000', binary=False, ybound='[0.0, pi/2.]',Rc='88.0*au', rdisk='44.0*au', mdisk='0.03*ms',scattering_mode_max='1', Tau='3.0', wbound='[0.0912, 7., 25., 1e4]', nw='[20,100,30]',epsilon = '0.6') 
#os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/datafiles/dustkappa_silicate.inp .') 
#os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/datafiles/mcmono_wavelength_micron.inp .')
#os.system('radmc3d mctherm')
os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/datafiles/main.f90 .')
os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/datafiles/userdef_module.f90 .')
os.system('cp -v /Users/pooneh/Library/Mobile\ Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/physical_model/Models/FUV/datafiles/Makefile .')
os.system('make')
os.system('./radmc3d mcmono nphot_mono 1000000')
