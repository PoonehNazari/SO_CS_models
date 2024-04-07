      SUBROUTINE MAIN(DENSMASS,TEMPGAS,TEMPDUST,
     1   GFUV,ZETACR,ZETAXR,FXR,DENS,R,Z,VISEXT,
     2   UVANG,UVFLUX,NUV,XREV,XRFLUX,NXR,PH,XR,
     3   REACFILE,SPECFILE,BINDFILE,GRAINFILE,RADFILE,SWITCHFILE,
     4   TSTART,TFINAL,NANA,DELTA,NREAC,
     5   GISM,TIME,ABUN,NTIME,SPEC,NTOT)   
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Declaration of variables
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      IMPLICIT NONE
      
      INTEGER I,J

C I,J = counters  

       INTEGER D,E,F,G,H            
       PARAMETER (D=1600)
       PARAMETER (E=50)
       PARAMETER (F=1000)
       PARAMETER (G=1000)
       PARAMETER (H=10000)

C D = max number of grid points
C E = max number of UV wavelength/X-ray energy grids
C F = max number of time steps
C G = max number of species
C H = max number of reactions
 
      DOUBLE PRECISION R,Z
      DOUBLE PRECISION DENSMASS,TEMPGAS,TEMPDUST
      DOUBLE PRECISION GFUV,ZETACR,ZETAXR,FXR
      DOUBLE PRECISION DENSNUC,DENS,M_MU
      DOUBLE PRECISION X_H2,X_HE,Y_H2,Y_HE,RATIO
      
C R = disk radius (AU)       
C Z = disk height (AU)
C DENSMASS = gas mass density (g cm-3)
C TEMPGAS =  gas temperature (K)
C TEMPDUST = dust temperature (K)
C GFUV = wavelength-integrated UV field (relative to ISRF) 
C ZETACR = cosmic-ray ionisation rate (s-1)
C ZETAXR = x-ray ionisation rates (s-1)
C FXR = energy-integrated x-ray flux (erg cm-2 s-1) 
C DENSNUC = H nuclei number density (cm-3)
C DENS = gas number density (cm-3)      
C M_MU = mean molecular mass (g)
C X_H2 = fractional abundance of H2 (relative to H nuclei density)
C X_HE = fractional abundance of He (relative to H nuclei density)
C Y_H2 = fractional abundance of H2 (relative to number density)
C Y_HE = fractional abundance of He (relative to number density)
C RATIO = ratio of total number density to H nuclei density

      INTEGER N,NUV,NXR,PH,XR
      DOUBLE PRECISION UVANG(E),UVFLUX(E)
      DOUBLE PRECISION XREV(E),XRFLUX(E)
      DOUBLE PRECISION PHRATE,XRRATE

C N = grid point for UV and X-ray fields
C NUV = number of UV radiation field grid points
C NXR = number of X-ray radiation field grid points
C PH = switch for spectrum-dependent photochemistry
C XR = switch for spectrum-dependent X-ray ionisation rates
C UVANG = UV wavelength grid (A)
C UVGRID = wavelength-dependent UV radiation field (photons cm-2 s-1 A-1)
C XREV = X-ray energy grid (eV) 
C XRGRID = energy-dependent X-ray radiation field (photons cm-2 s-1 eV-1) 
C PHRATE = photorate (s-1)
C XRRATE = x-ray ionisation rate (s-1)

       CHARACTER*90 REACFILE,SPECFILE,BINDFILE
       CHARACTER*90 GRAINFILE,RADFILE,SWITCHFILE

C REACFILE = reaction rates filename
C SPECFILE = species filename
C BINDFILE = binding energy filename
C GRAINFILE = grain parameters filename
C RADFILE = radiation field parameters filename
C SWITCHFILE = switches for chemistry (on/off)
                    
      INTEGER NREAC
      INTEGER LOWTEMP(H),UPTEMP(H),RTYPE(H)
      DOUBLE PRECISION ALPHA(H),BETA(H),GAMMA(H)
      DOUBLE PRECISION SCALE(H)
      DOUBLE PRECISION K(H)
      CHARACTER*10 RE1(H),RE2(H),RE3(H)
      CHARACTER*10 P1(H),P2(H),P3(H),P4(H),P5(H)
      CHARACTER*12 PRFILE(H)

C NREAC = number of reactions
C LOW/UPTEMP = lower/upper temperature bound for rate coefficient (K)
C RTYPE = reaction type (see README file)
C ALPHA,BETA,GAMMA = parameters for reaction rate coefficient (see README file)
C SCALE = scaling factor for photorate
C K = reaction rate coefficients (see README file)
C REx = reactants
C Px = products
C PRFILE = photoreaction cross section file
    
      INTEGER NTIME,NSPEC,NICE
      DOUBLE PRECISION Y0(G),Y(G)
      DOUBLE PRECISION TIME(F),ABUN(G,F)
      DOUBLE PRECISION MASS(G),BIND(G),YIELD(G)
      CHARACTER*10 SPEC(G),SICE(G)

C NTIME = number of time steps
C NSPEC = total number of species
C NICE = number of ice species
C Y0 = initial fractional abundances (relative to Hnuc density)
C Y = time-dependent abundances (relative to total density)
C ABUN = time and species-dependent abundances (relative to total density)
C MASS = species mass (amu)
C BIND = species binding energy (ice only; K)
C YIELD = photodesorption yield (ice only; molecules photon-1)
C SPEC = species names
C SICE = ice species names (for binding energies and grain reactions)

      INTEGER NCON,NELEM,NTOT
      PARAMETER(NELEM=10)
      DOUBLE PRECISION X(1),TOTAL,FRAC(NELEM)
      CHARACTER*2 ELEM(NELEM)

C NCON = number of conserved species (electrons only)
C NELEM = number of elements
C NTOT = total number of species
C X = initalise electron abundance
C TOTAL = total charge (conserved, hence, 0)  
C FRAC = fractional abundance of elements
C ELEM = element names (important for X-ray ionisation rates)    
           
      INTEGER NSWITCH,SWITCH(21)
      
C NSWITCH = number of switches for chemistry
C SWITCH = reaction types to switch off
          
      INTEGER NANA
      DOUBLE PRECISION TSTART,TFINAL,DELTA,T0
      DOUBLE PRECISION TBEGIN,TEND

C NANA = time step for analyse subroutine
C TSTART = first time step for integration (yr or sec)
C TFINAL = final time step (yr or sec)
C T0 = time counter (reset to 0)

      DOUBLE PRECISION GISM,HABING,STARUV,GEXT
      DOUBLE PRECISION LYMAN,LYWAV
      DOUBLE PRECISION ZISM,CRPHOT,UVPHOT,XRPHOT
      DOUBLE PRECISION UVEXT,VISEXT,SHIELD
      DOUBLE PRECISION N_H2,N_CO,N_N2

C GISM = interstellar radiation field (erg cm-2 s-1)
C HABING = habing field (photons cm-2 s-1)
C ZISM = cosmic-ray ionisation rate (s-1)
C CRPHOT = cosmic-ray-induced photon flux (photons cm-2 s-1)
C STARUV = total stellar FUV luminosity (erg s-1)
C GEXT = external radiation field (ISM field or external star)
C LYMAN = total Lyman-alpha luminosity (erg s-1)
C LYWAV =  Lyman-alpha wavelength (A)
C UVPHOT = total UV photon flux (photons cm-2 s-1)
C XRPHOT = total X-ray photon flux (photons cm-2 s-1)
C UVEXT = effective UV extinction (mag)
C VISEXT = effective visual extinction (mag)
C SHIELD = shielding function
C N_x = column densities for self and mutual shielding (cm2)

      DOUBLE PRECISION HLOSS,STICK
      DOUBLE PRECISION A1,A2,A3,B1,B2,ECHEM
      DOUBLE PRECISION GRAD,GFRAC,ALBEDO,BARRIER,NMONO
      DOUBLE PRECISION SITES,HOPRATIO,DUTY,CRTEMP
      DOUBLE PRECISION TOTSITES,DENSITES,DENSGRAIN
      DOUBLE PRECISION DUSTSUB,ICESUB
      
C HLOSS = rate of H destruction/H2 reformation
C STICK = sticking probability for H atoms
C Ax,Bx,Cx = constants for H-atom sticking coefficients (various units)
C ECHEM = H-atom chemisorption barrier (K)
C GRAD = grain radius (cm)
C GFRAC = fractional abundance of grain (relative to H nuclei density)
C ALBEDO = grain albedo in far UV
C BARRIER = barrier between grain surface sites (A)
C NMONO = number of chemically active monolayers
C SITES = density of grain-surface sites (cm-2)
C HOPRATIO = ratio of surface hopping to desorption energies
C DUTY = duty cycle of the the grain for cosmic-ray-induced thernal desorption
C CRTEMP = maximum temperature reached by grain per cosmic-ray impact (K)
C TOTSITES = total number of surface sites per grain
C DENSITES = number density of surface sites (cm-3)  
C DENSGRAIN = number density of grains (cm-3)  
C DUSTSUB = dust sublimation temperature (K)
C ICESUB = ice sublimation temperature (rough estimate - K)
            
      DOUBLE PRECISION A_RDIFF,B_RDIFF,A_DES,B_DES
      DOUBLE PRECISION A_MASS,B_MASS
      DOUBLE PRECISION A_RQUAN,B_RQUAN,KQUAN,BRANCH
      DOUBLE PRECISION A_FREQ,B_FREQ,FREQ,KAPPA

C x_RDIFF = hopping/diffusion rate (s-1)
C x_DES = binding energies (K)
C x_MASS = masses (g)
C x_QUAN = quantum tunnelling rate (s-1)
C KQUAN = quantum tunnelling reaction rate (cm3 s-1)
C BRANCH = branching ratio for reactive desorption     
C X_FREQ = vibrational frequency
C KAPPA = probability for grain-surface reaction
      
      DOUBLE PRECISION PI,HBAR,KERG,AMU,ECHARGE,EMASS,C
      DOUBLE PRECISION YR2SEC,AU2CM,ANG2CM

C PI = 3.14159 ...
C HBAR = Planck's constant (h/2pi; erg s)
C KERG = Boltzmann's constant (erg K-1)
C AMU = atomic mass unit
C ECHARGE = electron charge
C EMASS = electron mass
C C = speed of light (cm s-1)
C YR2SEC = conversion factor from years to seconds (sec yr-1)
C AU2CM = conversion factor from AU to cm (cm AU-1)
C ANG2CM = conversion factor from Angstroms to cm (cm A-1)

       DOUBLE PRECISION DUMMY          
     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Declaration of common blocks
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                
      COMMON/BLK1/K,X,TOTAL,HLOSS,DENSITES,NMONO
     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Data statements
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      DATA PI,AMU,KERG,HBAR/3.141593,1.6735E-24,1.3806E-16,1.0546E-27/
      DATA YR2SEC,ECHARGE /3.1536E7,4.80320425E-10/
      DATA C,EMASS /2.99792458E+10,5.48579909E-04/
      DATA AU2CM,ANG2CM /1.49597871E13,1.0E-08/
      DATA DUSTSUB,ICESUB /2300.0,150.0/
 
      DATA NCON,TOTAL/1,0.00/
        
      DATA (ELEM(I),I=1,NELEM)
     *   /'He','C','N','O','Si','S','Fe','Na','Mg','Cl'/
     
      DATA (FRAC(I),I=1,NELEM)
     *   /9.75e-02,1.40e-04,7.50e-05,3.20e-04,8.00e-09,8.00e-08,
     *    3.00e-09,2.00e-09,7.00e-09,4.00e-09 /
            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Open and read reaction file
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      OPEN(UNIT=1,FILE=REACFILE)
      
      I = 1
      NREAC = 0
      
100   READ(1,10,END=101) RE1(I),RE2(I),RE3(I),
     1   P1(I),P2(I),P3(I),P4(I),P5(I),
     2   ALPHA(I),BETA(I),GAMMA(I),LOWTEMP(I),UPTEMP(I),
     3   RTYPE(I),SCALE(I),PRFILE(I)
          
      I = I + 1
      NREAC = NREAC + 1
            
      GO TO 100
      
101   CONTINUE

10    FORMAT(5X,8(A10),E8.2,F9.2,F10.1,2(I5),X,I2,X,1PE8.2,X,A12) 

      CLOSE(UNIT=1)  
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Open and read species file
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      OPEN(UNIT=1,FILE=SPECFILE)
      
      I = 1
      NSPEC = 0

200   READ(1,*,END=201) DUMMY,SPEC(I),Y0(I),MASS(I)     

      I = I + 1
      NSPEC = NSPEC + 1
      
      GO TO 200
      
201   CONTINUE

      CLOSE(UNIT=1)
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Open and read binding energies file
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      OPEN(UNIT=1,FILE=BINDFILE)
      
      I = 1
      NICE = 0

300   READ(1,*,END=301) SICE(I),BIND(I),DUMMY,YIELD(I)

      I = I + 1
      NICE = NICE + 1
      
      GO TO 300
      
301   CONTINUE

      CLOSE(UNIT=1)
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Set total number of species (include conserved species)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      NTOT = NSPEC+NCON
           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Set start and end times
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           
C Convert output times to seconds

      TBEGIN = TSTART*YR2SEC
      TEND   = TFINAL*YR2SEC 
           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Open and read radiation parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      OPEN(UNIT=1,FILE=RADFILE)
      
      READ(1,*)
      READ(1,*)
      READ(1,*) GISM
      READ(1,*)
      READ(1,*) HABING
      READ(1,*)
      READ(1,*) ZISM
      READ(1,*)
      READ(1,*) CRPHOT
      READ(1,*)
      READ(1,*) GEXT
      READ(1,*)
      READ(1,*) STARUV
      READ(1,*)
      READ(1,*) LYMAN
      
      CLOSE(UNIT=1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Open and read grain parameters 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      OPEN(UNIT=1,FILE=GRAINFILE)

      READ(1,*)
      READ(1,*)
      READ(1,*) GRAD
      READ(1,*)
      READ(1,*) ALBEDO
      READ(1,*)
      READ(1,*) BARRIER
      READ(1,*)
      READ(1,*) SITES
      READ(1,*)
      READ(1,*) HOPRATIO
      READ(1,*)
      READ(1,*) NMONO
      READ(1,*)
      READ(1,*) DUTY
      READ(1,*)
      READ(1,*) CRTEMP
      READ(1,*)
      READ(1,*) BRANCH
      
      CLOSE(UNIT=1)

C Total number of surface sites per grain (no units)  

      TOTSITES = SITES*(4*PI*GRAD**2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Open and read reaction switches
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      OPEN(UNIT=1,FILE=SWITCHFILE)
      
      READ(1,*)
      READ(1,*)
      READ(1,*) NSWITCH
      READ(1,*) (SWITCH(I),I=1,NSWITCH)
      
      CLOSE(UNIT=1)
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     UV photon flux
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C UV photon flux i.e. in units of photons cm-2 s-1
C Can be calculated between 912 and 2000 Angstroms --> neglects Ly-alpha!
C G_FUV is in units of the interstellar radiation field
C     
C Scale UV photon flux by the Habing field (ISRF)

      UVPHOT = GFUV*HABING

C Calculate UV photon flux explicitly
C Restrict calculation to UV only i.e. lambda < 400 nm or 4000 A

      IF(PH.EQ.1) THEN

      UVPHOT = 0.0
      
      DO I=2,NUV 
         IF (UVANG(I).LE.2100) THEN    
         UVPHOT = UVPHOT + 0.5*(UVFLUX(I-1) + UVFLUX(I))*
     *      (UVANG(I) - UVANG(I-1))      
         END IF
      END DO

C Multiply by 4pi steradians
     
      UVPHOT = 4.0*PI*UVPHOT
      
      END IF

C Add on internally generated UV photons

      UVPHOT = UVPHOT + CRPHOT*(ZETACR/ZISM) 
      
C For self- and mutual shielding of H2, CO, and N2 photoreaction rates
C Calculate effective UV extinction
C Note factor of pi already taken into account

      UVEXT = GFUV*GISM/
     *    (4*(STARUV+LYMAN)/
     *    ((R*AU2CM)**2.0 + (Z*AU2CM)**2.0) 
     *   + GEXT)
                  
      UVEXT = -LOG(UVEXT)
      
C Correct for infinite extinction
      
      IF(GFUV.EQ.0.0) UVEXT = 1E+10
      
C Tau_{UV}/Av ~ 3.02

      VISEXT = UVEXT/3.02
      
C H_{nuc} ~ 1.59e+21 Av cm-2
C N(H2) ~ 0.5 x H_{nuc}
C N(CO) ~ 1e-5 x N(H2)
C N(N2) ~ 1e-5 x N(H2)
      
C Correct for large values
      
      IF(VISEXT.GT.1E+10) VISEXT = 1E+10
            
      N_H2 = 0.5*1.59E+21*VISEXT
      N_CO = 1e-5*N_H2
      N_N2 = 1e-5*N_H2
           
C Convert LYMAN from erg s-1 to photons cm-2 s-1

      LYWAV = 1215*ANG2CM
      
      IF(LYMAN.NE.0.0) LYMAN = (LYMAN/(HBAR*2*PI*C/LYWAV))/
     *   ((R*AU2CM)**2.0 + (Z*AU2CM)**2.0)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      X-ray photon flux
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Xray photon flux i.e. in units of x-rays cm-2 s-1
C Calculate over entire X-ray energy range

      XRPHOT = 0.0
      
      IF(XR.EQ.1) THEN
      
      DO I=1,NXR      
         XRPHOT = XRPHOT + 0.5*(XRFLUX(I) + XRFLUX(I+1))*
     *      (XREV(I+1) - XREV(I))      
      END DO

C Multiply by 4pi steradians
      
      XRPHOT = 4.0*PI*XRPHOT
      
      END IF
                        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Convert mass density to number density
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      X_H2 = 0.0
      X_HE = 0.0
      GFRAC = 0.0
      
      DO I=1,NSPEC      
         IF(SPEC(I).EQ.'H2') X_H2 = Y0(I)
         IF(SPEC(I).EQ.'He') X_HE = Y0(I)
         IF(SPEC(I).EQ.'GRAIN0') GFRAC = GFRAC + Y0(I)      
         IF(SPEC(I).EQ.'GRAIN-') GFRAC = GFRAC + Y0(I)      
      END DO

C RATIO = particle number density / hydrogen nuclei density 
      RATIO = 1 - X_H2 + X_HE

C H and H2 number densities relative to particle number density
      Y_H2 = X_H2/RATIO
      Y_HE = X_HE/RATIO
           
C Mean particle mass      
      M_MU = 2*Y_H2 + 4*Y_HE

C Particle number density
CC      DENS = DENSMASS/(AMU*M_MU)
      DENS = DENSMASS

C H nuclei density for initial abundances
      DENSNUC = DENS/RATIO

C Number density of grain-surface sites (cm-3)  
      DENSITES = TOTSITES*GFRAC*DENSNUC 

C Total number density of grains (cm-3)  
      DENSGRAIN = GFRAC*DENSNUC 
                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCC   BEGIN INITIALISING MODEL   CCCCCCCCCCCCCCCCCCCCCCCC    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Initial abundances
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DO I=1,NSPEC            
         Y(I) = Y0(I)*DENSNUC 
      END DO
                                         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Add conserved species to species list and assign initial abundance
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Add on electrons
      
      DO I=1,NCON
         X(I) = 0.0
     	 SPEC(NSPEC+I)='e-'         
      END DO	

      TOTAL = X(1)*DENSNUC
           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Recalculate X-ray ionisation rate if X-ray spectrum provided
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
      IF (XR.EQ.1) CALL ZETA_X (NXR,XREV,XRFLUX,NELEM,ELEM,FRAC,ZETAXR)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calculating rate coefficients
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           
      DO I=1,NREAC

C Reaction type 1: default two-body rate coefficient  (cm3 s-1)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      IF(RTYPE(I).EQ.1) THEN
      
      IF ((TEMPGAS.GE.LOWTEMP(I)).AND.(TEMPGAS.LT.UPTEMP(I))) THEN              
         K(I) = ALPHA(I)*((TEMPGAS/300.0)**BETA(I))
     *      *EXP(-GAMMA(I)/TEMPGAS)  
      ELSE
         K(I)=0.0
      END IF
      
      END IF

C Reaction type 2: direct cosmic-ray ionisation  (s-1)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      IF(RTYPE(I).EQ.2) THEN  
         
         K(I) = ALPHA(I)*(ZETACR+ZETAXR)/ZISM 
         
         IF (XR.EQ.1) THEN
         
            DO J=1,NELEM
            
            IF(RE1(I).EQ.ELEM(J)) K(I)= ALPHA(I)*(ZETACR/ZISM)
            
            END DO
         
         END IF            
            
      END IF	

C Reaction type 3: cosmic-ray-induced photoreaction (s-1)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF(RTYPE(I).EQ.3) THEN
         K(I)= ALPHA(I)*((ZETAXR+ZETACR)/ZISM)
     *      *((TEMPGAS/300.0)**BETA(I))
     *      *GAMMA(I)/(1.0-ALBEDO)
      END IF	

C Reaction type 4: photoreaction (s-1)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF(RTYPE(I).EQ.4) THEN

         IF((PH.EQ.0).OR.(PRFILE(I).EQ.'DEFAULT')) THEN
            K(I) = ALPHA(I)*GFUV
         ELSE
            CALL PHOTORATES(PRFILE(I),UVANG,NUV,UVFLUX,PHRATE)
            K(I) = PHRATE  
         END IF

C Self-shielding of H2

         IF(RE1(I).EQ.'H2') THEN
      
         CALL SHIELDING(RE1(I),N_H2,N_CO,N_N2,TEMPGAS,SHIELD) 

         IF(GFUV.EQ.0.0) SHIELD = 0.0

         K(I) = K(I)*SHIELD
      
C Self- and mutual-shielding of CO

         ELSE IF (RE1(I).EQ.'CO') THEN
      
         CALL SHIELDING(RE1(I),N_H2,N_CO,N_N2,TEMPGAS,SHIELD)

         IF(GFUV.EQ.0.0) SHIELD = 0.0

         K(I) = K(I)*SHIELD
      
C Self- and mutual-shielding of N2

         ELSE IF (RE1(I).EQ.'N2') THEN
      
         CALL SHIELDING(RE1(I),N_H2,N_CO,N_N2,TEMPGAS,SHIELD)

         IF(GFUV.EQ.0.0) SHIELD = 0.0

         K(I) = K(I)*SHIELD

         END IF
                  
      END IF

C Reaction type 5: direct X-ray ionisation (s-1) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Only call when X-ray spectrum provided

      IF(RTYPE(I).EQ.5) THEN
         IF(XR.EQ.1) THEN
         
            CALL XR_RATES(RE1(I),NXR,XREV,XRFLUX,XRRATE)         
        
            K(I) = XRRATE   
         
         ELSE                    
         
            K(I) = 0.0                    
         
         END IF          
      END IF
     
C Reaction type 6: grain-cation recombination rate (s-1)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Only allow recombination if dust temperature is lower than dust 
C sublimation temperature

C ALPHA = branching ratio
C GAMMA = mass (amu)

C For X+ + G- ---> X0 + G0 
C Enhancement in rate according to Draine & Sutin (1987)

C K = Racc * { 1 + e^2/(akT) } * { 1 + sqrt(2e^2/(akT + 2e^2)) }

      IF(RTYPE(I).EQ.6) THEN

         K(I) = ALPHA(I)*(PI*GRAD**2.0)
     1      *SQRT(8.0*KERG*TEMPGAS/(PI*AMU*GAMMA(I)))
     2      *(1.0 + (ECHARGE**2.0/(GRAD*KERG*TEMPGAS)))
     3      *(1.0 + SQRT(2.0*ECHARGE**2.0
     4      /((GRAD*KERG*TEMPGAS) + 2.0*ECHARGE**2.0)))
     
         IF(TEMPDUST.GT.DUSTSUB) K(I) = 0.0 
      
      END IF
      
C Reaction type 7: neutral grain accretion rate (s-1)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Only allow accretion if dust temperature is lower than ice 
C sublimation temperature - will strip ice mantle

C ALPHA = branching ratio
C GAMMA = mass (amu)

C Independent of grain charge state (using total grain density)
     
      IF(RTYPE(I).EQ.7) THEN
         
         K(I) = ALPHA(I)*(PI*GRAD**2.0)*GFRAC*DENSNUC
     1      *SQRT(8.0*KERG*TEMPGAS/(PI*AMU*GAMMA(I)))
     
         IF(TEMPDUST.GT.ICESUB) K(I) = 0.0 
      
       END IF

C Reaction type 8: thermal desorption rate (s-1)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C ALPHA = binding energy (K)
C GAMMA = mass (g)

C Use zeroth-order rate (cm-3 s-1)

C Look up binding energy

      IF(RTYPE(I).EQ.8) THEN

         ALPHA(I) = 0.0
         
         DO J=1,NICE              
            IF(SICE(J).EQ.RE1(I)) ALPHA(I) = BIND(J)
         END DO
         
         IF(ALPHA(I).EQ.0.0) STOP
         
         K(I)= SQRT((2*SITES*KERG*ALPHA(I))
     1      /((PI**2)*AMU*GAMMA(I)))
     2      *NMONO*DENSITES
     3      *EXP(-ALPHA(I)/TEMPDUST)
                       
      END IF	

C Reaction type 9: cosmic-ray-induced thermal desorption rate (s-1)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C ALPHA = binding energy (K)
C GAMMA = mass (g)
C DUTY = duty cycle of the grains
C CRTEMP = 70 K

C Look up binding energy

      IF(RTYPE(I).EQ.9) THEN

         ALPHA(I) = 0.0
         
         DO J=1,NICE              
            IF(SICE(J).EQ.RE1(I)) ALPHA(I) = BIND(J)
         END DO
         
         IF(ALPHA(I).EQ.0.0) STOP
           
         K(I)= (ZETACR/ZISM)*DUTY
     1      *SQRT((2*SITES*KERG*ALPHA(I))
     2      /((PI**2)*AMU*GAMMA(I)))
     3      *NMONO*DENSITES
     4      *EXP(-ALPHA(I)/CRTEMP)
          
      END IF

C Reaction type 10: photodesorption rate (s-1)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C ALPHA = branching ratio (no units)
C BETA = photodesorption yield (photon-1)

C Look up photodesorption yield

      IF(RTYPE(I).EQ.10) THEN

	 BETA(I)  = 0.0
         
         DO J=1,NICE              
            IF(SICE(J).EQ.RE1(I)) BETA(I) = YIELD(J)
         END DO
         
C         IF(ALPHA(I).EQ.0.0) STOP
          
         K(I) = (UVPHOT+XRPHOT)*ALPHA(I)*BETA(I)*
     *      NMONO*(4.0*PI*GRAD**2.0)*GFRAC*DENSNUC
               
      END IF

C Reaction type 11: grain-surface cosmic-ray-induced photoreaction rate (s-1)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Use same rates as for gas-phase process

      IF(RTYPE(I).EQ.11) THEN      
         K(I)= ALPHA(I)*((ZETAXR+ZETACR)/ZISM)
     *      *((TEMPGAS/300.0)**BETA(I))
     *      *GAMMA(I)/(1.0-ALBEDO)
     
C Direct process

         IF(GAMMA(I).EQ.0.0) K(I)=ALPHA(I)*((ZETAXR+ZETACR)/ZISM)
     
      END IF	

C Reaction type 12: grain-surface photoreaction rate (s-1)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Use same rates as for gas-phase process

      IF(RTYPE(I).EQ.12) THEN 
         IF((PH.EQ.0).OR.(PRFILE(I).EQ.'DEFAULT')) THEN
            K(I) = ALPHA(I)*GFUV
         ELSE
            CALL PHOTORATES(PRFILE(I),UVANG,NUV,UVFLUX,PHRATE)
            K(I) = PHRATE  
         END IF
      END IF   

C Reaction type 13: grain-surface two-body reaction rate (cm3 s-1)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C ALPHA = reaction barrier (K)
C BETA = branching ratio 

C Look up binding energies and masses

      IF((RTYPE(I).EQ.13).OR.(RTYPE(I).EQ.14)) THEN

         A_RDIFF = 0.0
         B_RDIFF = 0.0
             
C Select for first reactant

         DO J=1,NICE           
            IF (RE1(I).EQ.SICE(J)) A_DES = BIND(J)
         END DO
             
         DO J=1,NSPEC           
            IF (RE1(I).EQ.SPEC(J)) A_MASS = MASS(J)
         END DO
         
         IF(A_DES.EQ.0.0) STOP

C Calculate classical diffusion/hopping rate for first reactant
             
         A_FREQ = SQRT((2*SITES*KERG*A_DES)/((PI**2)*AMU*A_MASS))
         
         A_RDIFF = A_FREQ*EXP(-(A_DES*HOPRATIO)/TEMPDUST)

         IF ((RE1(I).EQ.'GH').OR.(RE1(I).EQ.'GH2')) THEN

C Calculate quantum diffusion rate (for H and H2 only)

            A_RQUAN = A_FREQ*EXP(-2*(BARRIER/HBAR)
     1         *SQRT(2*A_MASS*AMU*HOPRATIO*A_DES*KERG))

C Select quantum tunnelling if faster than classical diffusion/hopping

            IF(A_RQUAN.GT.A_RDIFF) A_RDIFF = A_RQUAN
                
         END IF

C Set vibrational frequency for reaction-diffusion competition

         FREQ = A_FREQ
             
C Select for second reactant
                           
         DO J=1,NICE           
            IF (RE2(I).EQ.SICE(J)) B_DES = BIND(J)
         END DO
             
         DO J=1,NSPEC           
            IF (RE2(I).EQ.SPEC(J)) B_MASS = MASS(J)
         END DO

         IF(B_DES.EQ.0.0) STOP

C Calculate classical diffusion/hopping rate for second reactant
                          
         B_FREQ = SQRT((2*SITES*KERG*B_DES)/((PI**2)*AMU*B_MASS))
         
         B_RDIFF = B_FREQ*EXP(-(B_DES*HOPRATIO)/TEMPDUST)
         
         IF ((RE2(I).EQ.'GH').OR.(RE2(I).EQ.'GH2')) THEN

C Calculate quantum diffusion rate (for H and H2 only)

            B_RQUAN = B_FREQ*EXP(-2*(BARRIER/HBAR)
     1         *SQRT(2*B_MASS*AMU*HOPRATIO*B_DES*KERG))

C Select quantum tunnelling if faster than classical diffusion/hopping

            IF(B_RQUAN.GT.B_RDIFF) B_RDIFF = B_RQUAN
          
         END IF

C Choose largest vibrational frequency for reaction-diffusion competition

         IF (B_FREQ.GT.FREQ) FREQ = B_FREQ
                                              
C Calculate classical reaction rate probability

         KAPPA = EXP(-ALPHA(I)/TEMPDUST)
                   
         IF((RE1(I).EQ.'GH').OR.(RE1(I).EQ.'GH2').OR.
     *      (RE2(I).EQ.'GH').OR.(RE2(I).EQ.'GH2')) THEN

C Calculation quantum tunnelling reaction probability (for H and H2 only)
         
         KQUAN = EXP(-2*(BARRIER/HBAR)*SQRT((2*AMU*KERG)
     1      *((A_MASS*B_MASS)/(A_MASS+B_MASS))*ALPHA(I)))
    
C Select quantum tunnelling rate if faster than classical reaction rate
              
         IF(KQUAN.GT.KAPPA) KAPPA = KQUAN
    
         END IF
         
C Calculate probability of reaction (reaction-diffusion competition)
C See e.g., Garrod & Pauly 2011 (equation 6)
C SWITCH OFF - NEEDS TESTED

C      KAPPA=FREQ*KAPPA/(FREQ*KAPPA+A_RDIFF+B_RDIFF)    

C Calculate reaction rate
      
      K(I) = BETA(I)*KAPPA*(A_RDIFF/TOTSITES+B_RDIFF/TOTSITES)*
     *(NMONO**2)*(DENSITES**2)/(GFRAC*DENSNUC)
          
      IF(TEMPDUST.GT.DUSTSUB) K(I) = 0.0 
               
      END IF             
                        
C Reaction type 14: reactive desorption (cm3 s-1)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C BRANCH = branching ratio for reactive desorption

      IF(RTYPE(I).EQ.14) K(I) = K(I)*BRANCH

C  Reaction type 15: three-body association (cm3 s-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C ALPHA = pre-exponential factor (cm6 s-1)
C BETA = temperature exponent
C GAMMA =  reaction barrier (K)
C Multiply by number density for quasi-two-body rate (cm3 s-1)

      IF(RTYPE(I).EQ.15) THEN

      IF ((TEMPGAS.GE.LOWTEMP(I)).AND.(TEMPGAS.LT.UPTEMP(I))) THEN              
          K(I) = ALPHA(I)*((TEMPGAS/300.0)**BETA(I))
     *       *EXP(-GAMMA(I)/TEMPGAS)*DENS
      END IF
         
      END IF

C  Reaction type 16: collisional dissociation (s-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C ALPHA = pre-exponential factor (cm3 s-1)
C BETA = temperature exponent
C GAMMA =  reaction barrier (K)

      IF(RTYPE(I).EQ.16) THEN

          K(I) = ALPHA(I)*((TEMPGAS/300.0)**BETA(I))
     *       *EXP(-GAMMA(I)/TEMPGAS)*DENS

      END IF

C  Reaction type 17: collisional de-excitation of H2*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Slightly different form for the rate coefficient

C ALPHA = pre-exponential factor (cm3 s-1)
C BETA = temperature exponent
C GAMMA =  reaction barrier (K)

      IF(RTYPE(I).EQ.17) THEN
          
          K(I) = ALPHA(I)*((TEMPGAS/300.0)**BETA(I))*
     *       EXP(-GAMMA(I)/(TEMPGAS+1200.0))  

      END IF

C  Reaction type 18: Lyman-alpha photorate
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Slightly different form for rate coefficient

C ALPHA = Lyman-alpha cross section (cm-2)
C SCALE = branching ratio

      IF(RTYPE(I).EQ.18) THEN

         K(I) = ALPHA(I)*SCALE(I)*LYMAN*EXP(-UVEXT)      
      
         IF (LYMAN.EQ.0.0) K(I) = 0.0         

      END IF 

C  Reaction type 19: Radiative de-excitation of H2*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF(RTYPE(I).EQ.19) THEN

         K(I) = ALPHA(I)         

      END IF 

C  Reaction type 20: Grain electron capture rate
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
      IF(RTYPE(I).EQ.20) THEN

         K(I) = (PI*GRAD**2.0)*SQRT(8.0*KERG*TEMPGAS/(PI*AMU*EMASS))
     
         IF(TEMPDUST.GT.DUSTSUB) K(I) = 0.0 
      
      END IF 

C  Reaction type 30: CH3OH co-desorption with CO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Special test case for co-desorption

C ALPHA = yield per CO molecule
C Rate is calculated by mulitplying yield by CO thermal desorption rate 
C   and fraction of CH3OH in top two monolayers of the ice mantle  
     
      IF(RTYPE(I).EQ.30) THEN

         K(I) = ALPHA(I)
           
      END IF 

C Switch off reactions here using switches
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
         DO J=1,NSWITCH         
            IF(RTYPE(I).EQ.SWITCH(J)) K(I) = 0.0         
         END DO
      
      END DO
                                                            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     HLOSS term to model H --> H2 formation on grain surface
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Old term

C      HLOSS = (5.2E-17)*DENSNUC*(TEMPGAS/300)**0.5

C Use temperature-dependent sticking coefficients from Sha et al. 2005 
C and Cuppen et al. 2010
C Need to consider both physisorption and chemisorption

C S_{ph} = (1 + b1*sqrt(Tgas + Tgrain) + b2*Tgas - b3*Tgas**2.0)**-1
C b1 = 4.2 x 10^{-2} K-0.5
C b2 = 2.3 x 10^{-3} K-1
C b3 = 1.3 x 10^{-7} K-2

C S_{ch} = exp{-Ech/Tgas}*(1 + c1*sqrt(Tgas + Tgrain) + c2*Tgas**4.0)**-1
C Ech = 0.15 eV = 1741 K 
C c1 = 5 x 10^{-2} K-0.5
C c2 = 1 x 10^{-14} K-4
      
      A1 = 4.2E-02
      A2 = 2.3E-03
      A3 = 1.3E-07
      
      STICK = (1 + A1*SQRT(TEMPGAS+TEMPDUST) + 
     *   A2*TEMPGAS - A3*(TEMPGAS**2.0))**(-1.0)
      
      ECHEM = 1741.0
      B1 = 5.0E-02
      B2 = 1.0E-14
      
      STICK = STICK + EXP(-ECHEM/TEMPGAS)/
     *   (1 + B1*SQRT(TEMPGAS+TEMPDUST) + B2*TEMPGAS**4.0)
            
      HLOSS = STICK*(PI*GRAD**2.0)*GFRAC*DENSNUC
     *      *SQRT(8.0*KERG*TEMPGAS/(PI*AMU))
     
C Turn off explicit H2 formation if HLOSS is switched on

      IF(HLOSS.NE.0.0) THEN
     
      DO I=1,NREAC      
        IF((RE1(I).EQ.'GH').AND.(RE2(I).EQ.'GH')) THEN 
           K(I) = 0.0
        END IF         

      END DO
     
      END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCC     BEGIN INTEGRATION    CCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
             
      NTIME = 1

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Call DVODE integrator
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Set time counter to zero      
      T0 = 0.0
             
20    CONTINUE

      CALL DRIVE(NSPEC,T0,Y,TBEGIN,DENS)
           
      TIME(NTIME)=TBEGIN/YR2SEC
                 
C Output abundances relative to total H nuclei density
          
      DO I=1,NSPEC
         ABUN(I,NTIME) = Y(I)/DENSNUC
         IF(ABUN(I,NTIME).LT.1D-50) ABUN(I,NTIME) = 0.0
      END DO

      DO I=1,NCON
      	ABUN(NSPEC+I,NTIME) = X(I)/DENSNUC
      END DO 
      
C Write time to standard output

      WRITE(*,*) 'Time = ', time(ntime)
      WRITE(*,*)
     
C If required call ANALYSE subroutine

      IF(NTIME.EQ.NANA) THEN 
      
         WRITE(*,*)
         WRITE(*,*) '           Analyse called!'
         WRITE(*,*)
            
         CALL ANALYSE(NSPEC,SPEC,NREAC,Y,X,K,TBEGIN,
     *      RE1,RE2,RE3,P1,P2,P3,P4) 
         WRITE(*,*)

      END IF   
         
C Set next values for recall to DVODE
      TBEGIN=DLOG10(TBEGIN)+DELTA
     
      TBEGIN=10.0**TBEGIN
            
      NTIME = NTIME + 1  
           
C Repeat until end time reached
      IF(TBEGIN.LE.TEND) GO TO 20
    
      NTIME = NTIME - 1
            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCC     RETURN    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      RETURN
      
      END SUBROUTINE MAIN
         
