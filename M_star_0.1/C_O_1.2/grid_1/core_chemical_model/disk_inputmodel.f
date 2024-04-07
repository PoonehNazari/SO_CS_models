
       PROGRAM Disk_Chemistry

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Declaration of variables
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       
       IMPLICIT NONE
       
       INTEGER I,J,L

C I,J,L = counters       
              
       INTEGER D,E,F,G,H             
       PARAMETER (D=2700)
       PARAMETER (E=50)
       PARAMETER (F=1000)
       PARAMETER (G=1000)
       PARAMETER (H=10000)

C D = max number of grid points
C E = max number of UV wavelength/X-ray energy grids
C F = max number of time steps
C G = max number of species
C H = max number of reactions

       INTEGER NGRID,NUV,NXR

C NGRID = number of disk model grid points
C NUV = number of UV radiation field grid points
C NXR = number of X-ray radiation field grid points
       
       CHARACTER*90 INPFILE
       CHARACTER*90 PHYSFILE,UVFILE,XRFILE
       CHARACTER*90 REACFILE,SPECFILE,BINDFILE
       CHARACTER*90 GRAINFILE,RADFILE,SWITCHFILE
       CHARACTER*90 RATESFILE,ANAFILE,OUTFILE

C INPFILE = master file containing filenames 
C PHYSFILE = disk physical input filename
C UVFILE = disk UV radiation field filename
C XRFILE = disk X-ray radiation field filename
C REACFILE = reaction rates filename
C SPECFILE = species filename
C BINDFILE = binding energy filename
C GRAINFILE = grain parameters filename
C RADFILE = radiation field parameters filename
C SWITCHFILE = switches for chemistry (on/off)
C RATESFILE = reaction rates output filename
C OUTFILE = molecular abundance output filename
C ANAFILE = analyse output filename
       
       DOUBLE PRECISION R(D),Z(D),DENSMASS(D),TEMPGAS(D),TEMPDUST(D)
       DOUBLE PRECISION GFUV(D),ZETACR(D),ZETAXR(D),FXR(D),AV

C R = disk radius (AU)       
C Z = disk height (AU)
C DENSMASS = disk mass density (g cm-3)
C TEMPGAS = disk gas temperature (K)
C TEMPDUST = disk dust temperature (K)
C GFUV = disk wavelength-integrated UV field (relative to ISRF)  
C ZETACR = disk cosmic-ray ionisation rate (s-1)
C ZETAXR = disk x-ray ionisation rates (s-1)
C FXR = energy-integrated x-ray flux (erg cm-2 s-1) 
C AV = visual extinction (mag)

       INTEGER PH,XR
       DOUBLE PRECISION UVANG(E),XREV(E) 
       DOUBLE PRECISION UVGRID(D,E),XRGRID(D,E)   
       DOUBLE PRECISION UVFLUX(E),XRFLUX(E)
       DOUBLE PRECISION GISM   
       
C PH = switch for spectrum-dependent photochemistry
C XR = switch for spectrum-dependent X-ray ionisation rates
C UVANG = UV wavelength grid (A)
C XREV = X-ray energy grid (eV) 
C UVGRID = wavelength-dependent UV radiation field (photons cm-2 s-1 A-1)
C XRGRID = energy-dependent X-ray radiation field (photons cm-2 s-1 eV-1) 
C UVFLUX = 1-D UV radiation field (photons cm-2 s-1 A-1)
C XRFLUX = 1-D X-ray radiation field (photons cm-2 s-1 eV-1) 
C GISM = interstellar radiation field (erg cm-2 s-1)

       INTEGER NTIME,NTOT,NREAC,NANA
       DOUBLE PRECISION TSTART,TFINAL,DELTA
       DOUBLE PRECISION TIME(F),ABUN(G,F),DENS(D)
       DOUBLE PRECISION K(H),X(1),TOTAL,HLOSS,DENSITES,NMONO
       CHARACTER*10 SPEC(G)

C NTIME = number of time steps
C NTOT = total number of species
C NREAC = number of reactions
C NANA = time step for analyse subroutine
C TSTART = first output time (years)
C TFINAL = final output time (years)
C TIME = output time steps (years)
C ABUN = output abundances (relative to total Hnuc density)
C DENS = output gas number density (cm-3)
C K = reaction rate coefficients (various units)
C NMONO = number of chemically active monolayers
C SPEC = species names               
       
       DOUBLE PRECISION DUMMY          

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Declaration of common blocks
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                
      COMMON/BLK1/K,X,TOTAL,HLOSS,DENSITES,NMONO
        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       
                           
       WRITE(*,*)
       WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++'
       WRITE(*,*) '++++++  PROTOPLANETARY DISK CHEMICAL MODEL  +++++++'
       WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++'
       WRITE(*,*)
       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Initialise model parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       NGRID = 0              
       NUV = 0
       NXR = 0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Read master input file from standard input
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              
      CALL GETARG(1,INPFILE)
      
      INPFILE = TRIM(INPFILE)
      
      WRITE(*,'(X,A22,X,A90)') 'Input file           =', INPFILE
            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Open and read file containing filenames
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            
      OPEN(UNIT=1,FILE=INPFILE)
      
      READ(1,*)
      READ(1,*)
      READ(1,'(A90)') PHYSFILE
      READ(1,*)
      READ(1,'(A90)') UVFILE
      READ(1,*)
      READ(1,'(A90)') XRFILE
      READ(1,*)
      READ(1,'(A90)') REACFILE
      READ(1,*)
      READ(1,'(A90)') SPECFILE
      READ(1,*)
      READ(1,'(A90)') BINDFILE
      READ(1,*)
      READ(1,'(A90)') GRAINFILE
      READ(1,*)
      READ(1,'(A90)') RADFILE
      READ(1,*)
      READ(1,'(A90)') SWITCHFILE
      READ(1,*)
      READ(1,'(A90)') OUTFILE
      READ(1,*)
      READ(1,'(A90)') RATESFILE
      READ(1,*)
      READ(1,'(A90)') ANAFILE
      READ(1,*)
      READ(1,*) TSTART
      READ(1,*)
      READ(1,*) TFINAL
      READ(1,*)
      READ(1,*) DELTA
      READ(1,*)
      READ(1,*) NANA
            
      CLOSE(UNIT=1)
      
      PHYSFILE   = TRIM(PHYSFILE)
      UVFILE     = TRIM(UVFILE)
      XRFILE     = TRIM(XRFILE)
      REACFILE   = TRIM(REACFILE)
      SPECFILE   = TRIM(SPECFILE)
      BINDFILE   = TRIM(BINDFILE)
      GRAINFILE  = TRIM(GRAINFILE)
      RADFILE    = TRIM(RADFILE)      
      SWITCHFILE = TRIM(SWITCHFILE)      

      OUTFILE    = TRIM(OUTFILE)      
      RATESFILE  = TRIM(RATESFILE)      
      ANAFILE     = TRIM(ANAFILE)      

      WRITE(*,*)
      WRITE(*,'(X,A22,X,A90)') 'Physical conditions  =', PHYSFILE
      WRITE(*,'(X,A22,X,A90)') 'UV field             =', UVFILE
      WRITE(*,'(X,A22,X,A90)') 'X-ray field          =', XRFILE
      WRITE(*,'(X,A22,X,A90)') 'Chemical network     =', REACFILE
      WRITE(*,'(X,A22,X,A90)') 'Species file         =', SPECFILE
      WRITE(*,'(X,A22,X,A90)') 'Binding energies     =', BINDFILE
      WRITE(*,'(X,A22,X,A90)') 'Grain parameters     =', GRAINFILE
      WRITE(*,'(X,A22,X,A90)') 'Radiation parameters =', RADFILE
      WRITE(*,'(X,A22,X,A90)') 'Reaction switches    =', SWITCHFILE
      WRITE(*,*)      
      WRITE(*,'(X,A22,X,A90)') 'Output file          = ',OUTFILE
      WRITE(*,'(X,A22,X,A90)') 'Reaction rates file  = ',RATESFILE
      WRITE(*,'(X,A22,X,A90)') 'Analyse output file  = ',ANAFILE
      WRITE(*,*)   
      WRITE(*,'(X,A22,X,1PE9.2,A3)') 'First output time    = ',TSTART, 
     *   ' YR'
      WRITE(*,'(X,A22,X,1PE9.2,A3)') 'Final output time    = ',TFINAL, 
     *   ' YR'
      WRITE(*,'(X,A22,X,I3)') 'Timestep for analyse = ',NANA
      WRITE(*,*)   
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Open and read physical conditions file
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      WRITE(*,*) 'Opening and reading physical conditions file ...'
      WRITE(*,*)
      
      OPEN(UNIT=1,FILE=PHYSFILE)
      
      I = 1
      NGRID = 0
      
100   READ(1,*,END=101) R(I),Z(I),DENSMASS(I),TEMPGAS(I),TEMPDUST(I),
     *   GFUV(I),ZETAXR(I),ZETACR(I),FXR(I)
           
      I = I + 1
      NGRID = NGRID + 1
      
      GO TO 100

101   CONTINUE

      CLOSE(UNIT=1)
                     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Read in UV and X-ray spectral data
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Set switch for wavelength-dependent photochemistry

      IF(UVFILE.NE.'DEFAULT') THEN
      
      WRITE(*,*) 'Opening and reading UV field grid ...'
      WRITE(*,*)
      
      PH = 1

C Open and read UV wavelength grid

      OPEN(UNIT=1, FILE='uv_wavelength.dat')

      I = 1
       
200   READ(1,*,END=201) UVANG(I)

      NUV = NUV + 1
      I = I + 1
        
      GO TO 200

201   CONTINUE
                
      CLOSE(UNIT=1)   
       
C Open and read wavelength-dependent UV field

      I = 1

      OPEN(UNIT=1,FILE=UVFILE)  

202   READ(1,*,END=203) (UVGRID(I,J), J=1,NUV)    

      I = I + 1
      
      GO TO 202

203   CONTINUE

      CLOSE(UNIT=1)  

      ELSE
      
      WRITE(*,*) 'Using integrated UV field ...'
      WRITE(*,*)
      
      PH = 0
      
      END IF
                                  
C Set switch for energy-dependent X-ray ionisation

      IF(XRFILE.NE.'DEFAULT') THEN
      
      WRITE(*,*) 'Opening and reading X-ray field grid ...'
      WRITE(*,*)
      
      XR = 1

C Open and read X-ray energy grid

      OPEN(UNIT=1, FILE='xray_energy.dat')
       
      I = 1
       
204   READ(1,*,END=205) XREV(I)

      NXR = NXR + 1
      I = I + 1
       
      GO TO 204                      

205   CONTINUE  

      CLOSE(UNIT=1)     
       
C Open and read energy-dependent X-ray field
       
      I = 1

      OPEN(UNIT=1,FILE=XRFILE)  
        
206   READ(1,*,END=207) (XRGRID(I,J), J=1,NXR)    

      I = I + 1
      
      GO TO 206

207   CONTINUE

      CLOSE(UNIT=1)  

      ELSE
      
      WRITE(*,*) 'Using integrated X-ray field ...'
      WRITE(*,*)
      
      XR = 0
      
      END IF
                          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Solve chemistry for each grid point
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       
      WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(*,*) 'Beginning to run model ...'
      WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(*,*)
      WRITE(*,*) 'Preparing output files ...'
      WRITE(*,*)
      
      OPEN(UNIT=9, FILE=OUTFILE)
      OPEN(UNIT=10,FILE=RATESFILE)
      OPEN(UNIT=11,FILE=ANAFILE)
      
      DO J=9,11
      
      WRITE(J,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(J,*) '++++++  PROTOPLANETARY DISK CHEMICAL MODEL  +++++++'
      WRITE(J,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(J,*)
      
      END DO

C Write progress bar to standard output
      
      WRITE(*,*) 'Progress = ...'

      DO I = 1, NGRID
        
      WRITE(*,'(X,A4,I4,A2,I4)') '... ', I, ' /', NGRID

      DO J=1,NUV
        UVFLUX(J) = UVGRID(I,J)
      END DO
      
      DO J=1,NXR
        XRFLUX(J) = XRGRID(I,J)
      END DO

C Write preamble for analyse routine

         WRITE(11,'(X,A9,1PE10.3,X,A2)')  'RADIUS = ',
     *      R(I), 'AU'
         WRITE(11,'(X,A9,1PE10.3,X,A2)')  'HEIGHT = ',
     *      Z(I), 'AU'
         WRITE(11,*)

C Call chemistry subroutine
       
       CALL MAIN(DENSMASS(I),TEMPGAS(I),TEMPDUST(I),
     1   GFUV(I),ZETACR(I),ZETAXR(I),FXR(I),DENS(I),R(I),Z(I),AV,
     2   UVANG,UVFLUX,NUV,XREV,XRFLUX,NXR,PH,XR,
     3   REACFILE,SPECFILE,BINDFILE,GRAINFILE,RADFILE,SWITCHFILE,
     4   TSTART,TFINAL,NANA,DELTA,NREAC,
     5   GISM,TIME,ABUN,NTIME,SPEC,NTOT)   

C Write preamble for grid points in output files

      DO J=9,10
      
         WRITE(J,'(X,A9,1PE10.3,X,A2)')  'RADIUS = ',
     *      R(I), 'AU'
         WRITE(J,'(X,A9,1PE10.3,X,A2)')  'HEIGHT = ',
     *      Z(I), 'AU'
         WRITE(J,'(X,A10,1PE10.3,X,A4)') 'DENSITY = ',
     *      DENS(I), 'CM-3' 
         WRITE(J,'(X,A18,1PE10.3,X,A2)') 'GAS TEMPERATURE = ',
     *      TEMPGAS(I), 'K' 
         WRITE(J,'(X,A19,1PE10.3,X,A2)') 'DUST TEMPERATURE = ',
     *      TEMPDUST(I), 'K' 
         WRITE(J,'(X,A10,1PE10.3,X,A12)') 'UV FLUX = ',
     *      GFUV(I)*GISM, 'ERG CM-2 S-1'
         WRITE(J,'(X,A11,1PE10.3,X,A12)') 'XRAY FLUX = ',
     *      FXR(I),'ERG CM-2 S-1'
         WRITE(J,'(X,A20,1PE10.3,X,A3)')  'VISUAL EXTINCTION = ',
     *      AV,'MAG'
         WRITE(J,*)

      END DO
       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Output routines
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              
	 WRITE(9,90) (TIME(J),J=1,NTIME)
	 WRITE(9,91) 
          
	 DO L=1,NTOT
            WRITE(9,92) SPEC(L),(ABUN(L,J), J=1,NTIME)  
	 END DO
         WRITE(9,*)
         WRITE(9,*)

90    FORMAT (/,1X,'TIME',6X,150(1X,1PE10.3))
91    FORMAT (1X,150('*'))
92    FORMAT ((1X,A10),150(1X,1PE10.3))
              
              
         DO J=1,NREAC
            WRITE(10,'(I5,X,1PE15.5)') J,K(J)
         END DO   
         
      WRITE(10,*) 
      WRITE(11,*) 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      END DO
                                   
      CLOSE(UNIT=9)
      CLOSE(UNIT=10)
      CLOSE(UNIT=11)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     END OF CODE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      WRITE(*,*)
      WRITE(*,*) '   ... DONE!'
      WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(*,*)
        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     END OF PROGRAM
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       END PROGRAM Disk_Chemistry
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
