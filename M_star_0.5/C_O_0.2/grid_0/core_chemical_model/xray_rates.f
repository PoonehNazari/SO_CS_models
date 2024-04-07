
      SUBROUTINE XR_RATES(SPEC,NE_KEV,E_KEV,XR_FLUX,RATE)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*8 SPEC
      
      PARAMETER (N=100)
      
      DIMENSION E_KEV(N),XR_FLUX(N)
            
C Initialise RATE      
      RATE = 0.0

C Calculate rate as a function of energy      
      DO I=1,NE_KEV-1
           
         CALL CALC_SIGMA(SPEC,E_KEV(I),SIGMA)
         SIGMA1=SIGMA
           
         CALL CALC_SIGMA(SPEC,E_KEV(I+1),SIGMA)
         SIGMA2=SIGMA
           
         RATE = RATE + 0.5*(SIGMA1*XR_FLUX(I)+SIGMA2*XR_FLUX(I+1))
     *      *(E_KEV(I+1)-E_KEV(I))
         
      END DO
      
      RETURN

      END SUBROUTINE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SUBROUTINE CALC_SIGMA(I,E,SIGMA)
C Calculate sigma for each species at given energy E
C Follow method of Verner et al., 1993
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CALC_SIGMA(SPEC,E,SIGMA)
      
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*8 SPEC


      IF (SPEC .EQ. 'H') THEN     
         SIGMA=0.0
         GO TO 200
      
      ELSE IF (SPEC .EQ. 'He') THEN
         L = 0
         P = 6.218
         YA = 9.648
         YW = 0.0
         SIGMA_0 = 2.578E-15
         E_0 = 2.024*1.00E-03  ! keV
         E_KSHELL = 2.342E-02
      
      ELSE IF (SPEC .EQ. 'C') THEN
         L = 0
         P = 1.503
         YA = 5.498E+01
         YW = 0.0
         SIGMA_0 = 7.421E-17
         E_0 = 8.655E+01*1.00E-03
         E_KSHELL = 0.2910
      
      ELSE IF (SPEC .EQ. 'N') THEN
         L = 0
         P = 1.252
         YA = 1.380E+02
         YW = 0.0
         SIGMA_0 = 4.748E-17
         E_0 = 0.1270
         E_KSHELL = 0.4048
      
      ELSE IF (SPEC .EQ. 'O') THEN
         L = 0
         P = 1.083
         YA = 3.812E+02
         YW = 0.0
         SIGMA_0 = 3.237E-17
         E_0 = 0.1774
         E_KSHELL = 0.5373
      
      ELSE IF (SPEC .EQ. 'Si') THEN
         L = 0
         P = 1.102
         YA = 2.580E+02
         YW = 0.0
         SIGMA_0 = 1.184E-17
         E_0 = 5.322E+02*1.00E-03  ! keV
         E_KSHELL = 1.828
      
      ELSE IF (SPEC .EQ. 'S') THEN
         L = 0
         P = 0.8646
         YA = 3.734E+03
         YW = 0.0
         SIGMA_0 = 6.649E-18
         E_0 = 0.8114
         E_KSHELL = 2.456
      
      ELSE IF (SPEC .EQ. 'Fe') THEN
         L = 0
         P = 2.118
         YA = 3.633E+01
         YW = 0.0
         SIGMA_0 = 2.055E-17
         E_0 = 8.044E+02*1.00E-03  ! keV
         E_KSHELL = 7.083
      
      ELSE IF (SPEC .EQ. 'Na') THEN
         L = 0
         P = 0.7736
         YA = 5.642E+07
         YW = 0.0
         SIGMA_0 = 1.119E-17
         E_0 = 0.4216
         E_KSHELL = 1.064
      
      ELSE IF (SPEC .EQ. 'Mg') THEN
         L = 0
         P = 1.952
         YA = 2.374E+01
         YW = 0.0
         SIGMA_0 = 3.561E-17
         E_0 = 2.711E+02*1.00E-03
         E_KSHELL = 1.294
      
      ELSE IF (SPEC .EQ. 'Cl') THEN
         L = 0
         P = 0.7888
         YA = 1.856E+06
         YW = 0.0
         SIGMA_0 = 5.255E-18
         E_0 = 0.9700
         E_KSHELL = 2.805
      
      ELSE IF (SPEC .EQ. 'P') THEN 
         L = 0
         P = 1.063
         YA = 2.562E+02
         YW = 0.0
         SIGMA_0 = 9.167E-18
         E_0 = 0.6472
         E_KSHELL = 2.130
      
      END IF

      IF (E.LT.E_KSHELL) THEN
      
         SIGMA=0.0
      
      ELSE

C Equations (1) and (2) in Verner et al., 1993
      
        Q=5.5+1.0*L-0.5*P
        Y=E/E_0
        
        F=((Y-1.0)**2+YW*YW)*(Y**(-1.0*Q))*(1+SQRT(Y/YA))**(-1.0*P)
        
        SIGMA=SIGMA_0*F
      
      END IF

 200  CONTINUE
      RETURN
      END
