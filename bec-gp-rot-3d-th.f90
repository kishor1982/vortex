!# File name : bec-gp-rot-3d-th.f90
!# Last modified : 20 JANUARY 2019
!# Fortran program for Gross-Pitaevskii equation in three-dimensional ROTATING
!# anisotropic trap by imaginary time propagation (Fortran 90/95 Version)
!#
!# BEC-GP-ROT-OMP programs are developed by:
!#
!# R. Kishor Kumar
!# (Instituto de Fisica, Universidade de Sao Paulo, Sao Paulo, Brazil)
!#
!# Vladimir Loncar, Antun Balaz
!# (Scientific Computing Laboratory, Center for the Study of Complex Systems, Institute of Physics Belgrade, Serbia)
!#
!# Paulsamy Muruganandam
!# (Department of Physics, Bharathidasan University, Tiruchirappalli, Tamil Nadu, India)
!#
!# Sadhan K. Adhikari
!# (Instituto de Fisica Teorica, UNESP - Sao Paulo State University, Brazil)
!#
!# Public use and modification of these codes are allowed provided that the following papers are cited:
!# [1] L. E. Young-S. et al., Comput. Phys. Commun. 220 (2017) 503.
!# [2] P. Muruganandam and S. K. Adhikari, Comput. Phys. Commun. 180 (2009) 1888.
!# [3] D. Vudragovic et al., Comput. Phys. Commun. 183 (2012) 2021.
!# [4] R. Kishor Kumar et al., Comput. Phys. Commun. 195 (2015) 117.
!# [5] B. Sataric et al., Comput. Phys. Commun. 200 (2016) 411.
!# [6] V. Loncar et al., Comput. Phys. Commun. 200 (2016) 406.
!# [7] L. E. Young-S. et al., Comput. Phys. Commun. 204 (2016) 209.
!# [8] V. Loncar et al., Comput. Phys. Commun. 209 (2016) 190.
!#
!# The authors would be grateful for all information and/or comments
!# regarding the use of the programs.
!
!# To compile :
!# (1) Intel Fortran Compiler
! ifort -O3 -openmp -w -mcmodel medium -shared-intel -V
!
!# (2) GNU Fortran (gfortran)
! gfortran -O3 -fopenmp -w
!
!# (3) PGI Fortran (pgfortran)
! pgfortran -O3 -fast -mp=allcores

MODULE COMM_DATA
! NX, NY, NZ : Number of space mesh points (X, Y and Z)
  INTEGER, PARAMETER :: NX = 256, NXX = NX-1, NX2 = NX/2 ! No. of X steps
  INTEGER, PARAMETER :: NY = 256, NYY = NY-1, NY2 = NY/2 ! No. of Y steps
  INTEGER, PARAMETER :: NZ = 32, NZZ = NZ-1, NZ2 = NZ/2 ! No. of Z steps
! NSTP : Number of iterations to introduce the nonlinearity.
! NSTP=0 reads the wave function, /=0 calculates the wave function.
! NPAS : Number of subsequent iterations with fixed nonlinearity.
! NRUN : Number of final time iterations with fixed nonlinearity.
! For real time set NSTP=0, and read imag time wave function from im3d-fin-wf.txt
  INTEGER, PARAMETER :: NSTP = 1, NPAS = 3000000, NRUN = 1, ITER = NPAS/10
  INTEGER, PARAMETER :: NUMBER_OF_THREADS = 0! sets the number of CPU cores to be used
! NUMBER_OF_THREADS = 0 deactivates the command and uses all available CPU cores
  REAL (8), PARAMETER :: PI = 3.14159265358979D0
  INTEGER :: NO_OF_THREADS
END MODULE COMM_DATA

MODULE GPE_DATA
  USE COMM_DATA, ONLY : PI
  REAL (8), PARAMETER :: AHO = 1.0D-6 ! Unit of length (l= 1 MICRON)
  REAL (8), PARAMETER :: Bohr_a0 =  5.2917720859D-11/AHO ! Bohr radius (scaled with AHO)
  COMPLEX (8), PARAMETER :: CI = (0.0D0,1.0D0)

  REAL (8), PARAMETER :: DX = 0.05D0, DY = 0.05D0,  DZ = 0.025D0 ! DX, DY, DZ : Space step and DT : Time step

! INTEGER, PARAMETER :: RANDOM=0 ! No random phase in the initial function
  INTEGER, PARAMETER :: RANDOM=1 ! Random phase included in the initial function


  ! OPTION_RE_IM switches between real and imaginary time propagation
  ! Select parameter  OPTION_RE_IM and time step DT for imaginary and real-time propagation
  INTEGER, PARAMETER :: OPTION_RE_IM = 1    !  Imaginary time propagation
  REAL (8), PARAMETER :: DT = 0.00025d0    !  Imaginary time propagation, Time step
 ! INTEGER, PARAMETER ::  OPTION_RE_IM = 2   !  Real time propagation
 ! REAL (8), PARAMETER :: DT = 0.0001d0    ! Real time propagation, Time step

!  SELECT INITIAL FUNCTION 
!  INTEGER, PARAMETER :: FUNCTION=0 !  GAUSSIAN function with no vortex
   INTEGER, PARAMETER :: FUNCTION=1 !  1VORTEX function with 1 vortex

  INTEGER, PARAMETER  :: NATOMS = 10000 ! Number of Atoms
  REAL (8), PARAMETER :: AS = 3.769458264D0*Bohr_a0 ! Scattering length (in units of Bohr_a0)
  REAL (8), PARAMETER :: GAMMA = 1.D0, NU = 1.D0 ! GAMMA and NU : Parameteres of Trap
  REAL (8), PARAMETER :: LAMBDA = 100.D0 ! Trap aspect ratio
! G_3D : Nonlinearity in the 3D GP equation
  REAL (8), PARAMETER :: G_3D =  4.D0 * PI * AS * NATOMS
  REAL (8), PARAMETER :: G_2D = G_3D * SQRT(LAMBDA/(2*pi))

  REAL (8), PARAMETER :: OMEGA = 0.8D0 ! Angular Frequency

! OPTION  decides which equation to be solved.
! OPTION=1 Solves -psi_xx-psi_yy+V(x,y)psi+G_3D|psi|^2 psi =i psi_t
! OPTION=2 Solves [-psi_xx-psi_yy+V(x,y)psi]/2+G_3D|psi|^2 psi =i psi_t
  INTEGER, PARAMETER :: OPTION = 2

! X(0:NX), Y(0:NY), Z(0:NZ): Space mesh, V(0:NX,0:NY) : Potential, CP(0:NX,0:NY): Wave function
  REAL (8), DIMENSION(:), ALLOCATABLE :: X, X2, Y, Y2, Z, Z2
  REAL (8), DIMENSION(:,:,:), ALLOCATABLE :: V, R2
  COMPLEX (8), DIMENSION(:,:,:), ALLOCATABLE :: CP
  REAL (8) :: G, XOP, GSTP
END MODULE GPE_DATA

MODULE CN_DATA
  USE COMM_DATA, ONLY : NX, NY, NZ
  COMPLEX (8), DIMENSION(0:NX) :: CBP, CBM
  COMPLEX (8), DIMENSION(0:NY) :: CAP, CAM
  COMPLEX (8), DIMENSION(:,:), ALLOCATABLE :: CALA, CALB, CGAA, CGAB
  COMPLEX (8), DIMENSION(:), ALLOCATABLE :: CALC, CGAC
  COMPLEX (8) :: CT0X, CT0Y, CT0, CT0Z, CTMPX, CTMPY
  COMPLEX (8) :: CA0, CB0, CA0R, CB0R, CC0, CC0R, C0
END MODULE CN_DATA

PROGRAM GROSS_PITAEVSKII_SSCN_3D
  USE OMP_LIB
  USE COMM_DATA, ONLY : NX, NY, NZ, NX2, NY2, NZ2, NSTP, NPAS, NRUN, NUMBER_OF_THREADS,ITER
  USE GPE_DATA

  IMPLICIT NONE
! Subroutine INTIIALIZE() used to initialize the space mesh, and the initial
! wave function. Subroutine CALCULATE_TRAP() initializezpotential V(I),
! while subroutine COEF() is used to generate the coefficients for the
! Crank-Nicholson Scheme. The routine NU() performs time progation for the
! non-derivative part and LU() performs time propagation of derivative part.
! NORM() calculates the norm and normalizes the wave function, CHEM() and
! RAD() are used to calculate the chemical potential, energy and the RMS
! radius, respectively. The functions DIFF(), used to calculate the space
! derivatives of the wave function used in CHEM(). The function SIMP() does
! the integration by Simpson's rule.

!------------------------ INTERFACE BLOCKS -----------------------
  INTERFACE
    SUBROUTINE ALLOCATE_VARIABLES()
    END SUBROUTINE ALLOCATE_VARIABLES

    SUBROUTINE FREE_VARIABLES()
    END SUBROUTINE FREE_VARIABLES

    SUBROUTINE INITIALIZE()
      IMPLICIT NONE
    END SUBROUTINE INITIALIZE

    SUBROUTINE CALCULATE_TRAP()
    END SUBROUTINE CALCULATE_TRAP

    SUBROUTINE COEF()
      IMPLICIT NONE
    END SUBROUTINE COEF

    SUBROUTINE CALCNU(CP,DT)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:NX, 0:NY, 0:NZ), INTENT(INOUT) :: CP
      REAL (8), INTENT(IN) :: DT
    END SUBROUTINE CALCNU

    SUBROUTINE LUX(CP)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(INOUT) :: CP
    END SUBROUTINE LUX

    SUBROUTINE LUY(CP)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(INOUT) :: CP
    END SUBROUTINE LUY

    SUBROUTINE LUZ(CP)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(INOUT) :: CP
    END SUBROUTINE LUZ

    SUBROUTINE CHEM(CP, MU, EN)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: CP
      REAL (8), INTENT(OUT) :: MU, EN
    END SUBROUTINE CHEM

    SUBROUTINE WRITE_DENSITY(FUNIT, U2)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FUNIT
      REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: U2
    END SUBROUTINE WRITE_DENSITY

    SUBROUTINE WRITE_3D(FUNIT, U2)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FUNIT
      REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: U2
    END SUBROUTINE WRITE_3D

    SUBROUTINE WRITE_WAVE_FUNCTION(FUNIT, CP)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FUNIT
      COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: CP
    END SUBROUTINE WRITE_WAVE_FUNCTION

    SUBROUTINE NORM(CP, ZNORM)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: ZNORM
    END SUBROUTINE NORM

    SUBROUTINE NORMCHEM(CPR, ZNORM)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(INOUT) :: CPR
      REAL (8), INTENT(OUT) :: ZNORM
    END SUBROUTINE NORMCHEM

    SUBROUTINE RAD(CP2, RMS)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: CP2
      REAL (8), DIMENSION(:), INTENT(OUT) :: RMS
    END SUBROUTINE RAD

    SUBROUTINE DEN1DX(FUNIT, U2, DE1DX)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FUNIT
      REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: U2
      REAL (8), DIMENSION(0:NX), INTENT(OUT) :: DE1DX
    END SUBROUTINE DEN1DX

    SUBROUTINE DEN1DY(FUNIT, U2, DE1DY)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FUNIT
      REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: U2
      REAL (8), DIMENSION(0:NY), INTENT(OUT) :: DE1DY
    END SUBROUTINE DEN1DY

    SUBROUTINE DEN1DZ(FUNIT, U2, DE1DZ)
      USE COMM_DATA, ONLY : NX, NY, NZ
      INTEGER, INTENT(IN) :: FUNIT
      REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: U2 
      REAL (8), DIMENSION(0:NZ), INTENT(OUT) :: DE1DZ
    END SUBROUTINE DEN1DZ
  END INTERFACE
!------------------------ END INTERFACE BLOCKS -------------------

  INTEGER :: I, J, K, H, NO_OF_THREADS, STEP
  REAL (8) :: ZNORM, T1, T2
  REAL (8) :: MU, EN
  REAL (8), DIMENSION(:,:,:), ALLOCATABLE :: CP2
  REAL (8), DIMENSION(4) :: RMS
  REAL(KIND(0.0D0)) :: START_D, END_D
  REAL (8), DIMENSION(0:NX) :: DEN1X
  REAL (8), DIMENSION(0:NY) :: DEN1Y
  REAL (8), DIMENSION(0:NZ) :: DEN1Z
  CHARACTER(256) :: ITER_FILENAME

  INTEGER :: CLCK_COUNTS_BEG, CLCK_COUNTS_END, CLCK_RATE,COUNT_MAX

  CALL SYSTEM_CLOCK ( CLCK_COUNTS_BEG, CLCK_RATE )
  CALL CPU_TIME (T1)

  IF (NUMBER_OF_THREADS /= 0) THEN
    CALL OMP_SET_NUM_THREADS(NUMBER_OF_THREADS)
  END IF

  !$OMP PARALLEL
  !$OMP MASTER
  NO_OF_THREADS = OMP_GET_NUM_THREADS()
  !$OMP END MASTER
  !$OMP END PARALLEL

  SELECT CASE (OPTION)
    CASE (1)
      XOP = 1.0D0
    CASE (2)
      XOP = 2.0D0
    CASE (:0,3:)
      PRINT *, 'ERROR: Wrong option', OPTION
      STOP
  END SELECT

  CALL ALLOCATE_VARIABLES()
  ALLOCATE(CP2(0:NX, 0:NY, 0:NZ))

  CALL INITIALIZE()

  IF (OPTION_RE_IM.EQ.2) THEN
    OPEN(7, FILE = 're3d-out.txt')
  ELSE
    OPEN(7, FILE = 'im3d-out.txt')
  END IF

  IF (OPTION_RE_IM.EQ.2) THEN
    OPEN(4, FILE = 're3d-rms.txt')
  ELSE
    OPEN(4, FILE = 'im3d-rms.txt')
  END IF

  IF (OPTION_RE_IM.EQ.1) THEN
    WRITE(7,900) OPTION, NO_OF_THREADS
    WRITE(4,900) OPTION, NO_OF_THREADS
  ELSE
    WRITE(7,800) OPTION, NO_OF_THREADS
    WRITE(4,800) OPTION, NO_OF_THREADS
  END IF

  WRITE(7,*)
  WRITE(4,*)
  WRITE(7,901) NATOMS, AHO
  WRITE(7,902) AS/Bohr_a0
  WRITE(7,903) G_3D , G_2D
  WRITE(7,904) GAMMA, NU, LAMBDA
  WRITE(7,911) OMEGA
  WRITE(7,*)
  WRITE(7,905) NX + 1, NY + 1, NZ + 1
  WRITE(7,906) DX, DY, DZ
  WRITE(7,907) NSTP, NPAS, NRUN
  WRITE(7,908) DT
  WRITE(7,*)

  900 FORMAT('Imaginary time propagation 3d,   OPTION = ',I0,', NUM_THREADS = ',I0)
  800 FORMAT('Real time propagation 3d,   OPTION = ',I0,', NUM_THREADS = ',I0)

  901 FORMAT('Number of Atoms N = ',I0,', Unit of length AHO = ',F10.8,' m')
  902 FORMAT('Scattering length a = ',F0.6,'*a0')
  903 FORMAT('Nonlinearity G_3D = ',F0.6, ', G_2D = ',F0.6)
  904 FORMAT('Parameters of trap: GAMMA = ',F0.2, ', NU = ',F0.2, ', LAMBDA = ',F0.2)
  911 FORMAT('Parameters of rotation: ANG VEL = ',F6.4)
  905 FORMAT('Space step: NX = ', I0, ', NY = ', I0, ', NZ = ', I0)
  906 FORMAT('            DX = ', F8.6, ', DY = ', F8.6, ', DZ = ', F8.6)
  907 FORMAT('Time step:  NSTP = ',I0,', NPAS = ',I0,', NRUN = ',I0)
  908 FORMAT('            DT = ', F8.6)

  CALL CALCULATE_TRAP() ! CALCULATE_TRAP() initializes the harmonic oscillator potential V.
  CALL COEF() ! COEF() defines the coefficients of Crank-Nicholson Scheme.

  CP2 = CP * CONJG(CP)
  CALL RAD(CP2, RMS) ! RAD() calculates the r.m.s radius RMS
  WRITE (7, 1001)
  WRITE (7, 1002)
  WRITE (7, 1001)

  WRITE (4, 1001)
  WRITE (4, 1012)
  WRITE (4, 1001)
  WRITE (4, 1013) RMS(2:4)
  1001 FORMAT (12X,'-----------------------------------------------------')
  1002 FORMAT (14X,'Iter',5x, 'Norm', 7X, 'Chem', 7X, 'Ener/N', 6X, '<rho>', 3X)
  1003 FORMAT ('Initial: ', 6X, F13.4, 2F12.5, 2F11.5, F11.4)
  1012 FORMAT ('RMS size:', 4x, 'Iter', 8x, '<x>', 13X, '<y>', 13X, '<z>')
  1013 FORMAT ('Initial:',9X, F14.5, 2F16.5)
  FLUSH(7)
  FLUSH(4)

! IF (OPTION_RE_IM.EQ.2) THEN
!   OPEN(21, FILE = 're3d-initial-den.txt')
! ELSE
!   OPEN  (21, FILE = 'im3d-initial-den.txt')
! END IF

! CALL WRITE_DENSITY(21, CP2)
! CLOSE(21)

  998 FORMAT(2E20.8)
  999 FORMAT (2F12.6, F16.8)

  IF (NSTP /= 0) THEN
    GSTP =  XOP * G_3D / DFLOAT(NSTP)
    G = 0.0D0
    CALL NORM(CP, ZNORM) ! NORM() calculates norm and restores normalization
    CALL CHEM(CP, MU, EN) ! CHEM() calculates the chemical potential MU and energy EN.
    CP2 = CP * CONJG(CP)
    CALL RAD(CP2, RMS) ! RAD() calculates the r.m.s radius RMS
    WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RMS(1) 

    DO K = 1, NSTP       ! NSTP iterations to introduce the nonlinearity
      G = G + GSTP
      CALL CALCNU(CP,DT) ! CALCNU() performs time propagation with non-derivative parts.
      CALL LUX(CP)       ! LUX() performs the time iteration with space (x) derivative alone.
      CALL LUY(CP)       ! LUY() performs the time iteration with space (y) derivative alone.
      CALL LUZ(CP)       ! LUZ() performs the time iteration with space (z) derivative alone.
      IF (OPTION_re_im.EQ.1) CALL NORM(CP, ZNORM)
    END DO

    IF (OPTION_re_im.EQ.2) CALL NORM(CP, ZNORM)
    CALL CHEM(CP, MU, EN)
    CP2 = CP * CONJG(CP)
    CALL RAD(CP2, RMS)
    WRITE (7, 1005) ZNORM, REAL(MU)/XOP, REAL(EN)/XOP, RMS(1)
    WRITE (4, 1015) RMS(2:4)
    1005 FORMAT('NSTP iter.:', 8x, F9.4, 2F12.5, 2F11.5, F11.4)
    1015 FORMAT('NSTP iter.:', 6x, F14.5, 2F16.5)
    FLUSH(7)
    FLUSH(4)
  ELSE
    G = XOP * G_3D
  END IF

  G = XOP * G_3D
  STEP = 1 
  DO H = 1, NPAS ! NPAS iterations transient
    CALL CALCNU(CP, DT)
    CALL LUX(CP)
    CALL LUY(CP)
    CALL LUZ(CP)
    IF (OPTION_re_im.EQ.1) CALL NORM(CP, ZNORM)

    IF (MOD(H, ITER) == 0) THEN
      IF (OPTION_RE_IM.EQ.1) WRITE(ITER_FILENAME,'(A,I0,A)') 'im3d-den-', STEP, '.txt'
      IF (OPTION_RE_IM.EQ.2) WRITE(ITER_FILENAME,'(A,I0,A)') 're3d-den-', STEP, '.txt'

      OPEN(101, FILE = ITER_FILENAME)

      CP2 = CP * CONJG(CP)
      CALL WRITE_DENSITY(101, CP2)   ! WRITES DENSITY FOR EVERY ITER
      CLOSE(101)

      CALL CHEM(CP, MU, EN)
      IF (OPTION_re_im.EQ.2) CALL NORM(CP, ZNORM)
      CALL RAD(CP2,  RMS)
      WRITE (7, 1007)  STEP, ZNORM, MU/XOP, EN/XOP, RMS(1)
      WRITE (4, 1017) STEP, RMS(2:4)
      1017 FORMAT('NPAS iter.:',I6, F14.5, 2F16.5)
      1007 FORMAT('NPAS iter.:',I6,2x,F9.4, 2F12.5, 2F11.5, F11.4)
      FLUSH(7)
      FLUSH(4)

      STEP = STEP + 1

      IF (OPTION_RE_IM.EQ.1) OPEN(113, FILE = 'im3d-fin-wf.txt')
      IF (OPTION_RE_IM.EQ.2) OPEN(113, FILE = 're3d-fin-wf.txt')

      CALL WRITE_WAVE_FUNCTION(113, CP)
      CLOSE(113)

    END IF
  END DO

  IF (NRUN /= 0) THEN
    DO K = 1, NRUN ! NRUN iterations to check convergence
      CALL CALCNU(CP,DT)
      CALL LUX(CP)
      CALL LUY(CP)
      CALL LUZ(CP)
      IF (OPTION_re_im.EQ.1) CALL NORM(CP, ZNORM)
    END DO
    IF (OPTION_re_im.EQ.2) CALL NORM(CP, ZNORM)
    CALL CHEM(CP, MU, EN)
    CP2 = CP * CONJG(CP)
    CALL RAD(CP2, RMS)

    WRITE (7, 1006) ZNORM, MU/XOP, EN/XOP, RMS(1)
    WRITE (4, 1019) RMS(2:4)
    1006 FORMAT('NRUN iter.:', 8x, F9.4, 2F12.5, 2F11.5, F11.4)
    1019 FORMAT('NRUN iter.:', 6x, F14.5, 2F16.5)
    FLUSH(7)
    FLUSH(4)
  END IF

!  IF (OPTION_RE_IM.EQ.1) OPEN(13, FILE = 'im3d-fin-den2d.txt')
!  IF (OPTION_RE_IM.EQ.2) OPEN(13, FILE = 're3d-fin-den2d.txt')

 ! CALL WRITE_DENSITY(13, CP2)
 ! CLOSE(13)

  IF (OPTION_RE_IM.EQ.1) OPEN(113, FILE = 'im3d-fin-wf.txt')
  IF (OPTION_RE_IM.EQ.2) OPEN(113, FILE = 're3d-fin-wf.txt')

  CALL WRITE_WAVE_FUNCTION(113, CP)
  CLOSE(113)

  ! IF (OPTION_RE_IM.EQ.1) OPEN(31, FILE = 'im3d-den3d.txt')
  ! IF (OPTION_RE_IM.EQ.2) OPEN(31, FILE = 're3d-den3d.txt')
  ! CALL WRITE_3D(31, CP2)
  ! CLOSE(31)

  IF (OPTION_RE_IM.EQ.1) THEN 
    OPEN(11, FILE = 'im3d-den1d_z.txt')
    OPEN(12, FILE = 'im3d-den1d_x.txt')
    OPEN(13, FILE = 'im3d-den1d_y.txt')
  ELSE
    OPEN(11, FILE = 're3d-den1d_z.txt')
    OPEN(12, FILE = 're3d-den1d_x.txt')
    OPEN(13, FILE = 're3d-den1d_y.txt')
  END IF

  CALL DEN1DZ(11, CP2, DEN1Z)
  CALL DEN1DX(12, CP2, DEN1X)
  CALL DEN1DY(13, CP2, DEN1Y)
  CLOSE(11)
  CLOSE(12)
  CLOSE(13)

  CALL FREE_VARIABLES()
  IF (ALLOCATED(CP2)) DEALLOCATE(CP2)

  CALL SYSTEM_CLOCK (CLCK_COUNTS_END, CLCK_RATE,COUNT_MAX) 
  CALL CPU_TIME(T2)
  END_D = DBLE(CLCK_COUNTS_END)
  START_D = DBLE(CLCK_COUNTS_BEG)
  IF (END_D < START_D) END_D = END_D + DBLE(COUNT_MAX)
  WRITE (7, 1001)
  WRITE (4, 1001)
  CLOSE (4)
  WRITE (7,*)
  WRITE (7,'(A,I8,A)') ' Clock Time: ', INT(( END_D  -  START_D )/DBLE(CLCK_RATE)), ' seconds'
  WRITE (7,'(A,I8,A)') '   CPU Time: ', INT(T2-T1), ' seconds'
  CLOSE (7)
END PROGRAM GROSS_PITAEVSKII_SSCN_3D

SUBROUTINE ALLOCATE_VARIABLES()
  USE COMM_DATA, ONLY : NX, NY, NZ
  USE GPE_DATA, ONLY : X, Y, Z, X2, Y2, Z2, V, CP, R2
  USE CN_DATA, ONLY : CALA, CGAA, CALB, CGAB, CALC, CGAC
  IMPLICIT NONE

  ALLOCATE(X(0:NX))
  ALLOCATE(X2(0:NX))
  ALLOCATE(Y(0:NY))
  ALLOCATE(Y2(0:NY))
  ALLOCATE(Z(0:NZ))
  ALLOCATE(Z2(0:NZ))
  ALLOCATE(V(0:NX, 0:NY, 0:NZ)) 
  ALLOCATE(CP(0:NX, 0:NY, 0:NZ))
  ALLOCATE(R2(0:NX, 0:NY, 0:NZ))
  ALLOCATE(CALA(0:NX,0:NY))
  ALLOCATE(CGAA(0:NX,0:NY))
  ALLOCATE(CALB(0:NX,0:NY))
  ALLOCATE(CGAB(0:NX,0:NY))
  ALLOCATE(CALC(0:NZ))
  ALLOCATE(CGAC(0:NZ))
END SUBROUTINE ALLOCATE_VARIABLES

SUBROUTINE FREE_VARIABLES()
  USE GPE_DATA, ONLY : X, Y, Z, X2, Y2, Z2, V, CP, R2
  USE CN_DATA, ONLY : CALA, CGAA, CALB, CGAB, CALC, CGAC
  IMPLICIT NONE

  IF (ALLOCATED(X)) DEALLOCATE(X)
  IF (ALLOCATED(X2)) DEALLOCATE(X2)
  IF (ALLOCATED(Y)) DEALLOCATE(Y)
  IF (ALLOCATED(Y2)) DEALLOCATE(Y2)
  IF (ALLOCATED(Z)) DEALLOCATE(Z)
  IF (ALLOCATED(Z2)) DEALLOCATE(Z2)
  IF (ALLOCATED(V)) DEALLOCATE(V)
  IF (ALLOCATED(CP)) DEALLOCATE(CP)
  IF (ALLOCATED(R2)) DEALLOCATE(R2)
  IF (ALLOCATED(CALA)) DEALLOCATE(CALA)
  IF (ALLOCATED(CGAA)) DEALLOCATE(CGAA)
  IF (ALLOCATED(CALB)) DEALLOCATE(CALB)
  IF (ALLOCATED(CGAB)) DEALLOCATE(CGAB)
  IF (ALLOCATED(CALC)) DEALLOCATE(CALC)
  IF (ALLOCATED(CGAC)) DEALLOCATE(CGAC)
END SUBROUTINE FREE_VARIABLES

SUBROUTINE INITIALIZE()
  ! Routine that initializes the constant and variables.
  ! Calculates the potential term V and the initial wave function CP
  USE COMM_DATA, ONLY : NX, NY, NZ, NX2, NY2, NZ2, PI, NSTP
  USE GPE_DATA, ONLY : GAMMA, NU, LAMBDA, DX, DY, DZ, &
                       CP, X, X2, Y, Y2, Z, Z2, CI,R2,  OPTION_RE_IM,CI,RANDOM,FUNCTION
  IMPLICIT NONE
  REAL (8) :: TMP, TX, TY
  INTEGER :: I, J, K
  REAL (8) ::  SP, PI34
  COMPLEX (8), DIMENSION(0:NX,0:NY) :: TMP2D
  INTEGER, DIMENSION (100) :: seed   
  real (8), dimension((1+NX)*(1+NY)) :: RANDNUM

  
  PI34 = SQRT(PI*SQRT(PI))
  SP = PI34 /SQRT(SQRT(NU * LAMBDA * GAMMA))

  FORALL (I=0:NX) X(I) = (I-NX2)*DX
  FORALL (J=0:NY) Y(J) = (J-NY2)*DY
  FORALL (K=0:NZ) Z(K) = (K-NZ2)*DZ
  X2 = X*X
  Y2 = Y*Y
  Z2 = Z*Z

  DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
    R2(I,J,K)=X(I)**2+Y(J)**2+Z(K)**2
  END DO; END DO; END DO

  IF (NSTP == 0) THEN
    IF (OPTION_RE_IM.EQ.1) WRITE(*,'(a)') "Run the program using the input file to read. e.g.: ./imag3d < im3d-fin-wf.txt"
    IF (OPTION_RE_IM.EQ.2) WRITE(*,'(a)') "Run the program using the input file to read. e.g.: ./real3d < im3d-fin-wf.txt"

    DO K = 0, NZ
      DO J = 0, NY
        DO I = 0, NX
          READ (*, 998) TX, TY
          CP(I,J,K) = CMPLX(TX,TY)
        END DO
        READ(*,*)
      END DO
      READ(*,*)
    END DO
  ELSE

   seed = 13
  IF(RANDOM.EQ.1)  CALL RANDOM_SEED(PUT=seed)
  
    IF(RANDOM.EQ.1)   call random_number(RANDNUM)
 


  
    DO J = 0, NY; DO I = 0, NX
    K=(I+1)*(J+1)
    
   IF(FUNCTION.EQ.0)     TMP2D(I,J)= EXP (2.d0*PI*CI*RANDNUM(K))
    IF(FUNCTION.EQ.1)     TMP2D(I,J)=(X(I)+CI*Y(J))* EXP (2.d0*PI*CI*RANDNUM(K))
 
  END DO; END DO
    !$OMP PARALLEL DO PRIVATE(I,J,K,TMP)
  DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
      TMP = ( GAMMA * X2(I) +  NU * Y2(J) + LAMBDA * Z2(K) )/2.0D0
      CP(I,J,K) =  TMP2D(I,J)* EXP(-TMP)/SP
    END DO; END DO; END DO
    !$OMP END PARALLEL DO
  END IF
  998 FORMAT(ES16.9, 1X, ES16.9)
END SUBROUTINE INITIALIZE

SUBROUTINE CALCULATE_TRAP()
  ! Calculates the harmonic oscillator potential term V.
  USE COMM_DATA, ONLY : NX, NY, NZ
  USE GPE_DATA, ONLY : XOP, V, X2, Y2, Z2, GAMMA, NU, LAMBDA
  IMPLICIT NONE

  INTEGER :: I, J, K
  REAL (8) :: GAMMA2, NU2, LAMBDA2

  GAMMA2 = GAMMA * GAMMA
  NU2 = NU * NU
  LAMBDA2 = LAMBDA*LAMBDA

  !$OMP PARALLEL DO PRIVATE(I,J,K)
  DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
    V(I,J,K) = XOP * ( GAMMA2 * X2(I) + NU2 * Y2(J) + LAMBDA2 * Z2(K) ) / 2.0D0
  END DO; END DO; END DO
  !$OMP END PARALLEL DO
END SUBROUTINE CALCULATE_TRAP

SUBROUTINE COEF()
  ! Calculates the coefficients needed in subroutine LUX and LUY.
  USE COMM_DATA, ONLY : NXX, NYY, NZZ
  USE GPE_DATA, ONLY : DX, DY, DZ, DT, X, Y, Z, XOP, OMEGA,  OPTION_RE_IM,CI
  USE CN_DATA
  IMPLICIT NONE

  INTEGER :: I, J, K
  REAL (8) :: DX2, DY2, DZ2
  REAL (8) :: DXX, DYY, DZZ
  COMPLEX (8) :: CDT, CIJ

  IF (OPTION_RE_IM.EQ.1) CIJ = CMPLX(1.D0,0.D0)
  IF (OPTION_RE_IM.EQ.2) CIJ = CMPLX(0.D0,1.D0)

  DX2 = DX*DX   ! Generating the Coefficients for the
  DY2 = DY*DY   ! C-N method (to solve the spatial part)
  DZ2 = DZ*DZ

  DXX = 1.0D0/DX2
  DYY = 1.0D0/DY2
  DZZ = 1.0D0/DZ2
  CDT = CIJ*DT

  CA0 = 1.0D0 + CDT*DXX
  CA0R = 1.0D0 - CDT*DXX
  CB0 = 1.0D0 + CDT*DYY
  CB0R = 1.0D0 - CDT*DYY
  CC0 = 1.D0 + CDT*DZZ
  CC0R = 1.D0 - CDT*DZZ

  C0 = CIJ*CI*DT*XOP*OMEGA/2.0D0

  CTMPX = C0/(2.0D0*DX)
  CT0 = CDT*DXX/2.0D0
  CAM = -(CT0 - CTMPX*Y)
  CAP = -(CT0 + CTMPX*Y)

  DO J = 0, NY
    CALA(NXX, J) = 0.0D0
    CGAA(NXX, J) = -1.0D0/CA0
    DO I = NXX, 1, -1
      CALA(I-1,J) = CAM(J)*CGAA(I,J)
      CGAA(I-1,J) = -1.0D0/(CA0+CAP(J)*CALA(I-1,J))
    END DO
  END DO

  CTMPY = C0/(2.0D0*DY)
  CT0 = CDT*DYY/2.0D0
  CBM = -(CT0 + CTMPY*X)
  CBP = -(CT0 - CTMPY*X)

  DO I = 0, NX
    CALB(I, NYY) = 0.0D0
    CGAB(I, NYY) = -1.0D0/CB0
    DO J = NYY, 1, -1
      CALB(I,J-1) = CBM(I)*CGAB(I,J)
      CGAB(I,J-1) = -1.0D0/(CB0+CBP(I)*CALB(I,J-1))
    END DO
  END DO

  CT0Z = -CDT*DZZ/2.0D0
  CALC(NZZ) = 0.0D0
  CGAC(NZZ) = -1.0D0/CC0
  DO K = NZZ, 1, -1
    CALC(K-1) = CT0Z*CGAC(K)
    CGAC(K-1) = -1.0D0/(CC0+CT0Z*CALC(K-1))
  END DO
END SUBROUTINE COEF

SUBROUTINE CALCNU(CP, DT) ! Exact solution
  ! Solves the partial differential equation with the potential and
  ! the nonlinear term.
  USE COMM_DATA, ONLY : NX, NY, NZ, PI
  USE GPE_DATA, ONLY : OPTION, V, G, OPTION_RE_IM, CI
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT

  REAL (8), DIMENSION(0:NX,0:NY,0:NZ) :: P2,  TMP3D
  COMPLEX (8) :: CIJ
  INTEGER :: I, J, K

  IF (OPTION_RE_IM.EQ.1) CIJ=CMPLX(1.D0,0.D0)
  IF (OPTION_RE_IM.EQ.2) CIJ=CMPLX(0.D0,1.D0)

  !$OMP PARALLEL DO PRIVATE(I, J, K)
  DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
    P2(I,J,K)  = CP(I,J,K) * CONJG(CP(I,J,K))
    TMP3D(I,J,K) = (V(I,J,K) + G * P2(I,J,K))
    CP(I,J,K) = CP(I,J,K) * EXP(-CIJ * DT * TMP3D(I,J,K))
  END DO; END DO; END DO
  !$OMP END PARALLEL DO
END SUBROUTINE CALCNU

SUBROUTINE LUX(CP)
  ! Solves the partial differential equation only with the X-space
  ! derivative term using the Crank-Nicholson method
  USE OMP_LIB
  USE COMM_DATA, ONLY : NX, NY, NZ, NXX
  USE GPE_DATA, ONLY : option_re_im
  USE CN_DATA, ONLY : CAP, CAM, CA0, CA0R, CALA, CGAA
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(INOUT) :: CP
  COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ) :: CBE
  COMPLEX (8) :: CXX
  INTEGER :: I, J, K

  !$OMP PARALLEL DO PRIVATE(I, J, K, CXX)
  DO K = 0, NZ ; DO J = 0, NY
    IF (option_re_im.eq.1) CBE(NXX,J,K) = 0.d0
    IF (option_re_im.eq.2) CBE(NXX,J,K) = CP(NX,J,K)

    DO I = NXX, 1, -1
      CXX = -CAP(J) * CP(I+1,J,K) + CA0R * CP(I,J,K) - CAM(J) * CP(I-1,J,K)
      CBE(I-1,J,K) = CGAA(I,J) * (CAP(J) * CBE(I,J,K)- CXX)
    END DO
    CP(0,J,K) = 0.0D0
    DO I = 0, NXX
      CP(I+1,J,K) = CALA(I,J) * CP(I,J,K) + CBE(I,J,K)
    END DO
    CP(NX,J,K) = 0.0D0
  END DO; END DO
  !$OMP END PARALLEL DO
!-----------------------------
! Boundary condition periodic:
!      CP(0,J,K) = CP(NX,J,K)
!      CP(1,J,K) = CP(NXX,J,K)
!-----------------------------
END SUBROUTINE LUX

SUBROUTINE LUY(CP)
  ! Solves the partial differential equation only with the Y-space
  ! derivative term using the Crank-Nicholson method
  USE OMP_LIB
  USE COMM_DATA, ONLY : NX, NY, NZ, NYY
  USE CN_DATA, ONLY : CBP, CBM, CB0, CB0R, CALB, CGAB
  USE GPE_DATA, ONLY : option_re_im
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(INOUT) :: CP
  COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ) :: CBE
  COMPLEX (8) :: CYY
  INTEGER :: I, J, K

  !$OMP PARALLEL DO PRIVATE(I, J, K, CYY)
  DO K = 0, NZ; DO I = 0, NX
    IF (option_re_im.eq.1) CBE(I,NYY,K) = 0.D0
    IF (option_re_im.eq.2) CBE(I,NYY,K) = CP(I,NY,K)

    DO J = NYY, 1, -1
      CYY = -CBP(I) * CP(I,J+1,K) + CB0R * CP(I,J,K) - CBM(I) * CP(I,J-1,K)
      CBE(I,J-1,K) = CGAB(I,J) * (CBP(I) * CBE(I,J,K) - CYY)
    END DO
    CP(I,0,K) = 0.0D0
    DO J = 0, NYY
      CP(I,J+1,K) = CALB(I,J) * CP(I,J,K) + CBE(I,J,K)
    END DO
    CP(I,NY,K) = 0.0D0
  END DO; END DO
  !$OMP END PARALLEL DO
!-----------------------------
! Boundary condition periodic:
!        CP(I,0,:) = CP(I,NY,:)
!        CP(I,1,:) = CP(I,NYY,:)
!-----------------------------
END SUBROUTINE LUY

SUBROUTINE LUZ(CP)
  ! Solves the partial differential equation only with the Z-space
  ! derivative term using the Crank-Nicholson method
  USE OMP_LIB
  USE COMM_DATA, ONLY : NX, NY, NZ, NZZ
  USE CN_DATA, ONLY : CC0R, CT0Z, CALC, CGAC
  USE GPE_DATA, ONLY : option_re_im
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(INOUT) :: CP
  COMPLEX (8) :: CBE(0:NX,0:NY,0:NZ)
  COMPLEX (8) :: CZZ
  INTEGER :: I, J, K

  !$OMP PARALLEL DO PRIVATE(I, J, K, CZZ)
  DO J = 0, NY; DO I = 0, NX
    IF (option_re_im.eq.1) CBE(I, J, NZZ) = 0.d0
    IF (option_re_im.eq.2) CBE(I, J, NZZ) = CP(I, J, NZ)

    DO K = NZZ, 1, -1
      CZZ = -CT0Z * CP(I,J,K+1) + CC0R * CP(I,J,K) - CT0Z * CP(I,J,K-1)
      CBE(I, J, K-1) = CGAC(K) * (CT0Z * CBE(I,J,K) - CZZ)
    END DO
    CP(I,J,0) = 0.0D0
    DO K = 0, NZZ
      CP(I,J,K+1) = CALC(K) * CP(I,J,K) + CBE(I,J,K)
    END DO
    CP(I,J,NZ) = 0.0D0
  END DO; END DO
  !$OMP END PARALLEL DO
!-----------------------------
! Boundary condition periodic:
! DO K = 0, NZZ
! CP(I,J,K+1) = CALC(K)*CP(I,J,K) + CBE(K)
! END DO
! CP(I,J,0) = CP(I,J,NZ)
! CP(I,J,1) = CP(I,J,NZZ)
!-----------------------------
END SUBROUTINE LUZ

SUBROUTINE NORM(CP, ZNORM)
  ! Calculates the normalization of the wave function and sets it to unity.
  USE OMP_LIB
  USE GPE_DATA, ONLY : DX, DY, DZ, option_re_im
  USE COMM_DATA, ONLY : NX, NY, NZ
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: ZNORM

  INTERFACE
    FUNCTION INTEGRATE(P2, DX, DY,DZ) RESULT(RES)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: P2
      REAL (8), INTENT(IN) :: DX, DY, DZ
      REAL (8) :: RES
    END FUNCTION INTEGRATE
  END INTERFACE

  REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: P2
  INTEGER :: I, J, K

  !$OMP PARALLEL DO PRIVATE(I, J, K)
  DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
    P2(I,J,K) = CP(I,J,K) * CONJG(CP(I,J,K))
  END DO; END DO; END DO
  !$OMP END PARALLEL DO

  ZNORM = SQRT(INTEGRATE(P2, DX, DY, DZ))

  IF (option_re_im.EQ.1) THEN
    !$OMP PARALLEL DO PRIVATE(I, J, K)
    DO K = 0, NZ; DO J = 0, NY;  DO I = 0, NX
      CP(I,J,K) = CP(I,J,K)/ZNORM
    END DO; END DO; END DO
    !$OMP END PARALLEL DO
  END IF
END SUBROUTINE NORM

SUBROUTINE NORMCHEM(CPR, ZNORM)
  ! Calculates the normalization of the wave function and sets it to unity.
  USE OMP_LIB
  USE GPE_DATA, ONLY : DX, DY, DZ
  USE COMM_DATA, ONLY : NX, NY, NZ
  IMPLICIT NONE
  REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(INOUT) :: CPR
  REAL (8), INTENT(OUT) :: ZNORM

  INTERFACE
    FUNCTION INTEGRATE(P2, DX, DY,DZ) RESULT(RES)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: P2
      REAL (8), INTENT(IN) :: DX, DY, DZ
      REAL (8) :: RES
    END FUNCTION INTEGRATE
  END INTERFACE

  REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: P2
  INTEGER :: I, J, K

  !$OMP PARALLEL DO PRIVATE(I, J, K)
  DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
    P2(I,J,K) = CPR(I,J,K) * CPR(I,J,K)
  END DO; END DO; END DO
  !$OMP END PARALLEL DO

  ZNORM = INTEGRATE(P2, DX, DY, DZ)
END SUBROUTINE NORMCHEM

SUBROUTINE RAD(CP2, R)
  ! Calculates the root mean square size RMS
  USE COMM_DATA, ONLY : NX, NY, NZ
  USE GPE_DATA, ONLY : DX, DY, DZ, X2, Y2, Z2, R2
  IMPLICIT NONE
  REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: CP2
  REAL (8), DIMENSION(:), INTENT(OUT) :: R

  INTERFACE
    FUNCTION INTEGRATE(U, DX, DY, DZ) RESULT(VALUE)
      USE COMM_DATA, ONLY : NX, NY, NZ
      REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: U
      REAL (8), INTENT (IN) :: DX, DY, DZ
      REAL (8) :: VALUE
    END FUNCTION INTEGRATE
  END INTERFACE

  INTEGER :: I, J, K
  REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: TMP3D

  TMP3D = R2*CP2
  R(1) = SQRT(INTEGRATE(TMP3D, DX, DY, DZ))

  FORALL(K=0:NZ, J=0:NY) TMP3D(:,J,K) = X2 * CP2(:,J,K)
  R(2) = SQRT(INTEGRATE(TMP3D, DX, DY, DZ))

  FORALL(K=0:NZ, I=0:NX) TMP3D(I,:,K) = Y2 * CP2(I,:,K)
  R(3) = SQRT(INTEGRATE(TMP3D, DX, DY, DZ))

  FORALL(J=0:NY, I=0:NX) TMP3D(I,J,:) = Z2 * CP2(I,J,:)
  R(4) = SQRT(INTEGRATE(TMP3D, DX, DY, DZ))
END SUBROUTINE RAD

SUBROUTINE CHEM(CP, MU, EN)
  ! Calculates the chemical potential MU and energy EN.  CP is the wave
  ! function, V is the potential and G is the nonlinearity.
  USE COMM_DATA, ONLY : NX, NY, NZ
  USE GPE_DATA, ONLY : DX, DY, DZ, V, G, OMEGA, X, Y, CI,XOP
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: CP
  REAL (8), INTENT(OUT) :: MU, EN

  INTERFACE
    PURE FUNCTION DIFF(P, DX) RESULT (DP)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: P
      REAL (8), INTENT(IN) :: DX
      REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
    END FUNCTION DIFF

    FUNCTION INTEGRATE(P2, DX, DY, DZ) RESULT(RES)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: P2
      REAL (8), INTENT(IN) :: DX, DY, DZ
      REAL (8) :: RES
    END FUNCTION INTEGRATE

    SUBROUTINE NORMCHEM(CPR, ZNORM)
      USE COMM_DATA, ONLY : NX, NY, NZ
      IMPLICIT NONE
      REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(INOUT) :: CPR
      REAL (8), INTENT(OUT) :: ZNORM
    END  SUBROUTINE NORMCHEM
  END INTERFACE

  INTEGER :: I, J, K
  REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: P2, GP2, DP2, CPR,CPI,P2R
  REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: TMP3D, EMP3D, DPLZ
  REAL (8), DIMENSION(0:NX, 0:NY,0:NZ) :: DPXR, DPYR, DPXI,DPYI, DPZR
  REAL (8) :: ZNORM

  MU = 0.0D0
  EN = 0.0D0

  CPR=CP
  CPI=AIMAG(CP)

  DO K = 0, NZ; DO I = 0, NX
    DPYR(I,0:NY,K) = DIFF(CPR(I,0:NY,K), DY)
    DPYI(I,0:NY,K) = DIFF(CPI(I,0:NY,K), DY)
  END DO; END DO

  DO K = 0, NZ; DO J = 0, NY
    DPXR(0:NX,J,K) = DIFF(CPR(0:NX,J,K), DX)
    DPXI(0:NX,J,K) = DIFF(CPI(0:NX,J,K), DX)
  END DO; END DO

  DO J = 0, NY; DO I = 0, NX
    DPZR(I,J,0:NZ) = DIFF(CPR(I,J,0:NZ), DZ)
  END DO; END DO

  DP2 = DPXR*DPXR + DPYR*DPYR + DPZR*DPZR

  DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
    DPLZ(I,J,K) = CPR(I,J,K) * ( X(I)*DPYI(I,J,K) - Y(J)*DPXI(I,J,K) )
  END DO; END DO; END DO

  P2 = CP*CONJG(CP)
  P2R=CPR*CPR
  GP2 = G*P2
  TMP3D = (V + GP2) * P2R + DP2 - XOP * (OMEGA * DPLZ)
  EMP3D = (V + 0.5D0 * GP2) * P2R + DP2 - XOP * (OMEGA * DPLZ)

  CALL NORMCHEM(CPR,ZNORM)

  MU = INTEGRATE(TMP3D, DX, DY, DZ)/ZNORM
  EN = INTEGRATE(EMP3D, DX, DY, DZ)/ZNORM
END SUBROUTINE CHEM

FUNCTION INTEGRATE(U, DX, DY, DZ) RESULT(RES)
  USE COMM_DATA, ONLY : NX, NY, NZ
  IMPLICIT NONE
  REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: U
  REAL (8), INTENT (IN) :: DX, DY, DZ
  REAL (8) :: RES

  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE

  REAL (8), DIMENSION(0:NX) :: TMPX
  REAL (8), DIMENSION(0:NY) :: TMPY
  REAL (8), DIMENSION(0:NZ) :: TMPZ
  INTEGER :: I, J, K

  !$OMP PARALLEL DO PRIVATE(I, J, TMPX, TMPY)
  DO K = 0, NZ
    DO J = 0, NY
      DO I = 0, NX
        TMPX(I) = U(I,J,K)
      END DO
      TMPY(J) = SIMP(TMPX, DX)
    END DO
    TMPZ(K) = SIMP(TMPY, DY)
  END DO
  !$OMP END PARALLEL DO
  RES = SIMP(TMPZ, DZ)
END FUNCTION INTEGRATE

PURE FUNCTION SIMP(F, DX)
  ! Does the spatial integration with Simpson's rule.
  ! N refer to the number of integration points, DX space step, and
  ! F is the function to be integrated.
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: F
  REAL (8), INTENT(IN) :: DX
  REAL (8) :: SIMP
  REAL (8) :: F1, F2
  INTEGER :: I, N

  N = SIZE(F) - 1
  F1 = F(1) + F(N-1) ! N EVEN
  F2 = F(2)
  DO I = 3, N-3, 2
     F1 = F1 + F(I)
     F2 = F2 + F(I+1)
  END DO
  SIMP = DX*(F(0) + 4.D0*F1 + 2.D0*F2 + F(N))/3.D0
END FUNCTION SIMP

PURE FUNCTION DIFF(P,DX) RESULT (DP)
  ! Computes the first derivative DP of P using
  ! Richardsonextrapolation formula. The derivative at the
  ! boundaries are assumed to be zero
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: P
  REAL (8), INTENT(IN) :: DX
  REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
  INTEGER :: I, N

  N = SIZE(P) - 1
  DP(0) = 0.D0
  DP(1) = (P(2) - P(0))/(2.D0*DX)
  FORALL(I=2:N-2)
    DP(I) = (P(I-2)-8.D0*P(I-1)+8.D0*P(I+1)-P(I+2))/(12.D0*DX)
  END FORALL
  DP(N-1) = (P(N) - P(N-2))/(2.D0*DX)
  DP(N) = 0.D0
END FUNCTION DIFF

SUBROUTINE WRITE_WAVE_FUNCTION(FUNIT, CP)
  ! Writes the value of wave function (poth real and imaginary parts)
  ! at every spatial point.
  USE COMM_DATA, ONLY : NX, NY, NZ
  INTEGER, INTENT(IN) :: FUNIT
  COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: CP
  INTEGER :: I, J, K

  DO K = 0, NZ
    DO J = 0, NY
      DO I = 0, NX
        WRITE (FUNIT, 998) REAL(CP(I,J,K)), AIMAG(CP(I,J,K))
      END DO
      WRITE(FUNIT,*)
    END DO
    WRITE(FUNIT,*)
  END DO
  998 FORMAT(ES16.9, 1X, ES16.9)
END SUBROUTINE WRITE_WAVE_FUNCTION

SUBROUTINE WRITE_3D(FUNIT, U2)
  ! Writes the 3D density at every spatial point.
  USE COMM_DATA, ONLY : NX, NY, NZ
  INTEGER, INTENT(IN) :: FUNIT
  REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: U2
  INTEGER :: I, J, K

  DO K = 0, NZ; DO J = 0, NY; DO I = 0, NX
    WRITE(FUNIT, 1000) U2(I,J,K)
  END DO; END DO; END DO
  1000 FORMAT(E17.6E3)
END SUBROUTINE WRITE_3D

SUBROUTINE WRITE_DENSITY(FUNIT, U2)
  ! Writes the 2D density cross-section at Z = 0.
  ! The format is (Y,X,den(X,Y)).
  USE COMM_DATA, ONLY : NX, NY, NZ, NZ2
  USE GPE_DATA, ONLY : X, Y
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FUNIT
  REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: U2

  INTEGER :: I, J

  DO J = 0, NY
     DO I = 0, NX
        WRITE(FUNIT, 999) Y(J), X(I), U2(I,J,NZ2)
     END DO
     WRITE(FUNIT, *)
  END DO
  999 FORMAT(ES13.6, 1X, ES13.6, 1X, ES13.6)
END SUBROUTINE WRITE_DENSITY

SUBROUTINE DEN1DX(FUNIT, U2, DE1DX)
  ! Writes the integrated 1D density over X dimension.
  USE COMM_DATA, ONLY : NX, NY, NZ
  USE GPE_DATA, ONLY : X, DY, DZ
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FUNIT
  REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: U2
  REAL (8), DIMENSION(0:NX), INTENT(OUT) :: DE1DX
  REAL (8), DIMENSION(0:NX, 0:NZ) :: TMP2D

  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE

  INTEGER :: I, K

  FORALL (K = 0:NZ, I = 0:NX) TMP2D(I,K) = SIMP(U2(I,0:,K), DY)
  DO I = 0, NX
    DE1DX(I) = SIMP(TMP2D(I,0:), DZ)
    WRITE(FUNIT, 1001) X(I), DE1DX(I)
  END DO
  1001 FORMAT(ES13.6, 1X, ES13.6)
END SUBROUTINE DEN1DX

SUBROUTINE DEN1DY(FUNIT, U2, DE1DY)
  ! Writes the integrated 1D density over Y dimension.
  USE COMM_DATA, ONLY : NX, NY, NZ
  USE GPE_DATA, ONLY : Y, DX, DZ
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FUNIT
  REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: U2
  REAL (8), DIMENSION(0:NY), INTENT(OUT) :: DE1DY
  REAL (8), DIMENSION(0:NY, 0:NZ) :: TMP2D

  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE

  INTEGER :: J, K

  FORALL (K = 0:NZ, J = 0:NY) TMP2D(J,K) = SIMP(U2(0:,J,K), DX)
  DO J = 0, NY
    DE1DY(J) = SIMP(TMP2D(J,0:), DZ)
    WRITE(FUNIT, 1001) Y(J), DE1DY(J)
  END DO
  1001 FORMAT(ES13.6, 1X, ES13.6)
END SUBROUTINE DEN1DY

SUBROUTINE DEN1DZ(FUNIT, U2, DE1DZ)
  ! Writes the integrated 1D density over Z dimension.
  USE COMM_DATA, ONLY : NX, NY, NZ
  USE GPE_DATA, ONLY : Z, DX, DY
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FUNIT
  REAL (8), DIMENSION(0:NX,0:NY,0:NZ), INTENT(IN) :: U2
  REAL (8), DIMENSION(0:NZ), INTENT(OUT) :: DE1DZ
  REAL (8), DIMENSION(0:NY, 0:NZ) :: TMP2D

  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE

  INTEGER :: J, K

  FORALL (K = 0:NZ, J = 0:NY) TMP2D(J,K) = SIMP(U2(0:,J,K), DX)
  DO K = 0, NZ
    DE1DZ(K) = SIMP(TMP2D(0:,K), DY)
    WRITE(FUNIT, 1001) Z(K), DE1DZ(K)
  END DO
  1001 FORMAT(ES13.6, 1X, ES13.6)
END SUBROUTINE DEN1DZ
