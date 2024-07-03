
MODULE MOD_VERIFY

  ! this module (almost) matches what is in the three modules
  ! mod_precision.f90, mod_constants.f90 and mod_files.f90
  ! this was done to better isolate verify_lib.f90 from the 
  ! calling program. This makes it more self contained.

  IMPLICIT NONE
  SAVE
  ! select working precision (WP)
  INTEGER, PARAMETER :: WP = SELECTED_REAL_KIND(15,30)
  ! verification file to match calling program
  INTEGER, PARAMETER ::  KVER = 88
  ! define general constants that are available for all routines

  REAL(KIND=WP), PARAMETER :: ZERO = 0.0_WP, HALF = 0.5_WP, ONE = 1.0_WP, &
       TWO = 2.0_WP, THREE = 3.0_WP, FOUR = 4.0_WP, FIVE = 5.0_WP, &
       SIX = 6.0_WP, SEVEN = 7.0_WP, EIGHT = 8._WP, NINE = 9.0_WP, &
       TEN = 10._WP, TWELVE = 12._WP
  REAL(KIND=WP), PARAMETER :: SMALL = 1.0D-06, V_SMALL = 1.0D-12
  REAL(KIND=WP), PARAMETER :: &
       PI = 3.141592653589793238462643383279502884197_WP,&
       PIO2 = PI/2._WP, TWOPI = TWO*PI, FOURPI = FOUR*PI
  REAL(KIND=WP), PARAMETER :: &
       SQRT2 = 1.41421356237309504880168872420969807856967_WP

END MODULE MOD_VERIFY


SUBROUTINE VERIFICATION(NUM_NODE,GEOM_TYP,N_TIME,TIME,R,X,Z,T,SDOT)

  ! driver subroutine for steady conduction, no source, 
  ! k vs. T is composed of two piecewise linear segments 
  ! verification problem

  ! input:
  ! NUM_NODE = number of nodes
  ! GEOM_TYP = planar, cylindrical or spherical
  ! N_TIME   = time step number
  ! TIME     = problem time (dummy variable for steady state problem)
    ! R        = radial coordinate (may not be important for planar)
    ! X        = array of x-coordinates
    ! Z        = array of z-coordinates
    ! T        = array of nodal temperatures from numerical solution

  USE MOD_VERIFY
  IMPLICIT NONE
  interface
     subroutine anal_klin(geom_typ,l,x,r,r_l,r_r,temp)
       use mod_verify
       implicit none
       character(len=11)          :: geom_typ
       real(kind=wp), intent(in)  :: l, x, r, r_l, r_r
       real(kind=wp), intent(out) :: temp
     end subroutine anal_klin
  end interface
  INTEGER, INTENT(IN) :: NUM_NODE, N_TIME
  REAL(KIND=WP),  OPTIONAL, INTENT(IN) :: TIME, SDOT
  REAL(KIND=WP), DIMENSION(*), INTENT(IN) :: R, T, X, Z
  CHARACTER (LEN=11), INTENT(IN) :: GEOM_TYP

  ! local variables

  INTEGER :: I
  REAL(KIND=WP) :: H, L_2_NORM, SUM
  INTEGER, PARAMETER :: MAX_NUM_NODE=101
  REAL(KIND=WP), DIMENSION(MAX_NUM_NODE) :: TEMP
  WRITE(KVER,"('ZONE')")
  WRITE(KVER,"('#  I      X(I)            TEMP_anal             TEMP_num&
       &           |ERROR|')")

  ! TEMP = array of analytical solution at computational nodes

  DO I=1,NUM_NODE
     CALL ANAL_KLIN(GEOM_TYP,X(NUM_NODE)-X(1),X(I),R(I),R(1),&
          R(NUM_NODE),TEMP(I))
     WRITE(KVER,"(I4,1P,E13.5,2E23.15,E13.5)") I,X(I),TEMP(I),T(I), &
          ABS(TEMP(I)-T(I))
  END DO

  ! calculate L_2 norm of temperature error

  SUM = 0.0_WP
  DO I=1,NUM_NODE
     SUM = SUM + (TEMP(I)-T(I))**2
  END DO
  L_2_NORM = SQRT(SUM/FLOAT(NUM_NODE))       ! L_2 norm of temperature error
  H = (X(NUM_NODE) - X(1))/FLOAT(NUM_NODE-1) ! characteristic mesh size
  WRITE(KVER,"('# verification: steady conduction, k = C, h(T)')")
  WRITE(KVER,"('ZONE')")
  WRITE(KVER,"('# grid convergence'/'# num elem       h         L_2 norm')")
  WRITE(KVER,"(4X,I3,1P,3x,2E13.5)") NUM_NODE-1,H,L_2_NORM

END SUBROUTINE VERIFICATION


SUBROUTINE ANAL_KLIN(GEOM_TYP,L,X,R,R_L,R_R,TEMP)

  ! evaluation of steady state temperature profile in planar slab 
  ! with constant k
  ! this solution arose from convection bc on front face 
  ! h = a*(T_inf - T_w)**1/4, a = 12, T_inf = 1300
  ! and specified T_w = 300 on back face
  ! the analytical solution is T linear in x

  ! input:
  ! x = position to evaluate temperature

  ! output:
  ! TEMP = temperature at position x

  USE MOD_VERIFY
  IMPLICIT NONE
  CHARACTER(LEN=11)          :: GEOM_TYP
  REAL(KIND=WP), INTENT(IN)  :: L, X, R, R_L, R_R
  REAL(KIND=WP), INTENT(OUT) :: TEMP
  ! local variables
  INTEGER :: I
  REAL(KIND=WP)               :: T_L, T_R, THETA, THETA_L, THETA_R
  REAL(KIND=WP)               :: BETA, XOL
  REAL(KIND=WP), DIMENSION(2) :: K, T
  REAL(KIND=WP) :: dummy
 
 ! set some parameters
  T_L = 675.E0_WP      ! temperature on left boundary
  T_R = 300.E0_WP     ! temperature on right boundary

  ! analytical solution for constant property planar slab

  temp = T_l + (T_r - T_l)*x/L 

END SUBROUTINE ANAL_KLIN
