
module MOD_VERIFY

  ! this module (almost) matches what is in the three modules
  ! mod_precision.f90, mod_constants.f90 and mod_files.f90
  ! this was done to better isolate verify_lib.f90 from the 
  ! calling program. This makes it more self contained.

  implicit none
  save
  ! select working precision (WP)
  integer, parameter :: WP = selected_real_kind(15,30)
  ! verification file to match calling program
  integer, parameter ::  KVER = 88
  ! define general constants that are available for all routines

  real(KIND=WP), parameter :: ZERO = 0.0_WP, HALF = 0.5_WP, ONE = 1.0_WP, &
       TWO = 2.0_WP, THREE = 3.0_WP, FOUR = 4.0_WP, FIVE = 5.0_WP, &
       SIX = 6.0_WP, SEVEN = 7.0_WP, EIGHT = 8._WP, NINE = 9.0_WP, &
       TEN = 10._WP, TWELVE = 12._WP
  real(KIND=WP), parameter :: SMALL = 1.0D-06, V_SMALL = 1.0D-12
  real(KIND=WP), parameter :: &
       PI = 3.141592653589793238462643383279502884197_WP,&
       PIO2 = PI/2._WP, TWOPI = TWO*PI, FOURPI = FOUR*PI
  real(KIND=WP), parameter :: &
       SQRT2 = 1.41421356237309504880168872420969807856967_WP

end module MOD_VERIFY


subroutine VERIFICATION(NUM_NODE,GEOM_TYP,N_TIME,TIME,R,X,Z,T,SDOT)

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

  
  use mod_precision
  use mod_constants
  implicit none
  integer,   intent(in)               :: num_node, n_time
  real(dp),  intent(in), optional     :: time, sdot
  real(dp),  intent(in), dimension(*) :: r, t, x, z
  character(len=11), intent(in)     :: geom_typ

  ! local variables

  integer :: I
  real(DP)               :: L_2_NORM, msq 
  real(DP), dimension(2) :: TEMP

  ! TEMP = array of analytical solution at computational nodes


  ! calculate L_2 norm of temperature error for enclosure

  write(6,"(' iteration temps', 2es13.5)") t(11),t(13)
  temp(1) = 373.3712021228251_dp       ! from analytical solution
  temp(2) = 880.7359878695708_dp       ! from analytical solutioin
  msq = pt5*((t(11) - temp(1))**2 + (t(13) - temp(2))**2)
  l_2_norm = sqrt(msq)
  write(6,"('# verification: steady conduction with rad gap')")
  write(6,"('ZONE')")
  write(6,"('# iterative convergence l_2_norm')")
  write(6,"(4x,es13.5)") l_2_norm

end subroutine VERIFICATION

