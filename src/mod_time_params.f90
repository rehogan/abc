module mod_time_params

  ! things related to time and time integration

  use mod_precision
  use mod_constants
  implicit none
  integer                :: n_time, first_index, last_index
  integer                :: nm1, n, np1 ! 3 levels of time, (n-1) (n) (n+1)
  real(dp)               :: time, time_ini, time_fin, theta, max_corr
  real(dp), dimension(3) :: dtime       ! 3 time steps, @ (n-1) (n) (n+1) 
  character(len=12)      :: time_integ  ! either implicit or trapezoidal
  character(len=3)       :: ss          ! ss = yes, no transient
end module mod_time_params
