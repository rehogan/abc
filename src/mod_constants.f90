module mod_constants

  ! define general constants that are available for all routines

  use mod_precision
  implicit none
  save

  real(dp), parameter :: &
       zero   = 0.0_dp,   pt125 = 0.125_dp, pt25  = 0.25_dp, &
       pt375  = 0.375_dp, pt5   = 0.5_dp,   pt75  = 0.75_dp, & 
       one    = 1.0_dp,  &
       two    = 2.0_dp,   three = 3.0_dp,   four  = 4.0_dp,  &
       five   = 5.0_dp,   six   = 6.0_dp,   seven = 7.0_dp,  &
       eight  = 8._dp,    nine  = 9.0_dp,   ten   = 10._dp,  &
       twelve = 12._dp

  real(dp), parameter :: &
       small = 1.0e-06_dp, v_small = 1.0e-12_dp, v_large = 1.e20_dp

  real(dp), parameter :: &
       pi = 3.141592653589793238462643383279502884197_dp,&
       pio2    = pi/2._dp, twopi = two*pi, fourpi = four*pi, &
       eightpi = eight*pi

  real(dp), parameter :: &
       sqrt2 = 1.41421356237309504880168872420969807856967_dp

  ! Stefan-Boltzman, Thermal Radiation, Howell, Menguc, and Daun
  ! Btu/hr-ft^2-R^4
  real(dp), parameter :: sigma_e  = 0.17123e-08_dp ! engr
  ! Stefan-Boltzman, W/m^2-K
  real(dp), parameter :: sigma_si = 5.6704e-08_dp  ! SI

  ! acceleration of gravity
  real(dp), parameter :: g_eng = 32.174_dp      ! ft/s^2
  real(dp), parameter :: g_si  = 9.80665_dp     ! m/s^2, Myers

  ! conversion factors
  ! convective heat transfer coefficient: 1 W/m^2-K = 0.17611 Btu/hr-ft^2-F
  real(dp), parameter :: conv_htc = 0.17611_dp
  ! length:   1 ft = 0.3048 m
  real(dp), parameter :: conv_ft_to_m = 0.3048_dp
  ! temp:  T(K) = T(C) + 273.15
  real(dp), parameter :: conv_c_to_k = 273.15_dp
  ! temp:  T(R) = T(F) + 459.67
  real(dp), parameter :: conv_f_to_r = 459.67_dp

end module mod_constants
