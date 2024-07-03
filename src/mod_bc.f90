
module mod_bc
  ! things related to boundary conditions
  use mod_precision
  implicit none

  ! convection boundary conditions
  integer                                     :: no_conv_bc
  integer,          parameter                 :: max_no_conv_bc = 5
  integer,          dimension(max_no_conv_bc) :: ne_bc_conv, ni_bc_conv, nj_bc_conv
  real(dp),         dimension(max_no_conv_bc) :: t_inf_cbc 
  character(len=7), dimension(max_no_conv_bc) :: h_geom
  character(len=3)                            :: first_time

  ! specified temperature boundary conditions
  integer                                     :: no_spec_t_bc
  integer,          parameter                 :: max_no_st_bc = 10
  integer,          dimension(max_no_st_bc)   :: ni_bc_spec_t
  real(dp),         dimension(max_no_st_bc)   :: t_st_bc 

end module mod_bc
