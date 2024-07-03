
module mod_global_parameters

  use mod_precision
  implicit none
  integer  :: max_nl_iter, num_extrap, glob_iter, max_glob_iter
  integer  :: num_sig_sv
  real(dp) :: len_conv, nl_iter_tol, sigma, unif_ini_temp, z_coord_shft
  real(dp) :: pconv, rbar, gsubc, jfac
  ! room dimensions: H = height, W = width, D = depth, pipes in D direction
  real(dp) :: rm_h, rm_w, rm_d
  character(len=3) :: anal_ini_prof, verif, debug
  ! units   = problem units,'si' and 'engineering' are valid key words
  character(len=11) :: units

end module mod_global_parameters
