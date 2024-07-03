
module mod_lin_interp
  implicit none

contains
real function lin_interp_fun(f, x_j, x_jp1)
  use mod_precision
  use mod_constants
  implicit none
  real(dp), intent(in) :: f, x_j, x_jp1
  lin_interp_fun = x_j*(one - f) + x_jp1*f
end function lin_interp_fun

end module mod_lin_interp
