

complex function c_prop(jlo,f,del_t,p,tt)
  ! take f and jlo from subroutine find, perform the interpolation step and
  ! convert the interpolated quantity to a complex variable
  ! p(jlo) <= p <= p(jlo + 1), dp/dT is sensitivity of p to T
  ! jlo   = index in tabular data where T was found in subroutine find
  ! f     = interpolation fraction 0 < f < 1
  ! del_t = size of complex step
  ! p(:)  = array dependent variables
  use mod_precision
  use mod_constants
  implicit none
  integer,  intent(in)               :: jlo
  real(dp), intent(in)               :: f, del_t
  real(dp), intent(in), dimension(*) :: p, tt
  integer                            :: jlop1
  real(dp)                           :: prop, dp_dt

  jlop1 = jlo + 1
  prop = p(jlo)*(one - f) + p(jlop1)*f              ! linear interpolation
  dp_dT = (p(jlop1) - p(jlo))/(tt(jlop1) - tt(jlo)) ! sensitivity of p to T
  c_prop = cmplx(prop,del_t*dp_dt,dp)               ! convert prop to complex

end function c_prop
