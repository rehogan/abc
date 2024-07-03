
function mat_vec_mult(A, v) result (w)
  ! Matrix-Vector multiplication, adapted from
  ! http://www.cs.emory.edu/~cheung/Courses/561/Syllabus/6-Fortran/Progs/operator-func03.f90
  use mod_precision
  implicit none

  real(kind=dp), dimension(:,:), intent(in) :: A
  real(kind=dp), dimension(:), intent(in)   :: v
  real(kind=dp), dimension( SIZE(A,1) )     :: w

  integer :: i, j
  integer :: N

  N = size(v)

  w = 0.0_dp       !! clear whole vector
  do i=1,n
     do j=1,n
        w(i) = w(i) + a(i,j)*v(j)
     end do
  end do

end function mat_vec_mult
