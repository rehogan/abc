
subroutine el_bc_matrix_to_global_matrix(nnpe, nn_b, el_jac, res)

  ! ne    = element no
  ! nnpe  = no nodes per element, 2 for boundary condition

  use mod_precision
  use mod_constants
  use mod_nodes_elements
  implicit none
  integer,        intent(in)                       :: nnpe
  integer,        intent(in), dimension(nnpe)      :: nn_b
  real(kind=dp),  intent(in), dimension(nnpe,nnpe) :: el_jac
  real(kind=dp),  intent(in), dimension(nnpe)      :: res
  integer                                          :: i, ii, j, jj, kk, mr

  mr(i,j,nhbp1) = j + nhbp1 - i   ! local to global no transformation

  ! merge elem bc contribution into global matrix
  ! the row index stays the same but the column index is shifted
  ! because of the banded storage structure of the [C] matrix
  do ii = 1,2         ! ii = local node no
     i = nn_b(ii)     ! node no of bc node
     b(i) = b(i) - res(ii) ! residual contribution
     do jj = 1,2
        j = nn_b(jj)
        kk = mr(i,j,nhbp1)
        c(i,kk) = c(i,kk) + el_jac(ii,jj)  ! global matrix contribution
     end do
  end do

end subroutine el_bc_matrix_to_global_matrix
