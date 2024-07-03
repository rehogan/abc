
subroutine el_matrix_to_global_matrix(ne, nnpe, el_jac, res)

  ! ne    = element no
  ! nnpe  = no nodes per element, 3 for triangle, 4 for quad

  use mod_precision
  use mod_constants
  use mod_nodes_elements
  implicit none
  integer,        intent(in)                       :: ne, nnpe
  real(kind=dp),  intent(in), dimension(nnpe,nnpe) :: el_jac
  real(kind=dp),  intent(in), dimension(nnpe)      :: res
  integer                                          :: i, ii, j, jj, kk, mr

  mr(i,j,nhbp1) = j + nhbp1 - i   ! local to global no transformation

  do ii = 1 , nnpe         ! ii = local node no
     i = nn(ii,ne)         ! i  = global node no for local node no ii and element no ne
     b(i) = b(i) - res(ii) ! residual contribution
     do jj = 1 , nnpe
        j = nn(jj,ne)
        kk = mr(i,j,nhbp1)
        c(i,kk) = c(i,kk) + el_jac(ii,jj)  ! global matrix contribution
     end do
  end do

end subroutine el_matrix_to_global_matrix
