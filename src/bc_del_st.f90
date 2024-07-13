
subroutine bc_del_st(node_no)
  ! specified temperature bc in delta form
  ! node_no = node no for which temp is known bc
  use mod_precision
  use mod_constants
  use mod_nodes_elements
  implicit none
  integer,        intent(in) :: node_no
  integer                    :: i
  i = node_no
  ! zero row i of global matrix [C]
  c(i,1:nfb) = zero
  ! set diagonal of row i to unity
  c(i,nhb+1) = one
  ! RHS in delta form is zero
  b(i) = zero

end subroutine bc_del_st
