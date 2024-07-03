
subroutine bc_del_nl_flux(b_el_no,nn_b,a,exp,t_inf)
  ! non-linear flux of the form q_dot^" = a(T_inf - T_w)**m
  ! this form is prompted by natural convection correlations
  ! BE CAREFUL: the units can get screwed up because of sloppiness 
  ! should start with diminsionless correlations so this is not
  ! an issue

  ! b_el_no = boundary elem no
  ! nn_b()  = node no that define side length S_(i,j) for this b_el_no
  ! T_inf   = temp to which surf convects to
  ! a       = coefficient with some screwy units
  ! exp     = exponent in correlation
  ! t_inf   = free stream temperature for convection

  use mod_precision
  use mod_constants
  use mod_nodes_elements
  implicit none
  integer,       intent(in)               :: b_el_no
  integer,       intent(in), dimension(2) :: nn_b
  real(kind=dp), intent(in)               :: a, exp, t_inf
  real(kind=dp)                 :: tb, t_i, t_j, sij, coef
  real(kind=dp), dimension(2,2) :: jac_nl_flux
  real(kind=dp), dimension(2)   :: q_dot, res
  integer                       :: i, j, ii, jj, kk, mr

  interface
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
     end subroutine el_bc_matrix_to_global_matrix

    end interface

  mr(i,j,nhbp1) = j + nhbp1 - i   ! local to global no transformation

  i = nn_b(1)     ! global node no i for bc
  j = nn_b(2)     ! global node no i for bc
  sij = sqrt((x(i) - x(j))**2 + (y(i) - y(j))**2) ! side length (area)
  t_i = t_li(i)
  t_j = t_li(j)
  tb = pt5*(t_i + t_j) ! avg temp of bc edge containing nodes i and j
  coef = -a*pt5*sij
  q_dot(1) = coef*(t_inf - (pt75*t_i + pt25*t_j))**exp
  q_dot(2) = coef*(t_inf - (pt25*t_i + pt75*t_j))**exp
  res(:) = q_dot(:)

  jac_nl_flux(1,1) = -coef*exp*pt75*(t_inf - (pt75*t_i + pt25*t_j))**(exp-one)
  jac_nl_flux(1,2) = -coef*exp*pt25*(t_inf - (pt75*t_i + pt25*t_j))**(exp-one)
  jac_nl_flux(2,1) = -coef*exp*pt25*(t_inf - (pt25*t_i + pt75*t_j))**(exp-one)
  jac_nl_flux(2,2) = -coef*exp*pt75*(t_inf - (pt25*t_i + pt75*t_j))**(exp-one)

  ! merge elem bc contribution into global matrix
  call el_bc_matrix_to_global_matrix(2, nn_b, jac_nl_flux, res)

end subroutine bc_del_nl_flux
