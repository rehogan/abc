
subroutine bc_del_spec_flux(b_el_no,nn_b,del_z,qdot_pp,dqdot_pp_dT)
  ! specified flux qdot"(T,t)
  ! energy balance statements are written as follows:
  ! outflow - inflow + storage = source
  ! in this statement, an outflow is counted as positive.
  ! in physical situations, engineers typically want a positive flux
  ! to increase the temperature, which is backwards from the overall
  ! conservation of energy statement. The equations that follow will 
  ! have some negative signs to allow for positive applied flux to
  ! produce an increase in temperature.

  ! b_el_no = boundary elem no
  ! nn_b()  = node no that define side length S_(i,j) for this b_el_no
  ! del_z   = depth of element in 3rd dimension, likely unity
  ! qdot_pp = qdot", W/m^2 
  ! dqdot_pp_dT = d(qdot")/dT, derivative of specified flux wrt T

  use mod_precision
  use mod_constants
  use mod_nodes_elements
  implicit none
  integer,       intent(in)               :: b_el_no
  integer,       intent(in), dimension(2) :: nn_b
  real(kind=dp), intent(in)               :: del_z, qdot_pp, dqdot_pp_dT
  real(kind=dp)                 :: tb, t_i, t_j, sij, coef
  real(kind=dp), dimension(2,2) :: jac_spec_flux
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
  coef = pt5*sij*del_z
  q_dot(1) = coef*qdot_pp ! heat flow, W
  q_dot(2) = q_dot(1)
  res(:) = -q_dot(:)  ! note the negative sign

  jac_spec_flux(1,1) = -coef*pt5*dqdot_pp_dT ! note the negative sign
  jac_spec_flux(1,2) = jac_spec_flux(1,1)
  jac_spec_flux(2,1) = jac_spec_flux(1,1)
  jac_spec_flux(2,2) = jac_spec_flux(1,1)

  ! merge elem bc contribution into global matrix
  call el_bc_matrix_to_global_matrix(2, nn_b, jac_spec_flux, res)

end subroutine bc_del_spec_flux
