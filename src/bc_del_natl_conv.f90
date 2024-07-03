
subroutine bc_del_natl_conv(geom,b_el_no,nn_b,t_inf)

  ! integrals resulting from integration of radiation along element   
  ! boundary. nomenclature is similar to that in JTHT paper
  ! B. F. Blackwell and R. E. Hogan, Numerical Solution of 
  ! Axisymmetric Heat Conduction Problems Using Finite control Volume
  ! Technique, JTHT, Vol. 7, No. 3, July-Sept. 1993

  ! b_el_no = boundary elem no
  ! nn_b()  = node no that define side length S_(i,j) for this b_el_no
  ! T_inf   = temp to which surf convects to

  use mod_precision
  use mod_constants
  use mod_bc
  use mod_nodes_elements
  implicit none
  character(len=7),       intent(in) :: geom
  integer,  intent(in)               :: b_el_no
  integer,  intent(in), dimension(2) :: nn_b
  real(dp), intent(in)               :: t_inf
  real(dp)                 :: tbar, t_i, t_j, sij, coef
  real(dp), dimension(2,2) :: jac_conv
  real(dp), dimension(2)   :: q_dot, res
  real(dp)                 :: qdot_pp,  h_conv, dh_conv_dt_w, dh_conv_dt_inf
  real(dp)                 :: h, dhdt

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

     subroutine flux_natl_conv(geom, t_w, t_inf, qdot_pp, &
          h_conv, dh_conv_dt_w, dh_conv_dt_inf)
       use mod_precision
       use mod_constants
       use mod_bc
       use mod_air_properties
       use mod_lin_interp
       implicit none
       save
       character(len=7), intent(in) :: geom
       real(dp),         intent(in) :: t_w, t_inf
       real(dp), intent(out)        :: qdot_pp,  h_conv, & 
                                       dh_conv_dt_w, dh_conv_dt_inf
     end subroutine flux_natl_conv

    end interface

  mr(i,j,nhbp1) = j + nhbp1 - i  ! local to global no transformation

  i = nn_b(1)
  j = nn_b(2)
  sij = sqrt((x(i) - x(j))**2 + (y(i) - y(j))**2) ! side length (area)
  t_i = t_li(i)
  t_j = t_li(j)
  tbar = pt5*(t_i + t_j) ! avg temp of bc edge containing nodes i and j

  call flux_natl_conv(geom, tbar, t_inf, qdot_pp, &
     h_conv, dh_conv_dt_w, dh_conv_dt_inf)
  coef = h_conv*pt5*sij
  q_dot(1) = coef*(pt75*(t_i - t_inf) + pt25*(t_j - t_inf))
  q_dot(2) = coef*(pt25*(t_i - t_inf) + pt75*(t_j - t_inf))
  res(:) = q_dot(:)

  h = h_conv
  dhdt = dh_conv_dt_w
  jac_conv(1,1) = coef*pt75 + pt5*(dhdt/h)*q_dot(1)
  jac_conv(1,2) = coef*pt25 + pt5*(dhdt/h)*q_dot(1)
  jac_conv(2,1) = coef*pt25 + pt5*(dhdt/h)*q_dot(2)
  jac_conv(2,2) = coef*pt75 + pt5*(dhdt/h)*q_dot(2)

  ! merge elem bc contribution into global matrix
  call el_bc_matrix_to_global_matrix(2, nn_b, jac_conv, res)

end subroutine bc_del_natl_conv
