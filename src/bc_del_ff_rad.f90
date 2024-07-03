subroutine bc_del_ff_rad(b_el_no,nn_b,eps,fij,t_r)
  ! far field radiation bc in delta form
  ! integrals resulting from integration of radiation along element   
  ! boundary. nomenclature is similar to that in JTHT paper
  ! B. F. Blackwell and R. E. Hogan, Numerical Solution of 
  ! Axisymmetric Heat Conduction Problems Using Finite control Volume
  ! Technique, JTHT, Vol. 7, No. 3, July-Sept. 1993

  ! b_el_no = boundary elem no
  ! i,j     = node no that define side length S_(i,j) for this b_el_no
  ! T_r     = temp to which surf radiates to

  use mod_precision
  use mod_constants
  use mod_nodes_elements
  implicit none
  integer,       intent(in)               :: b_el_no
  integer,       intent(in), dimension(2) :: nn_b
  real(kind=dp), intent(in)               :: eps, fij, t_r
  real(kind=dp)                 :: avg_t4, tb, t_i, t_j, sij, coef, tt_i4, tt_j4
  real(kind=dp), dimension(2,2) :: jac_ffr, ptt4pt  
  real(kind=dp), dimension(2)   :: q_dot, res
  integer                       :: i, j, ii, jj, kk, mr
  interface
     function avg_t4(tb,t_i)
       use mod_precision
       use mod_constants
       implicit none
       real(kind=dp), intent(in) :: tb, t_i
     end function avg_t4

     subroutine el_bc_matrix_to_global_matrix(nnpe, nn_b, el_jac, res)
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
  j = nn_b(2)     ! global node no j for bc
  sij = sqrt((x(i) - x(j))**2 + (y(i) - y(j))**2) ! side length (area)
  t_i = t_li(i)
  t_j = t_li(j)
  tb = pt5*(t_i + t_j) ! avg temp of bc edge containing nodes i and j
  tt_i4 = avg_t4(tb,t_i)
  tt_j4 = avg_t4(tb,t_j)
  coef = eps*sigma_si*pt5*sij*fij
  q_dot(1) = coef*(tt_i4 - t_r**4)
  q_dot(2) = coef*(tt_j4 - t_r**4)
  res(:) = q_dot(:)

  ! 
  ! partial T-~-i**4(Tbar,Ti)/partial Ti
  ! ptt4pt(1,1) = partial (T^~_i)**4/partial Ti
  ptt4pt(1,1) = 0.2_dp*(3._dp*tb**3 + 3.5_dp*tb**2*t_i + 4._dp*tb*t_i**2 + &
       4.5_dp*t_i**3)
  ! ptt4pt(2,2) = partial (T^~_j)**4/partial Tj
  ptt4pt(2,2) = 0.2_dp*(3._dp*tb**3 + 3.5_dp*tb**2*t_j + 4._dp*tb*t_j**2 + &
       4.5*t_j**3)
  ! ptt4pt(1,2) = partial (T^~_i)**4/partial Tj
  ptt4pt(1,2) = 0.2_dp*(2._dp*tb**3 + 1.5_dp*tb**2*t_i + tb*t_i**2 + &
       0.5_dp*t_i**3)
  ! ptt4pt(2,1) = partial (T^~_j)**4/partial Ti
  ptt4pt(2,1) = 0.2_dp*(2._dp*tb**3 + 1.5_dp*tb**2*t_j + tb*t_j**2 + &
       0.5_dp*t_j**3)

  jac_ffr(:,:) = coef*ptt4pt(:,:)


  ! merge elem bc contribution into global matrix
  call el_bc_matrix_to_global_matrix(2, nn_b, jac_ffr, res)

end subroutine bc_del_ff_rad

function avg_t4(tb,t_i)
  use mod_precision
  use mod_constants
  implicit none
  real(kind=dp), intent(in) :: tb, t_i
  real(kind=dp) :: avg_t4
  avg_t4 = 0.2_dp*(tb**4 + tb**3*t_i + tb**2*t_i**2 + &
       tb*t_i**3 + t_i**4)
end function avg_t4
