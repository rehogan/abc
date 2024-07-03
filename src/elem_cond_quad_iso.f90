
subroutine elem_cond_quad_iso(ne, x_g, y_g, t_gn, cond, dkdt, el_jac, res)
  ! calculate the Jacobian and its inverse for 4-node quad element and evaluate it at
  ! integration points (ip)
  ! input:  ne, x_g(), y_g() = elem #, global coord of elem # ne
  !         cond, dk/dT      = isotropic thermal conductivity
  ! output: Jac of Res Vec   = (4x4) element matrix
  !         res              = Res vector

  use mod_precision
  use mod_constants
  implicit none
  save

  integer,       intent(in)                  :: ne
  real(kind=dp), intent(in),  dimension(4)   :: x_g, y_g, t_gn ! global coord elem # ne
  real(kind=dp), intent(in)                  :: cond, dkdt     ! k of elem & dk/dT
  real(kind=dp), intent(out), dimension(4,4) :: el_jac          ! Jac matrix for elem e
  real(kind=dp), intent(out), dimension(4)   :: res            ! residual vector elem i
  integer                                    :: i, j, ip
  real(kind=dp)                              :: det_jac, dx_dxi, dy_dxi, dx_deta, dy_deta
  real(kind=dp),            dimension(2,2)   :: jac, jac_inv, norm_matrix
  real(kind=dp),            dimension(2,4)   :: deriv_sf_mat_ip
  real(kind=dp),            dimension(1,2)   :: a_vec
  real(kind=dp),            dimension(4,2)   :: a_mat
  real(kind=dp),            dimension(1,4)   :: el_cond_row ! elem cond row
  real(kind=dp),            dimension(4,4)   :: el_cond     ! elem cond matrix
  real(kind=dp),            dimension(4)     :: q_dot      ! heat flow vector
  real(kind=dp),            dimension(0:4,4) :: el_cond_ss ! sub-surface
  real(kind=dp)                              :: cond_sen_fac  ! 0.25*dkdt/cond
  real(kind=dp),            dimension(4,4)   :: g

  interface
     function mat_vec_mult(A, v) result (w)
       use mod_precision
       implicit none
       real(kind=dp), dimension(:,:), intent(in) :: A
       real(kind=dp), dimension(:),   intent(in) :: v
       real(kind=dp), dimension( SIZE(A,1) )     :: w
     end function mat_vec_mult
  end interface

  !      dN_i/dxi(i,j)
  ! local node # -| |- integration point # ip_j
  real(kind=dp), dimension(4,4) :: dn_dxi_ip(4,4) = &
       reshape ( (/&
       -pt375, pt375, pt125, -pt125, &
       -pt25,  pt25,  pt25,  -pt25,  &
       -pt125, pt125, pt375, -pt375, &
       -pt25,  pt25,  pt25,  -pt25/), (/4,4/) )

  !     dN_i/deta(i,j)
  ! local node # -| |- integration point # ip_j
  real(kind=dp), dimension(4,4) :: dn_deta_ip(4,4) = &
       reshape ( (/&
       -pt25,  -pt25,  pt25,  pt25,  &
       -pt125, -pt375, pt375, pt125, &
       -pt25,  -pt25,  pt25,  pt25,  &
       -pt375, -pt125, pt125, pt375/), (/4,4/) )

  ! matrix of path lengths in computational space and for each integration point ip
  real(kind=dp), dimension(4,4) :: path_len(4,2) = &
       reshape ( (/&
       one,  zero,  -one, zero, &            ! column 1
       zero, one,   zero, -one/), (/4,2/) )  ! column 2

  do ip = 1,4   ! loop over integration points ip
     dx_dxi  = dot_product(x_g,dn_dxi_ip(:,ip))
     jac(1,1) = dx_dxi
     dy_dxi  = dot_product(y_g,dn_dxi_ip(:,ip))
     jac(1,2) = dy_dxi
     dx_deta = dot_product(x_g,dn_deta_ip(:,ip))
     jac(2,1) = dx_deta
     dy_deta = dot_product(y_g,dn_deta_ip(:,ip))
     jac(2,2) = dy_deta
     det_jac = dx_dxi*dy_deta - dy_dxi*dx_deta

     jac_inv(1,1) =  dy_deta/det_jac
     jac_inv(1,2) = -dy_dxi/det_jac
     jac_inv(2,1) = -dx_deta/det_jac
     jac_inv(2,2) =  dx_dxi/det_jac

     norm_matrix(1,1) =  jac(2,2)
     norm_matrix(1,2) = -jac(2,1)
     norm_matrix(2,1) = -jac(1,2)
     norm_matrix(2,2) =  jac(1,1)
     a_mat = matmul(path_len,norm_matrix)
     a_vec(1,:) = a_mat(ip,:)

!     write(6,"(/' ip = ',i1,'  det[J] = ',es13.5)") ip,det_jac
!     write(6,"('             [J]                inverse of [J]')")
!     do i=1,2
!        write(6,"(4es13.5)") jac(i,1),jac(i,2),jac_inv(i,1),jac_inv(i,2)
!     end do

     ! put shape function derivative vectors (1,4) into (2x4) matrix
     deriv_sf_mat_ip(1,:) = dn_dxi_ip(:,ip)
     deriv_sf_mat_ip(2,:) = dn_deta_ip(:,ip)

     ! fill in the (1x2) matrix that is part of the area vector (a_vec)
     ! different for each integration point
     el_cond_row = -cond*matmul(a_vec,matmul(jac_inv,deriv_sf_mat_ip))
!     write(6,"(4es13.5)") (el_cond_row(1,j),j=1,4)
     el_cond_ss(ip,:) = el_cond_row(1,:)
  end do

!  write(6,"('   el_cond_ss')")
!  do i=1,4
!     write(6,"(4es13.5)") (el_cond_ss(i,j),j=1,4)
!  end do
  ! fill in row 0 of el_cond_ss for computational convenience
  el_cond_ss(0,:) = el_cond_ss(4,:)
  ! perform outflow - inflow for each sub-control volume
  do i=1,4
     el_cond(i,:) = el_cond_ss(i,:) - el_cond_ss(i-1,:)
  end do

  ! begin k(T) portion of residual form
  ! element heat flow vector for quad
  q_dot(1:4) = mat_vec_mult(el_cond(1:4,1:4),t_gn(1:4))

  ! calculate residual vector R and Jacobian matrix dR/dT
  cond_sen_fac = pt25*dkdt/cond                      ! coefficient
  do i=1,4
     res(i) = q_dot(i)
     do j=1,4
        el_jac(i,j) = el_cond(i,j) + cond_sen_fac*q_dot(i)
     end do
  end do

end subroutine elem_cond_quad_iso
