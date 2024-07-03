subroutine vol_of_scv_quad(ne,x_g,y_g)

  ! calculate the Jacobian and its inverse for 4-node quad element and evaluate it at
  ! sub-control volume centers
  ! input:  x_g, y_g  = global (x,y) coordinates for 4-node quad
  use mod_precision
  use mod_constants
  use mod_nodes_elements
  implicit none
  save
  integer,  intent(in)                  :: ne       ! elem no
  real(dp), intent(in),  dimension(4)   :: x_g, y_g ! global element coord
  integer                               :: i, j, nscv
  real(dp)                              :: dx_dxi, dy_dxi, dx_deta, dy_deta
  real(dp),            dimension(4)     :: det_jac
  real(dp),            dimension(2,2)   :: jac
  real(dp)                              :: area, area_1, area_2
  real(dp),            dimension(5)     :: xx, yy

  ! the shape function derivatives below are defined in one long
  ! vector that is then reshaped into a (4x4) matrix where the 1st
  ! index is the local node number and the 2nd is the scv center point
  ! for convenience in verifying the entries below, they appear as if
  ! they are the transpose of the desired (4x4) matrix

  ! shape function n_i evaluated at sub-control volume center (scvc) point
  !           N_i(i,j)
  ! local node # -| |- sub-control volume center point # scv_j
!  real(dp), dimension(4,4) ::    n_i_scv(4,4) =       &
!       reshape ( (/&
!       nine/16_dp,  three/16_dp, one/16_dp,   three/16_dp, &
!       three/16_dp, nine/16_dp,  three/16_dp, one/16_dp,   &
!       one/16_dp,   three/16_dp, nine/16_dp,  three/16_dp, &
!       three/16_dp, one/16_dp,   three/16_dp, nine/16_dp/), (/4,4/) )


  !      dN_i/dxi(i,j)
  ! local node # -| |- sub-control volume center point # scv_j
  real(dp), dimension(4,4) :: dn_dxi_scv(4,4) = &
       reshape ( (/&
       -pt375,  pt375,  pt125, -pt125, &
       -pt375,  pt375,  pt125, -pt125, &
       -pt125,  pt125,  pt375, -pt375, &
       -pt125,  pt125,  pt375, -pt375/), (/4,4/) )

  !     dN_i/deta(i,j)
  ! local node # -| |- sub-control volume center point # scv_j
  real(dp), dimension(4,4) :: dn_deta_scv(4,4) = &
       reshape ( (/&
       -pt375, -pt125,  pt125,  pt375, &
       -pt125, -pt375,  pt375,  pt125, &
       -pt125, -pt375,  pt375,  pt125, &
       -pt375, -pt125,  pt125,  pt375/), (/4,4/) )
!  write(6,"(' N_i evaluated at sub-control volume center point')")
!  do i=1,4
!     write(6,"(4es13.5)") (n_i_scv(i,j), j=1,4)
!  end do
!  write(6,"(/' dN^T/dxi/dN^T/deta at sub-control volume centers')")
!  do i=1,4
!     write(6,"(/4es13.3)") (dn_dxi_scv(i,j), j=1,4)
!     write(6,"(4es13.3)")  (dn_deta_scv(i,j), j=1,4)
!  end do

  do nscv = 1,4 ! loop on # of sub-control volumes
!     write(6,"(' dn_dxi'/4es13.5)") (dn_dxi_scv(j,nscv),j=1,4)
     ! calculate the components of the Jacobian and evaluate it
     ! at the sub-control volume centers
     dx_dxi  = dot_product(x_g,dn_dxi_scv(:,nscv))
     jac(1,1) = dx_dxi
     dy_dxi  = dot_product(y_g,dn_dxi_scv(:,nscv))
     jac(1,2) = dy_dxi
     dx_deta = dot_product(x_g,dn_deta_scv(:,nscv))
     jac(2,1) = dx_deta
     dy_deta = dot_product(y_g,dn_deta_scv(:,nscv))
     jac(2,2) = dy_deta
     det_jac(nscv) = dx_dxi*dy_deta - dy_dxi*dx_deta
     v_scv(nscv,ne) = dz*det_jac(nscv)

!     write(6,"(' nscv = ',i1,' det_jac(i) = ',es13.5)") nscv, det_jac(nscv)
  end do

  ! quad area formula from byjus.com/maths/area-of-quadrilateral 
  ! or https://www.cuemath.com/measurement/area-of-quadrilateral/
  xx(1:4) = x_g(1:4)
  yy(1:4) = y_g(1:4)
  xx(5) = x_g(1)
  yy(5) = y_g(1)

  area_1 = zero; area_2 = zero
  do i=1,4
     area_1 = area_1 + xx(i)*yy(i+1)
     area_2 = area_2 + yy(i)*xx(I+1) 
  end do
  area = pt5*(area_1 - area_2)      ! area of element, not scv
  !write(6,"(' *****area error***',es13.5)") area - sum(v_scv(1:4,ne))

  end subroutine vol_of_scv_quad
